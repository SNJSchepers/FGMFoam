# FGMFoam agent — a code-aware CLI assistant for the FGMFoam solver.
#
# Copyright (C) 2025-2026 Stijn Schepers, Eindhoven University of Technology
#
# This file is part of FGMFoam. FGMFoam is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# FGMFoam is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# FGMFoam. If not, see <https://www.gnu.org/licenses/>.

"""Campaigns: one variable swept over several values, run and plotted together.

A campaign is the unit of actual research use — "flame speed against unburnt
temperature", "flame length against plate temperature" — rather than a single
case. It builds one case per value, runs them under a CPU plan the **user** has
chosen, analyses each with :mod:`postrun`, and collects the results into a table
and a graph.

Two rules are enforced in code rather than left to the model:

* **Only conditions can be swept.** The mesh and boundary-condition types come
  from the template unchanged, so geometry (hole size, plate length) is not a
  campaign variable; and the FGM table fixes the equivalence ratio and pressure,
  so those are refused with the reason rather than silently accepted.
* **The CPU plan is never guessed.** How to spend the machine — one case at a
  time on all cores, or several cases side by side on fewer each — changes
  wall-clock time by a large factor and depends on what else the machine is
  doing. :func:`cpu_options` enumerates the sensible splits for the detected
  core count and the tool refuses to start until one is chosen.
"""

from __future__ import annotations

import csv
import json
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import case_templates
import of_writer
import postrun
import validator
from config import AGENT_DIR

CAMPAIGN_DIR: Path = AGENT_DIR / "campaigns"

# What may be swept, what rebuilds the mesh, what is plotted and how axes are
# labelled — all declared by each case plugin (agent/cases/<name>/case.py).
SWEEPABLE: dict[str, tuple[str, ...]] = {
    t.name: t.sweepable for t in case_templates.TEMPLATES
}
GEOMETRIC: frozenset[str] = frozenset(
    name for t in case_templates.TEMPLATES for name in t.geometric
)

# Parameters a user may reasonably ask to sweep that cannot be swept.
REFUSED: dict[str, str] = {
    "equivalence_ratio": (
        "the FGM table is generated at a single equivalence ratio (phi = 0.7). "
        "Sweeping it needs one table per value, regenerated with Chem1D or "
        "another 1D flame solver"
    ),
    "p_0": (
        "the FGM table is generated at a single pressure (101325 Pa). Sweeping "
        "pressure needs a table regenerated at each value"
    ),
}

# Words that indicate a geometric variable; the agent never edits the mesh.
GEOMETRY_WORDS: tuple[str, ...] = (
    "hole", "slot", "width", "height", "length", "thickness", "gap", "spacing",
    "mesh", "geometry", "domain", "plate_size", "diameter",
)

METRIC: dict[str, dict[str, str]] = {
    t.name: t.metric for t in case_templates.TEMPLATES if t.metric
}
AXIS_LABEL: dict[str, str] = {
    name: label for t in case_templates.TEMPLATES for name, label in t.axis_labels.items()
}


@dataclass
class Campaign:
    """One swept variable and the cases that realise it."""

    name: str
    case_type: str
    variable: str
    values: list[float]
    base_params: dict[str, Any] = field(default_factory=dict)
    n_procs_per_case: int = 0
    n_concurrent_cases: int = 0
    cases: list[dict[str, Any]] = field(default_factory=list)
    created_at: str = ""
    # "prepared" (built, not started) or "launched".
    status: str = "prepared"
    # Run length in flow-through times, if sim_time was derived per case.
    flow_throughs: Optional[float] = None

    @property
    def directory(self) -> Path:
        return CAMPAIGN_DIR / self.name

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name, "case_type": self.case_type,
            "variable": self.variable, "values": self.values,
            "base_params": self.base_params,
            "n_procs_per_case": self.n_procs_per_case,
            "n_concurrent_cases": self.n_concurrent_cases,
            "cases": self.cases, "created_at": self.created_at,
            "status": self.status, "flow_throughs": self.flow_throughs,
        }


def _slug(value: float) -> str:
    """Filesystem-safe token for a swept value ("0.25" -> "0p25")."""
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def safe_name(text: str) -> str:
    """A campaign name usable as a path segment."""
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "", str(text).strip())
    return cleaned or "campaign"


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def check_variable(case_type: str, variable: str) -> Optional[str]:
    """Reason ``variable`` cannot be swept for ``case_type``, or None."""
    name = (variable or "").strip()
    if not name:
        return "no variable given to sweep"
    if name in REFUSED:
        return f"{name} cannot be swept: {REFUSED[name]}"
    allowed = SWEEPABLE.get(case_type, ())
    template_geometry = [
        g for g in allowed
        if g in case_templates.TEMPLATE_BY_NAME[case_type].geometric
    ] if case_type in case_templates.TEMPLATE_BY_NAME else []
    if name in allowed:
        return None
    lowered = name.lower()
    if any(word in lowered for word in GEOMETRY_WORDS):
        return (
            f"{name} is not a geometric parameter this template exposes. "
            + (
                f"{case_type}'s mesh is parameterised by "
                + ", ".join(template_geometry)
                + " — use one of those"
                if template_geometry else
                f"{case_type}'s blockMeshDict is not parameterised, so its "
                "geometry cannot be swept; edit the template's blockMeshDict "
                "by hand for that"
            )
        )
    return (
        f"{name} is not a sweepable parameter of {case_type}. Available: "
        + ", ".join(allowed)
    )


def cpu_options(n_cases: int, cores: int) -> list[dict[str, Any]]:
    """Sensible ways to spend ``cores`` on ``n_cases`` cases.

    Every option keeps ``procs_per_case * concurrent_cases`` within the physical
    core count: oversubscribing physical cores makes a CFD run slower, not
    faster. Ordered from "one case at a time, all cores each" to "many cases at
    once, one core each".
    """
    options: list[dict[str, Any]] = []
    procs = cores
    while procs >= 1:
        concurrent = min(n_cases, max(1, cores // procs))
        entry = {
            "n_procs_per_case": procs,
            "n_concurrent_cases": concurrent,
            "cores_used": procs * concurrent,
            "description": (
                f"{procs} core{'s' if procs > 1 else ''} per case, "
                f"{concurrent} case{'s' if concurrent > 1 else ''} at a time"
                + (" (fully sequential)" if concurrent == 1 else "")
                + (" (each case serial)" if procs == 1 else "")
            ),
            "waves": math.ceil(n_cases / concurrent),
        }
        if not any(
            o["n_procs_per_case"] == procs and o["n_concurrent_cases"] == concurrent
            for o in options
        ):
            options.append(entry)
        procs //= 2
    return options


def validate_cpu_plan(
    n_procs_per_case: int, n_concurrent_cases: int, cores: int,
) -> Optional[str]:
    """Reason a CPU plan is unusable, or None."""
    if n_procs_per_case < 1 or n_concurrent_cases < 1:
        return "n_procs_per_case and n_concurrent_cases must both be at least 1"
    used = n_procs_per_case * n_concurrent_cases
    if used > cores:
        return (
            f"{n_procs_per_case} cores per case x {n_concurrent_cases} cases = "
            f"{used} cores, but only {cores} physical cores are available. "
            "Oversubscribing physical cores makes CFD runs slower, not faster — "
            "choose a split that fits, or say explicitly to oversubscribe"
        )
    return None


# --------------------------------------------------------------------------- #
# Planning and persistence
# --------------------------------------------------------------------------- #


def plan(
    name: str, case_type: str, variable: str, values: list[float],
    base_params: dict[str, Any], n_procs_per_case: int, n_concurrent_cases: int,
    flow_throughs: Optional[float] = None,
) -> Campaign:
    """Turn a sweep into one case specification per value.

    With ``flow_throughs``, each case's ``sim_time`` is that many flow-through
    times of *its own* conditions (domain length / inlet velocity). A single
    shared sim_time would over-run the fast cases and under-run the slow ones
    when inlet velocity or plate thickness is the swept variable.
    """
    template = case_templates.find_template(case_type)
    campaign = Campaign(
        name=safe_name(name), case_type=template.name, variable=variable,
        values=list(values), base_params=dict(base_params),
        n_procs_per_case=n_procs_per_case,
        n_concurrent_cases=n_concurrent_cases,
        created_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        flow_throughs=flow_throughs,
    )
    for value in campaign.values:
        params = dict(base_params)
        params[variable] = value
        params["n_procs"] = n_procs_per_case
        if flow_throughs:
            filled = validator.apply_defaults(dict(params), template)
            length = case_templates.domain_length(template.name, filled)
            velocity = filled.get("inlet_velocity")
            if not length or not velocity:
                raise ValueError(
                    f"flow_throughs needs an inlet-driven case with a known domain "
                    f"length; {template.name} has none — give sim_time instead"
                )
            params["sim_time"] = float(f"{flow_throughs * length / float(velocity):.3g}")
        case_name = f"{campaign.name}_{variable}_{_slug(value)}"
        campaign.cases.append({
            "value": value,
            "case_name": f"{template.prefix}{case_name}",
            "params": params,
        })
    return campaign


def save(campaign: Campaign) -> Path:
    campaign.directory.mkdir(parents=True, exist_ok=True)
    path = campaign.directory / "campaign.json"
    path.write_text(json.dumps(campaign.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path


def load(name: str) -> Optional[Campaign]:
    path = CAMPAIGN_DIR / safe_name(name) / "campaign.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return Campaign(**data)


def list_campaigns() -> list[str]:
    if not CAMPAIGN_DIR.is_dir():
        return []
    return sorted(
        d.name for d in CAMPAIGN_DIR.iterdir()
        if (d / "campaign.json").is_file()
    )


# --------------------------------------------------------------------------- #
# Results
# --------------------------------------------------------------------------- #


def collect(campaign: Campaign, threshold: float = postrun.FLAME_TEMPERATURE) -> dict[str, Any]:
    """Analyse every case in the campaign and assemble the results table."""
    metric = METRIC.get(campaign.case_type, {})
    key = metric.get("key", "")
    rows: list[dict[str, Any]] = []

    for entry in campaign.cases:
        case_dir = of_writer.case_path(entry["case_name"])
        row: dict[str, Any] = {
            campaign.variable: entry["value"],
            "case_name": entry["case_name"],
            "status": "not_built",
            key: None,
        }
        if case_dir.is_dir():
            analysis = postrun.analyse(case_dir, threshold=threshold)
            run = analysis.get("run", {})
            flame = analysis.get("flame", {})
            row["status"] = flame.get("status", analysis.get("status", "unknown"))
            row["time_reached"] = run.get("time_reached")
            # The settled value is what says whether the step was right; the
            # peak is the ignition transient and is kept only for reference.
            row["courant_max_settled"] = run.get("courant_max_settled")
            row["courant_peak"] = run.get("courant_max")
            if key and key in flame:
                row[key] = flame[key]
            # Keep the qualitative outcomes visible in the table: a flashback is
            # a result, not a missing value.
            if flame.get("status") in ("flashback", "blow_off", "longer_than_domain"):
                row["note"] = flame.get("message", "")[:160]
        rows.append(row)

    rows.sort(key=lambda r: r[campaign.variable])
    measured = [r for r in rows if r.get(key) is not None]
    return {
        "campaign": campaign.name,
        "case_type": campaign.case_type,
        "variable": campaign.variable,
        "metric": key,
        "n_cases": len(rows),
        "n_measured": len(measured),
        "rows": rows,
    }


def write_csv(campaign: Campaign, results: dict[str, Any]) -> Path:
    """Write the results table next to the campaign manifest."""
    campaign.directory.mkdir(parents=True, exist_ok=True)
    path = campaign.directory / "results.csv"
    columns: list[str] = []
    for row in results["rows"]:
        for name in row:
            if name not in columns:
                columns.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in results["rows"]:
            writer.writerow(row)
    return path


def plot(campaign: Campaign, results: dict[str, Any]) -> tuple[Optional[Path], Optional[str]]:
    """Graph the measurement against the swept variable.

    Returns ``(path, error)``. Plotting is optional: without matplotlib the CSV
    is still written, and the caller says so rather than failing the campaign.
    """
    metric = METRIC.get(campaign.case_type)
    if not metric:
        return None, f"no plot is defined for case type {campaign.case_type!r}"

    points = [
        (row[campaign.variable], row[metric["key"]])
        for row in results["rows"] if row.get(metric["key"]) is not None
    ]
    if len(points) < 2:
        return None, (
            f"only {len(points)} case(s) produced a measurement, which is not "
            "enough to plot. Check the per-case status in results.csv"
        )

    try:
        import matplotlib
        matplotlib.use("Agg")           # no display on a compute node
        import matplotlib.pyplot as plt
    except ImportError:
        return None, (
            "matplotlib is not installed, so only the CSV was written: "
            "agent/.venv/bin/pip install matplotlib"
        )

    # Cases that ran but produced no number are still information — a flashback
    # at the low-velocity end is the physics, not a gap in the data.
    failures = [
        (row[campaign.variable], row["status"])
        for row in results["rows"]
        if row.get(metric["key"]) is None
        and row.get("status") in ("flashback", "blow_off", "longer_than_domain")
    ]

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    figure, axes = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    axes.plot(xs, ys, marker="o", linewidth=1.6, color="#b4462a")
    for x_failed, status in failures:
        axes.axvline(x_failed, color="#7a7a7a", linestyle=":", linewidth=1.2)
        axes.annotate(
            status.replace("_", " "), xy=(x_failed, min(ys)),
            rotation=90, fontsize=7, color="#555555",
            ha="right", va="bottom",
        )
    axes.set_xlabel(AXIS_LABEL.get(campaign.variable, campaign.variable))
    axes.set_ylabel(metric["label"])
    axes.set_title(f"{metric['title']} - {campaign.case_type}")
    axes.grid(True, alpha=0.25, linewidth=0.6)
    figure.tight_layout()

    campaign.directory.mkdir(parents=True, exist_ok=True)
    path = campaign.directory / f"{campaign.name}.png"
    figure.savefig(path)
    plt.close(figure)
    return path, None
