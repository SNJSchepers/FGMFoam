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

"""Check a built case *before* the solver is launched.

A parallel FGMFoam run can occupy a machine for hours. Several ways of wasting
that time are detectable in milliseconds beforehand, and all of them are
deterministic — no model is involved in any judgement here.

What is checked:

* **Conditions inside the FGM table.** Every thermochemical quantity comes from
  the table, so an initial enthalpy or progress variable outside its tabulated
  range means `lookupFGM_2D` cannot return sensible values. The operating
  pressure must match the pressure the table was generated at.
* **The time step against the mesh.** The solver runs at a fixed step
  (`adjustTimeStep no`, and `maxCo` in controlDict is inert), so nothing
  throttles it at run time: too large a step does not self-correct, it diverges.
  A Courant estimate from the real mesh spacing turns the temperature rule of
  thumb into a number.
* **Bookkeeping that silently wastes a run.** A write interval smaller than the
  time step, a decomposition with almost no cells per rank, an implausible step
  count.

The mesh-dependent checks need `constant/polyMesh`, so the gate runs after
`blockMesh` (seconds) and before the solver (hours). Without a mesh they are
reported as skipped rather than guessed: the shipped `blockMeshDict`s use named
vertices, variable expansion and `#eval{}`, and a half-correct reimplementation
of that would produce a *misleading* Courant number, which is worse than none.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import fgm_table
import of_mesh
import of_writer

# A Courant number above this is flagged. The solver is implicit, so this is an
# accuracy/robustness guide rather than a hard stability limit.
COURANT_WARN: float = 1.0

# Below this many cells per MPI rank, communication tends to dominate.
MIN_CELLS_PER_RANK: int = 1000

# Flow-through times after which an inlet-driven steady flame is normally
# converged. Below one pass the result still reflects the initial condition.
RECOMMENDED_FLOW_THROUGHS: float = 3.0

# A step count beyond this is almost always a units slip (ms vs s).
MAX_PLAUSIBLE_STEPS: int = 10_000_000


@dataclass(frozen=True)
class Finding:
    """One check result. ``error`` blocks a run; ``warning`` and ``info`` do not."""

    severity: str  # "error" | "warning" | "info"
    code: str
    message: str
    details: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        out = {"severity": self.severity, "code": self.code, "message": self.message}
        if self.details:
            out["details"] = self.details
        return out


# --------------------------------------------------------------------------- #
# Minimal OpenFOAM readers (the agent only ever needs uniform fields)
# --------------------------------------------------------------------------- #

_SCALAR = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"


def _uniform_scalar(path: Path) -> Optional[float]:
    """``internalField uniform <x>;`` as a float, or None."""
    try:
        text = path.read_text()
    except OSError:
        return None
    match = re.search(rf"internalField\s+uniform\s+({_SCALAR})\s*;", text)
    return float(match.group(1)) if match else None


def _max_velocity_magnitude(path: Path) -> Optional[float]:
    """Largest |U| among the uniform vectors in a velocity field file.

    Covers the internal field and any patch ``value``/``inletValue`` entries, so
    a prescribed inlet is seen even when the interior starts at rest.
    """
    try:
        text = path.read_text()
    except OSError:
        return None
    vectors = re.findall(
        rf"uniform\s*\(\s*({_SCALAR})\s+({_SCALAR})\s+({_SCALAR})\s*\)", text
    )
    if not vectors:
        return None
    return max(math.sqrt(sum(float(c) ** 2 for c in v)) for v in vectors)


def _set_fields_values(case_dir: Path) -> list[float]:
    """Every ``volScalarFieldValue Yc <x>`` in setFieldsDict."""
    try:
        text = (case_dir / "system" / "setFieldsDict").read_text()
    except OSError:
        return []
    return [
        float(m) for m in re.findall(
            rf"volScalarFieldValue\s+Yc\s+({_SCALAR})", text
        )
    ]


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #


def _check_table(case_dir: Path, findings: list[Finding]) -> None:
    table_path = case_dir / "constant" / "lookUp" / "database.fgm"
    if not table_path.is_file():
        table_path = of_writer.TABLE_SRC
    info = fgm_table.read_table_safe(table_path)
    if info is None:
        findings.append(Finding(
            "warning", "table_unreadable",
            f"Could not read the FGM table at {table_path}; "
            "the range checks were skipped.",
        ))
        return

    pressure = _uniform_scalar(case_dir / "0" / "p")
    if pressure is not None and not math.isclose(
        pressure, info.pressure, rel_tol=1e-3
    ):
        findings.append(Finding(
            "error", "pressure_mismatch",
            f"The case operates at {pressure:g} Pa but the FGM table was "
            f"generated at {info.pressure:g} Pa. Every thermochemical quantity "
            "comes from that table, so the run would be inconsistent. Rebuild "
            "the case at the table's pressure, or regenerate the table.",
            {"case_pressure": pressure, "table_pressure": info.pressure},
        ))

    ht = _uniform_scalar(case_dir / "0" / "ht")
    if ht is not None and info.dimension >= 2:
        lo, hi = info.range_of(1)
        if not lo <= ht <= hi:
            findings.append(Finding(
                "error", "enthalpy_outside_table",
                f"The initial enthalpy {ht:.4g} J/kg lies outside the table's "
                f"tabulated range [{lo:.4g}, {hi:.4g}] J/kg. Lookups cannot "
                "return meaningful values there — the conditions need a table "
                "that covers them.",
                {"ht": ht, "table_min": lo, "table_max": hi},
            ))
        else:
            margin = min(ht - lo, hi - ht) / (hi - lo)
            if margin < 0.02:
                findings.append(Finding(
                    "warning", "enthalpy_near_table_edge",
                    f"The initial enthalpy {ht:.4g} J/kg is within "
                    f"{margin * 100:.1f}% of the edge of the table's range "
                    f"[{lo:.4g}, {hi:.4g}]. Heat loss or gain during the run "
                    "may push it outside.",
                    {"ht": ht, "table_min": lo, "table_max": hi},
                ))

    lo, hi = info.range_of(0)
    values = [v for v in (_uniform_scalar(case_dir / "0" / "Yc"),) if v is not None]
    values += _set_fields_values(case_dir)
    outside = [v for v in values if not lo - abs(lo) * 1e-6 <= v <= hi * (1 + 1e-6)]
    if outside:
        findings.append(Finding(
            "error", "progress_variable_outside_table",
            f"Progress-variable value(s) {', '.join(f'{v:.4g}' for v in outside)} "
            f"lie outside the table's range [{lo:.4g}, {hi:.4g}] kmol/kg. The "
            "burnt value must come from the same equilibrium the table was "
            "built at.",
            {"values": outside, "table_min": lo, "table_max": hi},
        ))


def _check_run_controls(case_dir: Path, findings: list[Finding]) -> None:
    entries = of_writer.read_control_entries(
        case_dir, ["endTime", "deltaT", "writeInterval"]
    )
    try:
        end_time = float(entries["endTime"])
        delta_t = float(entries["deltaT"])
        write_interval = float(entries["writeInterval"])
    except (KeyError, TypeError, ValueError):
        findings.append(Finding(
            "warning", "control_dict_unreadable",
            "Could not read endTime/deltaT/writeInterval from system/controlDict.",
        ))
        return

    if delta_t <= 0 or end_time <= 0:
        findings.append(Finding(
            "error", "nonpositive_run_controls",
            f"endTime={end_time:g} s and deltaT={delta_t:g} s must both be positive.",
        ))
        return

    steps = end_time / delta_t
    if steps > MAX_PLAUSIBLE_STEPS:
        findings.append(Finding(
            "warning", "implausible_step_count",
            f"endTime/deltaT is {steps:,.0f} time steps. That is almost always a "
            "units slip (milliseconds given as seconds) rather than an intended run.",
            {"n_timesteps": steps},
        ))
    if write_interval < delta_t:
        findings.append(Finding(
            "warning", "write_interval_below_timestep",
            f"writeInterval ({write_interval:g} s) is smaller than deltaT "
            f"({delta_t:g} s), so output is written every step at best. Reduce "
            "sim_write_steps, or raise the number of snapshots deliberately.",
            {"write_interval": write_interval, "delta_t": delta_t},
        ))


def _check_mesh(case_dir: Path, findings: list[Finding]) -> None:
    cells = of_mesh.cell_count(case_dir)
    spacing = of_mesh.min_spacing(case_dir)
    if cells is None or spacing is None:
        findings.append(Finding(
            "info", "mesh_not_built",
            "No mesh found yet (constant/polyMesh), so the Courant and "
            "decomposition checks were skipped. They run once blockMesh has.",
        ))
        return

    findings.append(Finding(
        "info", "mesh",
        f"Mesh: {cells:,} cells, smallest spacing {spacing:.3g} m.",
        {"n_cells": cells, "min_spacing": spacing},
    ))

    n_procs = _n_procs(case_dir)
    if n_procs and n_procs > cells:
        findings.append(Finding(
            "error", "more_ranks_than_cells",
            f"The case is decomposed over {n_procs} ranks but the mesh has only "
            f"{cells:,} cells.",
            {"n_procs": n_procs, "n_cells": cells},
        ))
    elif n_procs and cells / n_procs < MIN_CELLS_PER_RANK:
        findings.append(Finding(
            "warning", "thin_decomposition",
            f"{cells:,} cells over {n_procs} ranks is {cells / n_procs:,.0f} "
            "cells per rank; below roughly "
            f"{MIN_CELLS_PER_RANK:,}, MPI communication tends to cost more than "
            "the extra ranks gain. Fewer processes may run faster.",
            {"n_procs": n_procs, "cells_per_rank": cells / n_procs},
        ))

    entries = of_writer.read_control_entries(case_dir, ["deltaT"])
    try:
        delta_t = float(entries["deltaT"])
    except (KeyError, TypeError, ValueError):
        return
    velocity = _max_velocity_magnitude(case_dir / "0" / "U")
    if velocity is None:
        return
    if velocity <= 0.0:
        findings.append(Finding(
            "info", "courant_not_determinable",
            "The initial velocity field is at rest, so the convective Courant "
            "number cannot be bounded before the run: it is set by the "
            "propagating flame, not by the initial condition. Watch the Courant "
            "number reported in log.FGMFoam once the run starts.",
            {"delta_t": delta_t, "min_spacing": spacing},
        ))
        return

    courant = velocity * delta_t / spacing
    detail = {
        "courant": courant, "velocity": velocity,
        "delta_t": delta_t, "min_spacing": spacing,
    }
    if courant > COURANT_WARN:
        suggested = COURANT_WARN * spacing / velocity
        detail["suggested_timestep"] = suggested
        findings.append(Finding(
            "warning", "courant_high",
            f"Estimated Courant number {courant:.2f} (|U|={velocity:g} m/s, "
            f"dt={delta_t:g} s, smallest cell {spacing:.3g} m). The solver runs "
            "at a fixed step and will not reduce it, so consider a time step "
            f"around {suggested:.2g} s. Gas expansion across the flame raises "
            "the velocity further, so treat this as a lower bound.",
            detail,
        ))
    else:
        findings.append(Finding(
            "info", "courant_ok",
            f"Estimated Courant number {courant:.2f} from the initial velocity "
            "field. Gas expansion across the flame will raise it during the run.",
            detail,
        ))


def _n_procs(case_dir: Path) -> Optional[int]:
    """``numberOfSubdomains`` from the case's decomposeParDict."""
    try:
        text = (case_dir / "system" / "decomposeParDict").read_text()
    except OSError:
        return None
    match = re.search(r"numberOfSubdomains\s+(\d+)\s*;", text)
    return int(match.group(1)) if match else None


def _check_residence_time(case_dir: Path, findings: list[Finding]) -> None:
    """Is the run long enough for the flow to have crossed the domain?

    A steady, inlet-driven flame is only a result once the fresh gas has swept
    the domain and the initial condition has been flushed out. The flow-through
    time is the axial domain length over the inlet velocity — conservative,
    since the burnt gas downstream moves faster. Not applicable where nothing
    flows in (the flamelet): there the flame itself sets the time scale.
    """
    box = of_mesh.bounds(case_dir)
    velocity = _max_velocity_magnitude(case_dir / "0" / "U")
    entries = of_writer.read_control_entries(case_dir, ["endTime"])
    try:
        end_time = float(entries["endTime"])
    except (KeyError, TypeError, ValueError):
        return
    if box is None or velocity is None:
        return
    if velocity <= 0.0:
        findings.append(Finding(
            "info", "residence_time_not_applicable",
            "No inflow, so there is no flow-through time: the flame propagates "
            "into gas at rest and sets its own time scale.",
        ))
        return

    length = box[0][1] - box[0][0]
    flow_through = length / velocity
    passes = end_time / flow_through
    details = {
        "flow_through_time": flow_through, "end_time": end_time,
        "flow_throughs": passes, "domain_length": length, "inlet_velocity": velocity,
        "recommended_end_time": RECOMMENDED_FLOW_THROUGHS * flow_through,
    }
    if passes < 1.0:
        findings.append(Finding(
            "warning", "shorter_than_one_flow_through",
            f"endTime ({end_time:g} s) is {passes:.2f} flow-through times: the "
            f"domain is {length * 1e3:.1f} mm long and the inlet gas moves at "
            f"{velocity:g} m/s, so one pass takes {flow_through:.3g} s. The fresh "
            "gas has not yet crossed the domain, so the result still reflects the "
            "initial condition. For a converged steady flame run at least "
            f"{RECOMMENDED_FLOW_THROUGHS:g} flow-throughs "
            f"(~{RECOMMENDED_FLOW_THROUGHS * flow_through:.2g} s).",
            details,
        ))
    elif passes < RECOMMENDED_FLOW_THROUGHS:
        findings.append(Finding(
            "info", "few_flow_throughs",
            f"endTime covers {passes:.1f} flow-through times ({flow_through:.3g} s "
            f"each). A steady flame is usually converged after about "
            f"{RECOMMENDED_FLOW_THROUGHS:g}; check that the flame length has "
            "stopped changing between the last snapshots.",
            details,
        ))
    else:
        findings.append(Finding(
            "info", "residence_time_ok",
            f"endTime covers {passes:.1f} flow-through times ({flow_through:.3g} s each).",
            details,
        ))


def check_case(case_dir: Path) -> list[Finding]:
    """Run every pre-run check against a built case."""
    case_dir = Path(case_dir)
    findings: list[Finding] = []
    if not case_dir.is_dir():
        return [Finding("error", "case_missing", f"No such case: {case_dir}")]
    _check_table(case_dir, findings)
    _check_run_controls(case_dir, findings)
    _check_mesh(case_dir, findings)
    _check_residence_time(case_dir, findings)
    return findings


def summarise(findings: list[Finding]) -> dict[str, Any]:
    """Group findings by severity into a result the model can relay."""
    errors = [f for f in findings if f.severity == "error"]
    warnings = [f for f in findings if f.severity == "warning"]
    status = "errors" if errors else ("warnings" if warnings else "ok")
    return {
        "status": status,
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "findings": [f.as_dict() for f in findings],
    }
