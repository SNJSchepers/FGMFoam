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

"""The case templates the agent can set up, and what each one needs.

One entry per template case under ``testcases/``. A template declares three
things, which is everything the rest of the agent needs to treat it generically:

* **which parameters apply** — the 1D flamelet has no inlet and no plate, so
  asking for an inlet velocity or a plate temperature there would be nonsense;
* **which placeholder tokens its dictionaries carry** — the two cases spell the
  unburnt temperature differently (`INLETTEMPERATURE` vs `INITIALTEMPERATURE`);
* **how its generated copies are named**, so `.gitignore` covers them and the
  original template can never be overwritten or deleted.

Only the *conditions* differ between generated cases: the mesh
(`blockMeshDict`), the boundary-condition types, and the numerical schemes come
from the template unchanged. Adding a third case is a matter of adding an entry
here plus a `.gitignore` line — no changes to the writer, the runner, or the
tools.
"""

from __future__ import annotations

import math

from pathlib import Path
from typing import Iterable, Optional

# The case definitions themselves live in agent/cases/<name>/case.py (and on
# FGM_CASE_PATH); this module is the registry the rest of the agent queries.
from case_api import REPO_ROOT, TESTCASES_DIR, CaseTemplate, ParamSpec  # noqa: F401
from cases import discover

TEMPLATES: tuple[CaseTemplate, ...] = discover()

TEMPLATE_BY_NAME: dict[str, CaseTemplate] = {t.name: t for t in TEMPLATES}

DEFAULT_TEMPLATE: CaseTemplate = TEMPLATES[0]

# Every generated-case prefix, longest first so name resolution is unambiguous.
PREFIXES: tuple[str, ...] = tuple(
    sorted((t.prefix for t in TEMPLATES), key=len, reverse=True)
)

# The template directories themselves — never generated, never deleted.
PROTECTED_NAMES: frozenset[str] = frozenset(t.name for t in TEMPLATES)


def find_template(name: Optional[str]) -> CaseTemplate:
    """Resolve a case-type name, falling back to the default template."""
    if not name:
        return DEFAULT_TEMPLATE
    key = name.strip()
    if key in TEMPLATE_BY_NAME:
        return TEMPLATE_BY_NAME[key]
    lowered = {n.lower(): t for n, t in TEMPLATE_BY_NAME.items()}
    if key.lower() in lowered:
        return lowered[key.lower()]
    raise ValueError(
        f"unknown case type {name!r}; valid types are: "
        + ", ".join(TEMPLATE_BY_NAME)
    )


def template_for_case_name(case_name: str) -> Optional[CaseTemplate]:
    """Return the template a *generated* case name belongs to, if any.

    Used so `run_case`/`delete_cases` can act on a case by name alone, without
    being told which template it came from.
    """
    name = (case_name or "").strip()
    for template in sorted(TEMPLATES, key=lambda t: len(t.prefix), reverse=True):
        if name.startswith(template.prefix):
            return template
    return None


def domain_length(case_type: str, params: dict) -> Optional[float]:
    """Streamwise domain length [m] from the case's own hook, if it has one."""
    template = TEMPLATE_BY_NAME.get(case_type)
    if template is None or template.domain_length is None:
        return None
    return template.domain_length(params)


def discriminating_parameters(template: CaseTemplate) -> frozenset[str]:
    """Parameters that belong to ``template`` and to no other template."""
    others: set[str] = set()
    for other in TEMPLATES:
        if other.name != template.name:
            others.update(other.parameters)
    return frozenset(set(template.parameters) - others)


def infer_template(supplied: Iterable[str]) -> Optional[CaseTemplate]:
    """Infer the case type from the parameters the user actually gave.

    "set up a case, inlet 0.2 m/s, plate 300 K" can only mean the slot burner —
    a flamelet has neither — so asking which case type would be pedantic. Only
    returns a template when exactly one is compatible; ambiguous input (say,
    just a temperature and a time) returns None so the caller can ask.
    """
    names = set(supplied)
    matches = [
        template for template in TEMPLATES
        if names & discriminating_parameters(template)
    ]
    return matches[0] if len(matches) == 1 else None


def timestep_advice(
    template: CaseTemplate, params: dict[str, float]
) -> Optional[dict]:
    """Warn when the recommended time step is likely too large for the conditions.

    The solver runs at a FIXED step (`adjustTimeStep no`, and the `maxCo` entry
    in controlDict is inert), so nothing throttles it at run time: too large a
    step does not self-correct, it diverges. The recommended values were chosen
    at the reference temperature; a hotter unburnt mixture burns faster and
    thinner, so the same step resolves less of the front.

    Returns None when the step was given explicitly, or when the conditions are
    close enough to the reference for the recommendation to stand.
    """
    temperature = params.get("T_inlet")
    if temperature is None:
        return None
    reference = template.reference_temperature
    recommended = template.defaults.get("sim_timestep")
    if recommended is None or temperature <= reference * 1.15:
        return None
    # Conservative starting point, not a derived criterion: the laminar flame
    # speed of methane-air rises steeply with unburnt temperature (roughly
    # T^1.5-2), so the step is scaled down by the square of the ratio.
    suggested = recommended * (reference / temperature) ** 2
    # One significant figure — this is guidance, not a computed requirement.
    magnitude = 10 ** math.floor(math.log10(suggested))
    suggested = math.floor(suggested / magnitude) * magnitude
    return {
        "parameter": "sim_timestep",
        "recommended_at_reference": recommended,
        "reference_temperature": reference,
        "suggested": suggested,
        "message": (
            f"The default time step {recommended:g} s is the value this case "
            f"ships with at {reference:g} K. At {temperature:g} K the mixture "
            "burns faster and the flame front is thinner, so that step resolves "
            f"less of it. Consider {suggested:g} s as a starting point. The "
            "solver runs at a fixed step (adjustTimeStep no), so it will not "
            "adapt on its own — check the Courant number reported in "
            "log.FGMFoam and reduce the step further if it approaches 1."
        ),
    }


def describe_templates() -> str:
    """Render the template list for the user (the `/cases` REPL command)."""
    width = max(len(t.name) for t in TEMPLATES)
    lines = ["Case templates the agent can set up:", ""]
    for template in TEMPLATES:
        lines.append(f"  {template.name.ljust(width)}  {template.summary}")
        lines.append(f"  {' ' * width}  needs: {', '.join(template.required)}")
        recommended = ", ".join(
            f"{k}={v:g}" for k, v in template.defaults.items()
        )
        if recommended:
            lines.append(
                f"  {' ' * width}  recommended: {recommended} "
                f"(at {template.reference_temperature:g} K)"
            )
        lines.append("")
    lines.append("  Build one with: \"set up a <name> case with ...\"")
    lines.append("  Generated cases are written to testcases/<name>_<descriptor>/.")
    return "\n".join(lines)
