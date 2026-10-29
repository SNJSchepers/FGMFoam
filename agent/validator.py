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

"""Parameter schema, validation and bookkeeping for the simulation setup.

This module is the single source of truth for which parameters the agent
understands, their types, units, defaults and whether they are mandatory.
Both the tools layer (to describe the build_case parameters to the model) and
the case builder (to validate and merge values) rely on it.

Which parameters actually *apply* depends on the case template: the 1D flamelet
has no inlet and no plate, so `inlet_velocity` and `T_plate` are neither asked
for nor required there. The per-template parameter sets live in
:mod:`case_templates`; this module only knows how to validate against one.
"""

from __future__ import annotations

from typing import Any, Optional

from case_api import COMMON_PARAM_SPECS, ParamSpec
from case_templates import DEFAULT_TEMPLATE, TEMPLATES, CaseTemplate

# Friendly synonyms accepted for a parameter name. The unburnt-gas temperature
# is the inlet temperature in the slot burner but the initial temperature of a
# quiescent mixture in the flamelet, so both spellings map to one parameter.
ALIASES: dict[str, str] = {
    "T_initial": "T_inlet",
    "T_init": "T_inlet",
    "T_unburnt": "T_inlet",
    "T_0": "T_inlet",
}


def canonical_name(name: str) -> str:
    """Resolve a parameter alias to its canonical name."""
    return ALIASES.get(name, name)


def _all_specs() -> list[ParamSpec]:
    """The common parameters, then each case plugin's own, first definition wins."""
    specs: list[ParamSpec] = list(COMMON_PARAM_SPECS)
    seen = {spec.name for spec in specs}
    for template in TEMPLATES:
        for spec in template.param_specs:
            if spec.name not in seen:
                specs.append(spec)
                seen.add(spec.name)
    return specs


# Every parameter any case understands. A case's own parameter set
# (template.parameters) selects from these.
PARAM_SPECS: list[ParamSpec] = _all_specs()
# Fast lookup by name.
SPEC_BY_NAME: dict[str, ParamSpec] = {spec.name: spec for spec in PARAM_SPECS}


def required_names(template: CaseTemplate = DEFAULT_TEMPLATE) -> list[str]:
    """Return the mandatory parameter names for ``template``.

    The template decides: a parameter is only required if the template declares
    it, so the flamelet never demands an inlet velocity it has no inlet for.
    """
    return [spec.name for spec in PARAM_SPECS if spec.name in template.required]


def applicable_specs(template: CaseTemplate = DEFAULT_TEMPLATE) -> list[ParamSpec]:
    """Return the specs of every parameter that applies to ``template``."""
    return [spec for spec in PARAM_SPECS if spec.name in template.parameters]


def find_missing(
    params: dict[str, Any], template: CaseTemplate = DEFAULT_TEMPLATE
) -> list[str]:
    """Return the names of required parameters that are absent or null."""
    return [
        name
        for name in required_names(template)
        if params.get(name) is None
    ]


def missing_specs(
    params: dict[str, Any], template: CaseTemplate = DEFAULT_TEMPLATE
) -> list[ParamSpec]:
    """Return the :class:`ParamSpec` of each still-missing required parameter."""
    return [SPEC_BY_NAME[name] for name in find_missing(params, template)]


def merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    """Overlay ``updates`` onto ``base``, ignoring unknown keys and nulls.

    Values already present in ``base`` are only overwritten when the update
    provides a non-null value, so successive follow-up answers accumulate.
    """
    merged = dict(base)
    for name, value in updates.items():
        name = canonical_name(name)
        if name in SPEC_BY_NAME and value is not None:
            merged[name] = value
    return merged


def coerce_types(params: dict[str, Any]) -> dict[str, Any]:
    """Cast values to the type declared in the schema.

    Keys that cannot be coerced are dropped so that they are treated as
    "not yet provided" and get asked about again, rather than crashing later.
    """
    coerced: dict[str, Any] = {}
    for name, value in params.items():
        name = canonical_name(name)
        spec = SPEC_BY_NAME.get(name)
        if spec is None or value is None:
            continue
        try:
            if spec.py_type is int:
                coerced[name] = int(round(float(value)))
            else:
                coerced[name] = float(value)
        except (TypeError, ValueError):
            continue
    return coerced


def apply_defaults(
    params: dict[str, Any], template: CaseTemplate = DEFAULT_TEMPLATE
) -> dict[str, Any]:
    """Fill in defaults for optional parameters that apply to ``template``.

    Parameters the template does not use are left out entirely rather than
    defaulted, so a flamelet case never carries a meaningless inlet velocity.
    """
    complete = dict(params)
    for spec in applicable_specs(template):
        if complete.get(spec.name) is not None:
            continue
        # The template's own recommended value wins over the generic default:
        # the two cases need different time steps (1e-5 s vs 2e-6 s).
        fallback = template.defaults.get(spec.name, spec.default)
        if fallback is not None:
            complete[spec.name] = fallback
    return {
        name: value
        for name, value in complete.items()
        if name in template.parameters
    }


# Unit slips worth naming when a value is outside its usual range: the factor
# that brings it back in, and what the user probably meant.
_UNIT_HINTS: dict[str, tuple[tuple[float, str], ...]] = {
    "m/s": ((0.01, "cm/s"), (1 / 3.6, "km/h"), (0.001, "mm/s")),
    "m": ((0.001, "mm"), (0.01, "cm"), (1e-6, "micrometres"), (1000.0, "mm converted twice")),
}


def unusual_values(params: dict[str, Any], template: Optional[CaseTemplate],
                   specs: Optional[list[ParamSpec]] = None) -> list[dict[str, Any]]:
    """Values outside the usual range of their parameter in this case.

    Not errors: the user may mean them. Each entry carries the usual range and,
    when a common unit slip would explain the value, what was probably meant.
    ``specs`` replaces the template's parameters (before a case type is known).
    """
    out = []
    for spec in specs if specs is not None else applicable_specs(template):
        value = params.get(spec.name)
        if spec.typical is None or value is None:
            continue
        values = value if isinstance(value, list) else [value]
        low, high = spec.typical
        for v in values:
            try:
                x = float(v)
            except (TypeError, ValueError):
                continue
            if low <= x <= high:
                continue
            entry: dict[str, Any] = {"parameter": spec.name, "value": x,
                                     "unit": spec.unit, "usual_range": [low, high]}
            for factor, meaning in _UNIT_HINTS.get(spec.unit, ()):
                if low <= x * factor <= high:
                    entry["maybe_meant"] = f"{x * factor:g} {spec.unit} (if the value was in {meaning})"
                    break
            out.append(entry)
    return out

