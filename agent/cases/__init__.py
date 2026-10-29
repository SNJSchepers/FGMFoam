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

"""Find and validate the case plugins.

Every ``<dir>/case.py`` defining ``CASE = CaseTemplate(...)`` is a case, where
``<dir>`` is a folder under ``agent/cases/`` or under any directory listed in
the ``FGM_CASE_PATH`` environment variable (``:``-separated), so your own cases
can live outside this repository.

A declaration that is internally inconsistent — a placeholder mapped to a
parameter nobody declared, a sweepable variable the case does not have — is
rejected here with the plugin's path, rather than surfacing later as a
half-built case.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

from case_api import COMMON_PARAMETERS, DERIVED_VALUES, CaseTemplate

BUILTIN_DIR: Path = Path(__file__).resolve().parent


class CasePluginError(ValueError):
    """A case plugin is missing, malformed, or inconsistent."""


def _search_dirs() -> list[Path]:
    dirs = [BUILTIN_DIR]
    for entry in os.environ.get("FGM_CASE_PATH", "").split(os.pathsep):
        if entry.strip():
            dirs.append(Path(entry).expanduser().resolve())
    return dirs


def _load(path: Path) -> CaseTemplate:
    spec = importlib.util.spec_from_file_location(f"fgm_case_{path.parent.name}", path)
    if spec is None or spec.loader is None:
        raise CasePluginError(f"{path}: cannot be imported")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        raise CasePluginError(f"{path}: failed to import: {exc}") from exc
    case = getattr(module, "CASE", None)
    if not isinstance(case, CaseTemplate):
        raise CasePluginError(f"{path}: must define CASE = CaseTemplate(...)")
    return case


def validate(case: CaseTemplate, where: str) -> None:
    """Raise CasePluginError describing every inconsistency in ``case``."""
    problems: list[str] = []
    if not case.name or "_" in case.name or "/" in case.name:
        problems.append(
            f"name {case.name!r} must be non-empty with no '_' or '/': generated "
            "cases are named '<name>_<descriptor>'"
        )
    if not case.template_dir.is_dir():
        problems.append(f"template directory not found: {case.template_dir}")
    declared = COMMON_PARAMETERS | {s.name for s in case.param_specs}
    for name in case.parameters:
        if name not in declared:
            problems.append(f"parameter {name!r} has no ParamSpec (add it to param_specs)")
    for token, source in case.placeholders.items():
        if source not in case.parameters and source not in DERIVED_VALUES:
            problems.append(f"placeholder {token} maps to {source!r}, which is not a "
                            "parameter of this case")
    for name in (*case.sweepable, *case.geometric):
        if name not in case.parameters:
            problems.append(f"{name!r} is sweepable/geometric but not a parameter")
    if case.analyse is not None and not (case.metric and case.metric.get("key")):
        problems.append("analyse is set but metric['key'] is not")
    if problems:
        raise CasePluginError(f"{where}:\n  - " + "\n  - ".join(problems))


def discover() -> tuple[CaseTemplate, ...]:
    """All valid cases, ordered by ``order`` then name."""
    found: dict[str, CaseTemplate] = {}
    for directory in _search_dirs():
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*/case.py")):
            case = _load(path)
            validate(case, str(path))
            if case.name in found:
                raise CasePluginError(f"{path}: case name {case.name!r} is already defined")
            found[case.name] = case
    if not found:
        raise CasePluginError(f"no cases found under {BUILTIN_DIR}")
    return tuple(sorted(found.values(), key=lambda c: (c.order, c.name)))
