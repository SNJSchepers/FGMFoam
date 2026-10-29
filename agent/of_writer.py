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

"""Create and patch the OpenFOAM case copy for the laminar slot burner.

This mirrors what the manual bash ``run`` script did, but in Python: it builds a
fresh working copy of the template case, drops in the FGM lookup table, and
substitutes every placeholder token with the resolved parameter values using
``str.replace``. The original template under ``testcases/laminarSlotBurner`` is
never modified.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any, Optional

import case_templates
from cantera_tool import Thermochemistry
from case_templates import CaseTemplate

# Repository root, derived from this file's location (agent/ lives at the root).
REPO_ROOT: Path = Path(__file__).resolve().parent.parent

TESTCASES_DIR: Path = REPO_ROOT / "testcases"
TABLE_SRC: Path = REPO_ROOT / "tables" / "database.fgm"

# Which template a case comes from decides its paths, placeholders and naming —
# see case_templates. These module-level constants describe the default template
# and exist so callers that never mention a case type keep working unchanged.
DEFAULT_TEMPLATE: CaseTemplate = case_templates.DEFAULT_TEMPLATE
TEMPLATE_DIR: Path = DEFAULT_TEMPLATE.template_dir
# Generated cases are named "<template>_<descriptor>" so one .gitignore rule per
# template (testcases/laminarSlotBurner_*, testcases/laminarFlamelet_*) covers
# all of them.
CASE_PREFIX: str = DEFAULT_TEMPLATE.prefix
DEFAULT_CASE_NAME: str = DEFAULT_TEMPLATE.default_case_name

# Cruft and run-time artifacts not carried into the copy. The bash ``run``/
# ``clean`` scripts are excluded on purpose (the agent now does the parameter
# substitution they did). ``foam.foam`` IS copied so the case opens in ParaView.
_IGNORE = shutil.ignore_patterns(
    "run", "clean", "log.*", "*.swp", ".DS_Store",
    "processor*", "0", "system",
)

# Case-relative files that contain placeholders to be patched. Per template, so
# a case without an inlet is not expected to carry an inlet-velocity file.
TARGET_FILES: tuple[str, ...] = DEFAULT_TEMPLATE.target_files


def _normalise_case_name(
    case_name: str, template: Optional[CaseTemplate] = None
) -> str:
    """Validate a case name and ensure it carries a generated-case prefix.

    The name must be a single safe path segment. A name that already carries any
    template's prefix is kept as-is — that is how ``run_case``/``delete_cases``
    address a case by name alone. A bare descriptor like ``lowPhi`` is prefixed
    with ``template``'s name, so the .gitignore rule covers it and the original
    template case can never be the target.
    """
    default_name = (template or DEFAULT_TEMPLATE).default_case_name
    name = (case_name or "").strip() or default_name
    if "/" in name or "\\" in name or name in (".", ".."):
        raise ValueError(f"invalid case name (no path separators): {case_name!r}")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        raise ValueError(
            f"case name may contain only letters, digits, '.', '_', '-': {case_name!r}"
        )
    if not any(name.startswith(prefix) for prefix in case_templates.PREFIXES):
        name = (template or DEFAULT_TEMPLATE).prefix + name
    return name


def case_path(
    case_name: str = "", template: Optional[CaseTemplate] = None
) -> Path:
    """Return the absolute path of a (normalised) case directory."""
    return TESTCASES_DIR / _normalise_case_name(case_name, template)


def generated_cases() -> list[str]:
    """Names of all agent-generated case directories (prefixed; never templates)."""
    if not TESTCASES_DIR.is_dir():
        return []
    return sorted(
        path.name
        for path in TESTCASES_DIR.iterdir()
        if path.is_dir()
        and any(path.name.startswith(prefix) for prefix in case_templates.PREFIXES)
    )


def delete_case(case_name: str) -> Path:
    """Delete a generated case directory and return its path.

    Refuses anything that is not a prefixed, agent-generated case directly under
    ``testcases/``. The original template cases (``laminarSlotBurner``,
    ``laminarFlamelet``) can never be a target: ``case_path`` always adds a
    template prefix, so a request to delete a template resolves to a
    (non-existent) prefixed name instead. The explicit check below makes that
    guarantee independent of how names are normalised.
    """
    target = case_path(case_name)
    is_generated = any(
        target.name.startswith(prefix) for prefix in case_templates.PREFIXES
    )
    if (
        target.parent != TESTCASES_DIR
        or not is_generated
        or target.name in case_templates.PROTECTED_NAMES
    ):
        raise ValueError(f"refusing to delete a non-generated path: {target}")
    if target.exists():
        shutil.rmtree(target)
    return target


def prepare_case(
    case_name: str = "", template: Optional[CaseTemplate] = None
) -> Path:
    """Create a clean working copy of ``template`` (default: the slot burner).

    Deletes any existing case of the same name, copies the template, materialises
    the working ``0/`` and ``system/`` folders from their ``.orig`` sources,
    installs the FGM lookup table, and ensures an empty ``foam.foam`` marker.

    Returns:
        The path to the prepared case directory.
    """
    template = template or DEFAULT_TEMPLATE
    if not template.template_dir.is_dir():
        raise FileNotFoundError(f"Template case not found: {template.template_dir}")
    if not TABLE_SRC.is_file():
        raise FileNotFoundError(f"FGM table not found: {TABLE_SRC}")

    target = case_path(case_name, template)

    # Start from scratch every time.
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(template.template_dir, target, ignore=_IGNORE)

    # Build the working field/system folders from the .orig templates.
    _materialise(target, "0")
    _materialise(target, "system")

    # Install the lookup table into constant/lookUp/.
    lookup_dir = target / "constant" / "lookUp"
    lookup_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(TABLE_SRC, lookup_dir / TABLE_SRC.name)

    # Empty ParaView marker (copied from the template; guaranteed regardless).
    (target / "foam.foam").touch()

    return target


def build_replacements(
    params: dict[str, Any],
    thermo: Thermochemistry,
    template: Optional[CaseTemplate] = None,
) -> dict[str, str]:
    """Map each of ``template``'s placeholder tokens to its formatted value.

    The token set is per template — the slot burner spells the unburnt
    temperature ``INLETTEMPERATURE``, the flamelet ``INITIALTEMPERATURE`` — but
    both draw from the same parameters and the same Cantera result.
    """
    template = template or DEFAULT_TEMPLATE
    sources: dict[str, Any] = {
        **params,
        "Yc_u": thermo.Yc_u,
        "Yc_b": thermo.Yc_b,
        "ht": thermo.ht,
        # A quotient: 12 significant figures drop the float noise (6e-06, not
        # 5.999999999999999e-06) and keep all the precision that matters.
        "write_time": float(f"{params['sim_time'] / params['sim_write_steps']:.12g}"),
    }
    missing = [
        key for key in template.placeholders.values() if key not in sources
    ]
    if missing:
        raise KeyError(
            f"{template.name}: no value for {', '.join(sorted(missing))}"
        )
    return {
        token: _format(sources[key])
        for token, key in template.placeholders.items()
    }


def patch_placeholders(
    case_dir: Path,
    replacements: dict[str, str],
    template: Optional[CaseTemplate] = None,
) -> None:
    """Replace every placeholder token in the case's templated files."""
    template = template or DEFAULT_TEMPLATE
    # Apply longer tokens first so no token can be a prefix of another.
    ordered = sorted(replacements.items(), key=lambda kv: len(kv[0]), reverse=True)
    for rel_path in template.target_files:
        file_path = case_dir / rel_path
        if not file_path.is_file():
            continue  # a template need not use every templated file
        text = file_path.read_text()
        for placeholder, value in ordered:
            text = text.replace(placeholder, value)
        file_path.write_text(text)


def read_control_entries(
    case_dir: Path, keys: list[str]
) -> dict[str, Optional[str]]:
    """Return the current value text of the given system/controlDict keys.

    A key not present in the file maps to ``None``. The value is the text between
    the key and its terminating ``;`` (stripped), e.g. ``endTime`` -> ``"0.1"``.
    """
    text = (case_dir / "system" / "controlDict").read_text()
    out: dict[str, Optional[str]] = {}
    for key in keys:
        match = re.search(rf"^\s*{re.escape(key)}\s+([^;\n]+);", text, re.MULTILINE)
        out[key] = match.group(1).strip() if match else None
    return out


def set_control_entries(
    case_dir: Path, entries: dict[str, Any]
) -> dict[str, dict[str, Optional[str]]]:
    """Set keys in an existing case's system/controlDict; return per-key old/new.

    A float/int value is rendered with :func:`_format`; a string (e.g. the word
    ``latestTime``) is written verbatim. Each key must already exist in the file
    (these are run-control edits to a built case, not new entries) — a missing key
    raises ``KeyError``. The mesh, fields and thermochemistry are untouched.
    """
    path = case_dir / "system" / "controlDict"
    text = path.read_text()
    changes: dict[str, dict[str, Optional[str]]] = {}
    for key, value in entries.items():
        new = value if isinstance(value, str) else _format(value)
        match = re.search(
            rf"^(\s*{re.escape(key)}\s+)([^;\n]+)(;)", text, re.MULTILINE
        )
        if match is None:
            raise KeyError(f"controlDict has no '{key}' entry to set")
        changes[key] = {"old": match.group(2).strip(), "new": new}
        text = text[: match.start()] + match.group(1) + new + match.group(3) + text[match.end() :]
    path.write_text(text)
    return changes


def _materialise(case_dir: Path, name: str) -> None:
    """Create ``case_dir/<name>`` from ``case_dir/<name>.orig``."""
    source = case_dir / f"{name}.orig"
    work = case_dir / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(source, work)


def _format(value: Any) -> str:
    """Render a numeric value as an OpenFOAM-friendly token."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    # repr gives the shortest round-trippable form, preserving full precision.
    return repr(float(value))
