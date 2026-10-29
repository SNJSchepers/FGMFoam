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

"""Write a provenance record into every generated case.

A built case is a directory of numbers with no memory of where they came from.
Months later — or to a reviewer, or a co-author — "an assistant set this up" is
not an answer. `agent-provenance.json` makes a case self-describing: the exact
parameters and **which of them the user actually chose**, the thermochemistry and
the mechanism it came from, the identity (content hash) of the FGM table, the
agent's own version and git commit, and which model routed the request.

Written as JSON rather than an OpenFOAM dictionary or YAML on purpose: no extra
dependency, no quoting or escaping to get subtly wrong, diffable between cases,
and directly loadable by whatever post-processing comes later. OpenFOAM ignores
it — it is a record, not an input.

Nothing here may break a build: every lookup that can fail is best-effort and
degrades to a null field.
"""

from __future__ import annotations

import json
import platform
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

import fgm_table
from case_templates import CaseTemplate

# Bump when the *shape* of the record changes, so later readers can adapt.
SCHEMA: str = "fgmfoam-agent-provenance/1"

FILENAME: str = "agent-provenance.json"


def _git(repo_root: Path, *args: str) -> Optional[str]:
    """Run a read-only git command, returning None if git or the repo is absent."""
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root), *args],
            capture_output=True, text=True, timeout=5, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return out.stdout.strip() or None


def _agent_identity(repo_root: Path, version: str) -> dict[str, Any]:
    commit = _git(repo_root, "rev-parse", "--short", "HEAD")
    status = _git(repo_root, "status", "--porcelain")
    return {
        "version": version,
        "git_commit": commit,
        # A dirty tree means the commit alone does not identify the code that ran.
        "git_dirty": None if status is None else bool(status),
        "python": platform.python_version(),
    }


def _thermochemistry(thermo: Any) -> dict[str, Any]:
    try:
        import cantera

        cantera_version = cantera.__version__
    except Exception:  # noqa: BLE001
        cantera_version = None
    return {
        "Yc_u": thermo.Yc_u,
        "Yc_b": thermo.Yc_b,
        "ht": thermo.ht,
        "units": {"Yc": "kmol/kg", "ht": "J/kg"},
        "definition": "Yc = Y_CO2 / M_CO2",
        "mechanism": "gri30.yaml",
        "cantera_version": cantera_version,
    }


def _relative(path: Path, repo_root: Path) -> str:
    """Repo-relative path where possible — a record may be shared or published,
    and an absolute path leaks the author's home directory."""
    try:
        return str(Path(path).resolve().relative_to(Path(repo_root).resolve()))
    except ValueError:
        return str(path)


def _table(table_path: Path, repo_root: Path) -> dict[str, Any]:
    info = fgm_table.read_table_safe(table_path)
    rel = _relative(table_path, repo_root)
    if info is None:
        return {"path": rel, "sha256": None, "note": "table unreadable"}
    return {
        "path": rel,
        "sha256": info.sha256,
        "pressure": info.pressure,
        "dimension": info.dimension,
        "datasize": list(info.datasize),
        "control_variables": list(info.cv_names),
        "control_variable_ranges": {
            name: {"min": lo, "max": hi}
            for name, (lo, hi) in zip(info.cv_names, info.cv_ranges)
        },
        "variables": list(info.variables),
    }


def _parameter_sources(
    params: dict[str, Any],
    supplied: set[str],
    template: CaseTemplate,
    estimated: Optional[set[str]] = None,
    inherited: Optional[set[str]] = None,
    shipped: Optional[set[str]] = None,
) -> dict[str, str]:
    """Label each parameter "user", "inherited", "shipped", "estimated",
    "recommended" or "default".

    The distinction is the point of the record: a reviewer needs to see which
    numbers were a deliberate choice and which the agent filled in.
    """
    estimated = estimated or set()
    inherited = inherited or set()
    shipped = shipped or set()
    sources: dict[str, str] = {}
    for name in params:
        if name in supplied:
            sources[name] = "user"
        elif name in inherited:
            # Chosen by the user for the previous case and not changed since.
            sources[name] = "inherited"
        elif name in shipped:
            # The template's own run-script value, used because the user asked
            # for the default conditions.
            sources[name] = "shipped"
        elif name in estimated:
            # Computed from the mesh and the flow for this case, not a default.
            sources[name] = "estimated"
        elif name in template.defaults:
            sources[name] = "recommended"
        else:
            sources[name] = "default"
    return sources


def build_record(
    *,
    case_dir: Path,
    template: CaseTemplate,
    params: dict[str, Any],
    supplied: set[str],
    thermo: Any,
    estimated: Optional[set[str]] = None,
    inherited: Optional[set[str]] = None,
    table_path: Path,
    repo_root: Path,
    version: str,
    model: Optional[dict[str, Any]] = None,
    warnings: Optional[list[dict[str, Any]]] = None,
    shipped: Optional[set[str]] = None,
) -> dict[str, Any]:
    """Assemble the provenance record for a freshly built case."""
    return {
        "schema": SCHEMA,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "agent": _agent_identity(repo_root, version),
        # Which model routed the request. The deterministic pipeline below does
        # not depend on it — recorded so an unexpected case can be traced anyway.
        "model": model or {"provider": None, "chat_model": None},
        "case": {
            "name": case_dir.name,
            "type": template.name,
            "template": _relative(template.template_dir, repo_root),
        },
        "parameters": dict(params),
        "parameter_sources": _parameter_sources(
            params, supplied, template, estimated, inherited, shipped
        ),
        "derived": {
            "write_interval": params["sim_time"] / params["sim_write_steps"],
            "n_timesteps": round(params["sim_time"] / params["sim_timestep"]),
        },
        "thermochemistry": _thermochemistry(thermo),
        "fgm_table": _table(table_path, repo_root),
        "warnings": warnings or [],
    }


def write(case_dir: Path, record: dict[str, Any]) -> Path:
    """Write ``record`` into the case directory; returns the file path."""
    path = case_dir / FILENAME
    path.write_text(json.dumps(record, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return path


def read(case_dir: Path) -> Optional[dict[str, Any]]:
    """Read a case's provenance record, or None when it has none."""
    path = Path(case_dir) / FILENAME
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
