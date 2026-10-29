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

"""Is the FGMFoam build on the path the one this repository's cases need?

The solver and its libraries are found through ``PATH`` and ``LD_LIBRARY_PATH``,
so a different FGMFoam build earlier on those paths (another checkout sharing
``FOAM_USER_LIBBIN``) is what a run loads. Verified: such a build's table reader
stopped every case of a campaign at ``decomposePar`` with "Couldn't find
keyword DATASIZE", a second after launch.

The check is what actually failed, asked directly: does the ``liblookUp.so`` a
run would load provide the C functions this repository's solver calls into it?
Those are read from this repository's own ``FGMlib.h`` and the sources that use
it, so the check follows the code. No solver is run.
"""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path
from typing import Optional

from config import REPO_ROOT

SOLVER_DIR = REPO_ROOT / "solver"
FGM_HEADER = SOLVER_DIR / "src" / "lookUp" / "FGM" / "FGMlib.h"
LIBRARIES = ("liblookUp.so", "libcombustionModelFGM.so", "liblookUpBoundaryCondition.so")

_DECLARATION = re.compile(r"^\s*[A-Za-z_][\w\s\*]*?\b([A-Za-z_]\w*)\s*\(", re.MULTILINE)
_CACHE: dict[tuple, Optional[str]] = {}


def required_functions() -> list[str]:
    """The FGMlib.h functions this repository's solver calls (readFGM, ...)."""
    try:
        declared = set(_DECLARATION.findall(FGM_HEADER.read_text(errors="replace")))
    except OSError:
        return []
    used: set[str] = set()
    for path in (SOLVER_DIR / "src").rglob("*"):
        if path.suffix in (".C", ".H") and "lookUp/FGM" not in str(path):
            text = path.read_text(errors="replace")
            used |= {name for name in declared if re.search(rf"\b{name}\s*\(", text)}
    return sorted(used)


def _find_library(name: str) -> Optional[Path]:
    """The file the dynamic linker would load for ``name`` (LD_LIBRARY_PATH order)."""
    for directory in os.environ.get("LD_LIBRARY_PATH", "").split(":"):
        candidate = Path(directory) / name
        if directory and candidate.is_file():
            return candidate
    return None


def _exports(library: Path, names: list[str]) -> set[str]:
    """Which of ``names`` the library defines (its dynamic symbol table)."""
    nm = shutil.which("nm")
    if nm:
        import subprocess
        try:
            out = subprocess.run([nm, "-D", "--defined-only", str(library)],
                                 capture_output=True, text=True, timeout=20).stdout
            defined = {line.split()[-1] for line in out.splitlines() if line.strip()}
            return {n for n in names if n in defined}
        except (OSError, subprocess.SubprocessError):
            pass
    # Without nm: a defined C symbol's name is in the dynamic string table.
    data = library.read_bytes()
    return {n for n in names if b"\0" + n.encode() + b"\0" in data}


def problem() -> Optional[str]:
    """Why a run would not use a solver build matching this repository, or None."""
    key = (os.environ.get("PATH", ""), os.environ.get("LD_LIBRARY_PATH", ""))
    if key in _CACHE:
        return _CACHE[key]
    _CACHE[key] = result = _problem()
    return result


def _problem() -> Optional[str]:
    fix = ("Put this repository's build first on the paths before starting the agent, "
           "e.g. export FOAM_USER_LIBBIN=<its lib dir> FOAM_USER_APPBIN=<its bin dir> "
           "LD_LIBRARY_PATH=<its lib dir>:$LD_LIBRARY_PATH PATH=<its bin dir>:$PATH, "
           "or build it: cd solver && source bashrc && ./Allwmake.")
    solver = shutil.which("FGMFoam")
    if solver is None:
        return ("FGMFoam is not on the PATH, so no case can run. Source OpenFOAM and "
                "build the solver (README: Build).")
    libraries = {name: _find_library(name) for name in LIBRARIES}
    missing = [name for name, path in libraries.items() if path is None]
    if missing:
        return (f"{', '.join(missing)} not found on LD_LIBRARY_PATH, so the solver "
                f"cannot load its FGM model. {fix}")
    needed = required_functions()
    lookup = libraries["liblookUp.so"]
    absent = sorted(set(needed) - _exports(lookup, needed))
    if absent:
        return (f"The FGM table library a run would load ({lookup}) is from a different "
                f"FGMFoam build: it lacks {', '.join(absent)}, which this repository's "
                f"solver calls, so the table cannot be read (the run stops at "
                f"decomposePar). FGMFoam on the PATH: {solver}. {fix}")
    return None
