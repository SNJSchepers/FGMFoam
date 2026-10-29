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

"""Read what the agent needs from a built `constant/polyMesh`.

Only point coordinates and the cell count are needed — for the pre-run Courant
estimate and to place the sampling line used by the post-run analysis. Cell
connectivity is deliberately not reconstructed: OpenFOAM's own `postProcess`
does the interpolation, and a partial reimplementation would only be a source of
subtly wrong numbers.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

_SCALAR = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
_POINT = re.compile(rf"\(\s*({_SCALAR})\s+({_SCALAR})\s+({_SCALAR})\s*\)")

# Parsed point lists, keyed by (path, mtime_ns, size).
_CACHE: dict[tuple[str, int, int], list[tuple[float, float, float]]] = {}


def points_file(case_dir: Path) -> Path:
    return Path(case_dir) / "constant" / "polyMesh" / "points"


def read_points(case_dir: Path) -> Optional[list[tuple[float, float, float]]]:
    """Every mesh vertex, or None when the mesh has not been built."""
    path = points_file(case_dir)
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
    cached = _CACHE.get(key)
    if cached is not None:
        return cached
    try:
        text = path.read_text()
    except OSError:
        return None
    points = [
        (float(a), float(b), float(c)) for a, b, c in _POINT.findall(text)
    ]
    if not points:
        return None
    _CACHE[key] = points
    return points


def bounds(case_dir: Path) -> Optional[tuple[tuple[float, float], ...]]:
    """Mesh bounding box as ``((xmin, xmax), (ymin, ymax), (zmin, zmax))``."""
    points = read_points(case_dir)
    if points is None:
        return None
    return tuple(
        (min(p[axis] for p in points), max(p[axis] for p in points))
        for axis in range(3)
    )


def min_spacing(case_dir: Path) -> Optional[float]:
    """Smallest gap between distinct coordinates, over all three axes.

    Exact for the axis-aligned block meshes these cases use, and it reflects the
    real grading — which matters, because both shipped meshes are strongly
    graded towards the flame and an average would flatter any Courant estimate.
    """
    points = read_points(case_dir)
    if points is None:
        return None
    smallest: Optional[float] = None
    for axis in range(3):
        values = sorted({round(p[axis], 12) for p in points})
        for a, b in zip(values, values[1:]):
            gap = b - a
            if gap > 0 and (smallest is None or gap < smallest):
                smallest = gap
    return smallest


def cell_count(case_dir: Path) -> Optional[int]:
    """Cell count from the polyMesh owner header's ``note`` line."""
    try:
        with (Path(case_dir) / "constant" / "polyMesh" / "owner").open() as handle:
            head = handle.read(2048)
    except OSError:
        return None
    match = re.search(r"nCells:\s*(\d+)", head)
    return int(match.group(1)) if match else None
