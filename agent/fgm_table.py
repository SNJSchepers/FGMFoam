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

"""Read an FGM table's header and control-variable ranges from Python.

The solver parses `.fgm` files in C (`solver/src/lookUp/FGM/FGMlib.c`); this is a
small read-only counterpart for the agent, used to record a table's identity in a
case's provenance and to check *before* a run that the case's conditions fall
inside what the table actually covers.

The format is keyword-delimited ASCII:

    [FGM] / NDCL            header marker
    [PRESSURE]   101325     operating pressure the table was generated at [Pa]
    [DIMENSION]  2          number of control variables
    [DATASIZE]   250 100    grid points per control variable
    [VARIABLES]  16 + names the columns of each data row, CV1 and CV2 first
    [END]
    [DATA]                  one row per grid point, row-major with CV1 fastest

Only the header and the first two columns matter here, but the data block is
parsed in full (25k x 16 floats is trivial) and cached per file identity, since
both callers want the CV ranges.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class TableInfo:
    """Everything the agent needs to know about an FGM table."""

    path: Path
    sha256: str
    pressure: float
    dimension: int
    datasize: tuple[int, ...]
    variables: tuple[str, ...]
    # Range of each control variable, in the table's own order (CV1, CV2, ...).
    cv_ranges: tuple[tuple[float, float], ...]
    n_rows: int

    @property
    def cv_names(self) -> tuple[str, ...]:
        return self.variables[: self.dimension]

    def range_of(self, index: int) -> tuple[float, float]:
        return self.cv_ranges[index]


# Parsed tables, keyed by (path, mtime_ns, size) so an edited table is re-read.
_CACHE: dict[tuple[str, int, int], TableInfo] = {}


def sha256_of(path: Path) -> str:
    """Content hash, so a provenance record pins the exact table used."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _section(lines: list[str], keyword: str) -> int:
    """Index of the line after ``[KEYWORD]``, or -1 when absent."""
    for i, line in enumerate(lines):
        if line.strip() == keyword:
            return i + 1
    return -1


def read_table(path: Path) -> TableInfo:
    """Parse an FGM table's header and control-variable ranges (cached)."""
    path = Path(path)
    stat = path.stat()
    key = (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
    cached = _CACHE.get(key)
    if cached is not None:
        return cached

    text = path.read_text()
    lines = text.splitlines()

    def scalar(keyword: str) -> str:
        index = _section(lines, keyword)
        if index < 0:
            raise ValueError(f"{path.name}: missing {keyword} section")
        return lines[index].strip()

    pressure = float(scalar("[PRESSURE]"))
    dimension = int(scalar("[DIMENSION]"))
    datasize = tuple(int(n) for n in scalar("[DATASIZE]").split())

    start = _section(lines, "[VARIABLES]")
    if start < 0:
        raise ValueError(f"{path.name}: missing [VARIABLES] section")
    n_vars = int(lines[start].strip())
    variables = tuple(line.strip() for line in lines[start + 1 : start + 1 + n_vars])

    data_start = _section(lines, "[DATA]")
    if data_start < 0:
        raise ValueError(f"{path.name}: missing [DATA] section")
    # The data block runs until the next bracketed keyword — the shipped tables
    # close with a trailing [END] after the rows, which numpy would choke on.
    rows: list[str] = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("["):
            break
        rows.append(stripped)
    if not rows:
        raise ValueError(f"{path.name}: [DATA] section is empty")
    # Only the control-variable columns are needed, but parsing them for every
    # row is cheap and catches a truncated table.
    data = np.loadtxt(rows, usecols=range(dimension), ndmin=2)
    expected = 1
    for n in datasize:
        expected *= n
    if data.shape[0] != expected:
        raise ValueError(
            f"{path.name}: [DATASIZE] says {expected} points but [DATA] has "
            f"{data.shape[0]} rows"
        )
    cv_ranges = tuple(
        (float(data[:, i].min()), float(data[:, i].max())) for i in range(dimension)
    )

    info = TableInfo(
        path=path,
        sha256=sha256_of(path),
        pressure=pressure,
        dimension=dimension,
        datasize=datasize,
        variables=variables,
        cv_ranges=cv_ranges,
        n_rows=int(data.shape[0]),
    )
    _CACHE[key] = info
    return info


def read_table_safe(path: Path) -> Optional[TableInfo]:
    """:func:`read_table`, returning None instead of raising on any problem."""
    try:
        return read_table(path)
    except Exception:  # noqa: BLE001 - provenance/checks must never break a build
        return None
