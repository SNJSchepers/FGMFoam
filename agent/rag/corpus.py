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

"""Collect the source files that make up the agent's knowledge corpus.

The corpus is the FGMFoam solver source, the project/solver docs, the example
case directories (text files only), and the text header of the FGM lookup table
(its bulk numeric data is excluded). Build artifacts, the large polyMesh
geometry, logs, and editor cruft are excluded.

Everything stays text-based and small (a few hundred KB), so the whole corpus is
read into memory at index-build time.
"""

from __future__ import annotations

import case_templates
import hashlib
from dataclasses import dataclass
from pathlib import Path

from config import REPO_ROOT

# Glob patterns (relative to the repository root) that may contribute to the
# corpus. Exclusions below are applied on top of these.
SOURCE_GLOBS: list[str] = [
    "solver/**/*.C",
    "solver/**/*.H",
    "solver/**/*.c",  # lowercase: e.g. the FGM lookup library FGMlib.c
    "solver/**/*.h",
    "solver/**/README",
    "README.md",
    "agent/README.md",
    "agent/knowledge/**/*.md",
    # Notes shipped with each case plugin, and the declarations themselves.
    "agent/cases/*/*.md",
    "agent/cases/*/case.py",
    "testcases/**/*",
]

# Directory names that disqualify a path anywhere in its parts.
EXCLUDE_DIR_PARTS: frozenset[str] = frozenset({
    ".git", ".venv", "__pycache__", "Make", "lnInclude", ".vscode",
})

# Generated cases (build_case output) are never indexed: they are rebuilt from
# scratch on every build, and indexing them would fill the corpus with near
# duplicates of the templates and mark the index stale after each one. The
# prefixes come from the template registry so a new case type is covered
# automatically — hardcoding one prefix here previously let generated
# laminarFlamelet cases slip into the index.
GENERATED_CASE_PREFIXES: tuple[str, ...] = case_templates.PREFIXES

# polyMesh geometry files are huge numeric data; keep only the small "boundary".
EXCLUDE_NAMES: frozenset[str] = frozenset({
    "points", "owner", "neighbour", "faces",
    ".DS_Store", "foam.foam",
})

EXCLUDE_SUFFIXES: frozenset[str] = frozenset({
    ".fgm", ".o", ".dep", ".swp", ".so", ".png", ".jpg", ".pdf",
})

# FGM lookup tables: only their text header is indexed (see _read_fgm_header).
FGM_GLOBS: list[str] = ["tables/*.fgm", "testcases/**/*.fgm"]
FGM_HEADER_MAX_LINES: int = 100


@dataclass
class SourceFile:
    """A single text file in the corpus."""

    rel_path: str   # path relative to the repository root
    text: str
    kind: str       # "cpp" | "markdown" | "openfoam" | "text" | "fgm"
    sha256: str


def _read_fgm_header(path: Path, max_lines: int = FGM_HEADER_MAX_LINES) -> str | None:
    """Return an .fgm table's text header, up to and including the [DATA] marker.

    The numeric data after [DATA] (hundreds of thousands of values) is skipped;
    only the metadata — pressure, dimensions, and tabulated variables — is kept.
    """
    lines: list[str] = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for i, line in enumerate(handle):
                stripped = line.rstrip("\n")
                lines.append(stripped)
                if stripped.strip() == "[DATA]" or i + 1 >= max_lines:
                    break
    except OSError:
        return None
    header = "\n".join(lines).strip()
    return header or None


def _kind_of(rel_path: str) -> str:
    """Classify a file so the chunker can pick an appropriate strategy."""
    lower = rel_path.lower()
    if lower.endswith((".c", ".h")):  # OpenFOAM C++ uses .C / .H
        return "cpp"
    if lower.endswith(".md"):
        return "markdown"
    if rel_path.startswith("testcases/"):
        return "openfoam"
    return "text"


def _in_excluded_dir(path: Path) -> bool:
    """True if any path component is an excluded directory or a generated case."""
    return any(
        part in EXCLUDE_DIR_PARTS
        or any(part.startswith(prefix) for prefix in GENERATED_CASE_PREFIXES)
        for part in path.parts
    )


def _is_excluded(path: Path) -> bool:
    if _in_excluded_dir(path):
        return True
    if path.name in EXCLUDE_NAMES:
        return True
    if path.suffix in EXCLUDE_SUFFIXES:
        return True
    if path.name.startswith("log."):
        return True
    return False


def collect_corpus(repo_root: Path = REPO_ROOT) -> list[SourceFile]:
    """Gather all corpus files, sorted by path for reproducible ordering."""
    seen: set[Path] = set()
    files: list[SourceFile] = []

    candidates: list[Path] = []
    for pattern in SOURCE_GLOBS:
        candidates.extend(repo_root.glob(pattern))

    for path in sorted(set(candidates)):
        if not path.is_file() or path in seen or _is_excluded(path):
            continue
        seen.add(path)
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            # Skip anything that is not valid UTF-8 text (stray binaries).
            continue
        if not text.strip():
            continue
        rel_path = path.relative_to(repo_root).as_posix()
        files.append(
            SourceFile(
                rel_path=rel_path,
                text=text,
                kind=_kind_of(rel_path),
                sha256=hashlib.sha256(text.encode("utf-8")).hexdigest(),
            )
        )

    # FGM lookup tables: index only the header. The master copy and the
    # identical per-case copies are de-duplicated by header content.
    seen_headers: set[str] = set()
    fgm_candidates: list[Path] = []
    for pattern in FGM_GLOBS:
        fgm_candidates.extend(repo_root.glob(pattern))
    for path in sorted(set(fgm_candidates)):
        if not path.is_file() or _in_excluded_dir(path):
            continue
        header = _read_fgm_header(path)
        if not header or header in seen_headers:
            continue
        seen_headers.add(header)
        rel_path = path.relative_to(repo_root).as_posix()
        files.append(
            SourceFile(
                rel_path=rel_path,
                text=header,
                kind="fgm",
                sha256=hashlib.sha256(header.encode("utf-8")).hexdigest(),
            )
        )

    return files
