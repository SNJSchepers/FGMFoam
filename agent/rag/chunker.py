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

"""Split source files into retrievable chunks with citable metadata.

Code, OpenFOAM dictionaries, and plain text are chunked with overlapping line
windows. Markdown is split on headings first, then any oversized section is
windowed the same way. Every chunk records its file and line range so the agent
can cite sources as ``path:start-end``.
"""

from __future__ import annotations

from dataclasses import dataclass

from rag.corpus import SourceFile

# Default windowing parameters (in lines).
WINDOW: int = 40
OVERLAP: int = 10


@dataclass
class Chunk:
    """A retrievable piece of one source file."""

    rel_path: str
    start_line: int  # 1-based, inclusive
    end_line: int    # 1-based, inclusive
    kind: str
    text: str

    @property
    def citation(self) -> str:
        """Human-readable source reference, e.g. ``solver/.../FGM.C:120-159``."""
        return f"{self.rel_path}:{self.start_line}-{self.end_line}"


def _window_lines(
    rel_path: str,
    kind: str,
    lines: list[str],
    base_offset: int,
    window: int,
    overlap: int,
) -> list[Chunk]:
    """Chunk a list of lines into overlapping windows.

    Args:
        base_offset: 0-based index of ``lines[0]`` within the original file,
            used to compute absolute line numbers.
    """
    chunks: list[Chunk] = []
    if not lines:
        return chunks
    step = max(1, window - overlap)
    for start in range(0, len(lines), step):
        block = lines[start : start + window]
        if not any(line.strip() for line in block):
            continue
        abs_start = base_offset + start
        chunks.append(
            Chunk(
                rel_path=rel_path,
                start_line=abs_start + 1,
                end_line=abs_start + len(block),
                kind=kind,
                text="\n".join(block),
            )
        )
        if start + window >= len(lines):
            break
    return chunks


def _chunk_markdown(
    sf: SourceFile, window: int, overlap: int
) -> list[Chunk]:
    """Split markdown on headings, windowing any oversized section."""
    lines = sf.text.splitlines()
    chunks: list[Chunk] = []
    section_start = 0  # 0-based index of the current section's first line

    def flush(start: int, end: int) -> None:
        # end is exclusive
        section = lines[start:end]
        if not any(line.strip() for line in section):
            return
        if len(section) <= window:
            chunks.append(
                Chunk(
                    rel_path=sf.rel_path,
                    start_line=start + 1,
                    end_line=end,
                    kind=sf.kind,
                    text="\n".join(section),
                )
            )
        else:
            chunks.extend(
                _window_lines(sf.rel_path, sf.kind, section, start, window, overlap)
            )

    for i, line in enumerate(lines):
        if line.lstrip().startswith("#") and i > section_start:
            flush(section_start, i)
            section_start = i
    flush(section_start, len(lines))
    return chunks


def chunk_file(
    sf: SourceFile, window: int = WINDOW, overlap: int = OVERLAP
) -> list[Chunk]:
    """Chunk one source file according to its kind."""
    if sf.kind == "markdown":
        return _chunk_markdown(sf, window, overlap)
    return _window_lines(
        sf.rel_path, sf.kind, sf.text.splitlines(), 0, window, overlap
    )


def chunk_corpus(
    files: list[SourceFile], window: int = WINDOW, overlap: int = OVERLAP
) -> list[Chunk]:
    """Chunk every file in the corpus."""
    chunks: list[Chunk] = []
    for sf in files:
        chunks.extend(chunk_file(sf, window, overlap))
    return chunks
