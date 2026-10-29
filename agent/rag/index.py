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

"""Build, persist, load, and search the local vector index.

The corpus is tiny, so the index is a flat NumPy matrix searched by cosine
similarity -- no vector database needed. The index directory holds:

    embeddings.npy   float32 matrix (n_chunks, dim), L2-normalised
    chunks.jsonl     one chunk (text + metadata) per line, row-aligned
    manifest.json    embed model, timestamp, and per-source content hashes

The manifest lets the agent detect when the indexed sources have changed and
warn that a rebuild is needed.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from config import INDEX_DIR, Config
from rag.chunker import Chunk, chunk_corpus
from rag.corpus import collect_corpus
from rag.embeddings import embed_query, embed_texts

EMBEDDINGS_FILE = "embeddings.npy"
CHUNKS_FILE = "chunks.jsonl"
MANIFEST_FILE = "manifest.json"

# -- hybrid (dense + lexical) retrieval ------------------------------------- #
#
# Small embedding models (e.g. nomic-embed-text) miss exact-term matches that
# live in identifiers and file paths: asking about the "enthalpy boundary
# condition" should surface lookUpEnthalpyFvPatchScalarField.C, but dense
# similarity ranks prose docs and unrelated .H files above it (the right source
# can sit past rank 100). The fix is a lexical bonus added to the cosine score:
# query terms appearing in a chunk's FILE PATH (where class/file names live) are
# the strong, discriminative signal, so they are weighted well above terms found
# only in the chunk text. The corpus is tiny, so scoring every chunk per query
# is cheap. Tune via the two weights below; with 0 lexical weight this reduces to
# pure cosine.
_TOKEN_RE = re.compile(r"[a-z0-9]+")
_LEXICAL_PATH_WEIGHT: float = 0.35
_LEXICAL_TEXT_WEIGHT: float = 0.05
# A query word found in at least this share of chunk paths earns no path bonus.
_LEXICAL_PATH_MAX_SHARE: float = 0.3
# Filler / intent words that carry no topic, dropped from the lexical query so
# they neither dilute the score nor match paths by accident.
_LEXICAL_STOPWORDS: frozenset[str] = frozenset({
    "the", "and", "for", "with", "this", "that", "these", "those", "from",
    "how", "does", "doe", "what", "which", "why", "when", "where", "who",
    "are", "is", "was", "were", "been", "being", "use", "used", "using",
    "can", "could", "would", "should", "you", "your", "explain", "describe",
    "tell", "about", "work", "works", "working", "case", "cases", "please",
    "set", "give", "show", "into", "its", "their", "have", "has", "had",
})


def _lexical_query_terms(query: str) -> list[str]:
    """Distinct content tokens from the query for lexical matching."""
    terms: list[str] = []
    seen: set[str] = set()
    for tok in _TOKEN_RE.findall(query.lower()):
        if len(tok) < 3 or tok in _LEXICAL_STOPWORDS or tok in seen:
            continue
        seen.add(tok)
        terms.append(tok)
    return terms


def _lexical_scores(query: str, chunks: list[Chunk]) -> np.ndarray:
    """Per-chunk lexical bonus: query terms found in the path (weighted high)
    plus terms found in the chunk text (weighted low), each as a fraction of the
    query's content terms. Returns zeros when the query has no content terms."""
    terms = _lexical_query_terms(query)
    scores = np.zeros(len(chunks), dtype=np.float32)
    if not terms:
        return scores
    n = len(terms)
    paths = [chunk.rel_path.lower() for chunk in chunks]
    # A word in most paths names no file: "solver" is in every solver/ path (and
    # solver_notes.md), so "what is the solver algorithm?" lifted all source code
    # above the knowledge note that answers it. Only discriminating words earn
    # the path bonus; they still count in the text bonus below.
    path_terms = [t for t in terms
                  if sum(t in path for path in paths) < _LEXICAL_PATH_MAX_SHARE * len(paths)]
    for i, chunk in enumerate(chunks):
        path = paths[i]
        path_hits = sum(1 for t in path_terms if t in path)
        text_tokens = set(_TOKEN_RE.findall(chunk.text.lower()))
        text_hits = sum(1 for t in terms if t in text_tokens)
        scores[i] = (
            _LEXICAL_PATH_WEIGHT * (path_hits / n)
            + _LEXICAL_TEXT_WEIGHT * (text_hits / n)
        )
    return scores


@dataclass
class Hit:
    """A retrieved chunk and its similarity score."""

    chunk: Chunk
    score: float


def _normalise(matrix: np.ndarray) -> np.ndarray:
    """L2-normalise rows so a dot product equals cosine similarity."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return matrix / norms


class VectorIndex:
    """An in-memory, file-backed cosine-similarity index."""

    def __init__(
        self,
        embeddings: np.ndarray,
        chunks: list[Chunk],
        manifest: dict,
    ) -> None:
        self.embeddings = embeddings  # assumed already normalised
        self.chunks = chunks
        self.manifest = manifest

    # -- search ------------------------------------------------------------- #
    def _dense_scores(self, query_vec: np.ndarray) -> np.ndarray:
        """Cosine similarity of every chunk to ``query_vec`` (rows are normalised)."""
        q = query_vec.astype(np.float32)
        norm = np.linalg.norm(q)
        if norm:
            q = q / norm
        return self.embeddings @ q

    @staticmethod
    def _top_hits(scores: np.ndarray, chunks: list[Chunk], k: int) -> list["Hit"]:
        k = min(k, scores.shape[0])
        top = np.argpartition(-scores, k - 1)[:k]
        top = top[np.argsort(-scores[top])]
        return [Hit(chunk=chunks[i], score=float(scores[i])) for i in top]

    def search(self, query_vec: np.ndarray, k: int) -> list[Hit]:
        """Return the top-k chunks most similar to ``query_vec`` (dense only)."""
        if self.embeddings.shape[0] == 0:
            return []
        return self._top_hits(self._dense_scores(query_vec), self.chunks, k)

    def search_hybrid(
        self, query_vec: np.ndarray, query_text: str, k: int
    ) -> list[Hit]:
        """Top-k by cosine similarity plus a lexical path/text bonus.

        The lexical term lifts chunks whose file path or text contains the query's
        content words — recovering exact-name matches (class/file identifiers) that
        a small embedding model ranks too low. ``Hit.score`` is the combined score
        that determined the ranking. See the lexical helpers above.
        """
        if self.embeddings.shape[0] == 0:
            return []
        scores = self._dense_scores(query_vec) + _lexical_scores(
            query_text, self.chunks
        )
        return self._top_hits(scores, self.chunks, k)

    # -- persistence -------------------------------------------------------- #
    def save(self, index_dir: Path = INDEX_DIR) -> None:
        index_dir.mkdir(parents=True, exist_ok=True)
        np.save(index_dir / EMBEDDINGS_FILE, self.embeddings)
        with (index_dir / CHUNKS_FILE).open("w", encoding="utf-8") as fh:
            for chunk in self.chunks:
                fh.write(json.dumps(asdict(chunk)) + "\n")
        (index_dir / MANIFEST_FILE).write_text(
            json.dumps(self.manifest, indent=2), encoding="utf-8"
        )

    @classmethod
    def load(cls, index_dir: Path = INDEX_DIR) -> "VectorIndex":
        embeddings = np.load(index_dir / EMBEDDINGS_FILE)
        chunks: list[Chunk] = []
        with (index_dir / CHUNKS_FILE).open(encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    chunks.append(Chunk(**json.loads(line)))
        manifest = json.loads((index_dir / MANIFEST_FILE).read_text(encoding="utf-8"))
        return cls(embeddings=embeddings, chunks=chunks, manifest=manifest)


def build_index(config: Config) -> VectorIndex:
    """Collect the corpus, chunk it, embed it, and return a normalised index."""
    files = collect_corpus()
    chunks = chunk_corpus(files)
    matrix = embed_texts([c.text for c in chunks], config=config)
    manifest = {
        # The embedding side is recorded (not the chat provider) because that is
        # what the stored vectors depend on: config reuses it so switching chat
        # models never forces a rebuild. "provider" is kept for older readers.
        "embed_model": config.embed_model,
        "embed_provider": config.embed_provider,
        "embed_base_url": config.embed_base_url or "",
        "provider": config.embed_provider,
        "chat_provider": config.provider,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_files": len(files),
        "n_chunks": len(chunks),
        "dim": int(matrix.shape[1]) if matrix.size else 0,
        "sources": {sf.rel_path: sf.sha256 for sf in files},
    }
    return VectorIndex(
        embeddings=_normalise(matrix) if matrix.size else matrix,
        chunks=chunks,
        manifest=manifest,
    )


def index_exists(index_dir: Path = INDEX_DIR) -> bool:
    return all(
        (index_dir / name).is_file()
        for name in (EMBEDDINGS_FILE, CHUNKS_FILE, MANIFEST_FILE)
    )


def staleness(embed_model: str, index_dir: Path = INDEX_DIR) -> list[str]:
    """Return human-readable reasons the on-disk index is out of date (if any).

    Only needs the embedding model name and the corpus hashes, so it requires no
    API key.
    """
    if not index_exists(index_dir):
        return ["no index found"]
    manifest = json.loads((index_dir / MANIFEST_FILE).read_text(encoding="utf-8"))
    reasons: list[str] = []
    if manifest.get("embed_model") != embed_model:
        reasons.append(
            f"embedding model changed "
            f"({manifest.get('embed_model')} -> {embed_model})"
        )
    current = {sf.rel_path: sf.sha256 for sf in collect_corpus()}
    stored = manifest.get("sources", {})
    if set(current) != set(stored):
        reasons.append("corpus file set changed")
    elif any(current[p] != stored.get(p) for p in current):
        reasons.append("corpus file contents changed")
    return reasons


def retrieve(config: Config, query: str, index: VectorIndex) -> list[Hit]:
    """Embed a query and return the top-k hits (hybrid dense + lexical) from ``index``."""
    query_vec = np.asarray(embed_query(query, config=config), dtype=np.float32)
    return index.search_hybrid(query_vec, query, config.top_k)
