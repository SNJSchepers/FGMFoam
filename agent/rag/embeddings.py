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

"""Embed text via the configured EMBEDDING provider, with simple
exponential-backoff retries for transient errors.

Used both at index-build time (the whole corpus) and per query at chat time (one
string). The call goes through ``llm.make_embed_client``, which resolves the embedding
provider independently of the chat provider (``config.embed_provider``) — the
Anthropic API has no embeddings endpoint, so a Claude chat model embeds via
Mistral or a local OpenAI-compatible server instead.
"""

from __future__ import annotations

import time

import numpy as np

from config import Config
from llm import make_embed_client

# Keep batches modest (some endpoints reject very large ones).
BATCH_SIZE: int = 64
MAX_RETRIES: int = 4
BACKOFF_BASE: float = 1.5  # seconds


def _batches(items: list[str], size: int):
    for start in range(0, len(items), size):
        yield items[start : start + size]


def _embed_batch(client, model: str, batch: list[str]) -> list[list[float]]:
    """Embed one batch, retrying transient failures with backoff."""
    last_error: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            return client.embed(texts=batch, model=model)
        except Exception as exc:  # noqa: BLE001 - retry any transient API error
            last_error = exc
            if attempt < MAX_RETRIES - 1:
                time.sleep(BACKOFF_BASE ** attempt)
    raise RuntimeError(
        f"Embedding request failed after {MAX_RETRIES} attempts: {last_error}"
    )


def embed_texts(
    texts: list[str], *, config: Config, batch_size: int = BATCH_SIZE
) -> np.ndarray:
    """Embed a list of texts into a float32 matrix of shape (len(texts), dim)."""
    if not texts:
        return np.zeros((0, 0), dtype=np.float32)

    client = make_embed_client(config)
    vectors: list[list[float]] = []
    for batch in _batches(texts, batch_size):
        vectors.extend(_embed_batch(client, config.embed_model, batch))
    return np.asarray(vectors, dtype=np.float32)


def embed_query(query: str, *, config: Config) -> np.ndarray:
    """Embed a single query string into a 1-D float32 vector."""
    return embed_texts([query], config=config)[0]
