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

"""CLI to build (or rebuild) the RAG index over the solver code and docs.

Usage (from the repository root):

    agent/.venv/bin/python agent/index_build.py            # build if missing/stale
    agent/.venv/bin/python agent/index_build.py --force    # always rebuild
    agent/.venv/bin/python agent/index_build.py --status   # report status (no API key)
    agent/.venv/bin/python agent/index_build.py --embed-model mistral-embed
    agent/.venv/bin/python agent/index_build.py --model local   # use a preset's embedder
"""

from __future__ import annotations

import argparse
import sys

import config as config_module
from config import INDEX_DIR, ConfigError
from rag.index import build_index, index_exists, staleness


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the FGMFoam agent RAG index.")
    parser.add_argument("--force", action="store_true", help="rebuild even if up to date")
    parser.add_argument("--status", action="store_true",
                        help="report whether the index is up to date and exit (no API key needed)")
    parser.add_argument("--embed-model", default=None, help="override FGM_EMBED_MODEL")
    parser.add_argument("--model", default=None,
                        help="model preset to take the embedder from (e.g. claude, local)")
    args = parser.parse_args()

    # Status check needs no API key, so resolve config without requiring one.
    if args.status:
        cfg = config_module.load_config(
            require_key=False, embed_model=args.embed_model, model=args.model
        )
        reasons = staleness(cfg.embed_model)
        print("Index status:", "; ".join(reasons) if reasons else f"up to date ({INDEX_DIR})")
        return 0

    try:
        cfg = config_module.load_config(embed_model=args.embed_model, model=args.model)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 1

    if index_exists() and not args.force:
        reasons = staleness(cfg.embed_model)
        if not reasons:
            print(f"Index is up to date ({INDEX_DIR}). Use --force to rebuild.")
            return 0
        print("Rebuilding index because:", "; ".join(reasons))

    print(f"Building index with embedding model '{cfg.embed_model}' ...")
    index = build_index(cfg)
    index.save()

    m = index.manifest
    print(
        f"Done. Indexed {m['n_files']} files into {m['n_chunks']} chunks "
        f"(dim {m['dim']}). Saved to {INDEX_DIR}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
