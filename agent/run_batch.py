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

"""Detached orchestrator: run several built cases with bounded concurrency.

Launched in the background by the run_cases tool. Runs each case's full solver
pipeline (via runner.run_case_tracked, so per-case run_status keeps working), at
most ``--max-parallel`` at a time, logging progress to ``--log``. The concurrency
limit is decided by the caller (run_cases) to fit the physical-core budget, so a
plain thread pool of that size never oversubscribes the machine.
"""

from __future__ import annotations

import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import of_writer
import runner


def _log(path: str, message: str) -> None:
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run built cases as a batch.")
    parser.add_argument("--max-parallel", type=int, default=1)
    parser.add_argument("--log", required=True)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()

    case_dirs = {name: of_writer.case_path(name) for name in args.cases}
    max_parallel = max(1, args.max_parallel)
    _log(
        args.log,
        f"Batch started: {len(args.cases)} case(s), {max_parallel} at a time "
        f"-> {', '.join(args.cases)}",
    )

    def run_one(name: str) -> tuple[str, str]:
        case_dir = case_dirs[name]
        if (case_dir / runner.CANCEL_MARKER).exists():
            _log(args.log, f"SKIP  {name}: cancelled before it started")
            return name, "cancelled"
        _log(args.log, f"START {name}")
        result = runner.run_case_tracked(case_dir)
        status = str(result.get("status"))
        _log(args.log, f"DONE  {name}: {status}")
        return name, status

    results: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=max_parallel) as pool:
        futures = {pool.submit(run_one, name): name for name in args.cases}
        for future in as_completed(futures):
            name = futures[future]
            try:
                _, status = future.result()
                results[name] = status
            except Exception as exc:  # noqa: BLE001 - record and continue
                results[name] = f"error: {exc}"
                _log(args.log, f"ERROR {name}: {exc}")

    summary = ", ".join(f"{name}={status}" for name, status in results.items())
    _log(args.log, f"Batch finished: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
