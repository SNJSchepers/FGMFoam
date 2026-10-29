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

"""Run the OpenFOAM solver pipeline for a prepared case.

Executes the standard sequence from the original bash ``run`` script -- for a
parallel case: ``blockMesh`` -> ``setFields`` -> ``decomposePar`` ->
``mpirun -n N FGMFoam -parallel`` -> ``reconstructPar`` -- each logged to
``log.<step>`` in the case directory. Steps run synchronously and stop at the
first failure. Requires the OpenFOAM environment and the built ``FGMFoam``
solver on PATH (i.e. run the agent from an OpenFOAM-sourced shell).
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import signal
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

import precheck

DEFAULT_TIMEOUT_S: int = 900
_TAIL_LINES: int = 15

# Marker files written in the case directory to track a background run.
ACTIVE_MARKER: str = ".run_active"
DONE_MARKER: str = ".run_done"
# Written by kill_jobs; a queued batch case with this marker is skipped.
CANCEL_MARKER: str = ".cancel"


def _pipeline_steps(parallel: bool, n_procs: int):
    """Return ``(pre_steps, solver_step, post_steps)`` for the run."""
    if parallel:
        pre = [
            ("blockMesh", ["blockMesh"]),
            ("setFields", ["setFields"]),
            ("decomposePar", ["decomposePar", "-force"]),
        ]
        solver = ("FGMFoam", ["mpirun", "-n", str(n_procs), "FGMFoam", "-parallel"])
        post = [("reconstructPar", ["reconstructPar"])]
    else:
        pre = [("blockMesh", ["blockMesh"]), ("setFields", ["setFields"])]
        solver = ("FGMFoam", ["FGMFoam"])
        post = []
    return pre, solver, post


def read_n_procs(case_dir: Path) -> int:
    """Read ``numberOfSubdomains`` from the case's decomposeParDict (default 1)."""
    try:
        text = (case_dir / "system" / "decomposeParDict").read_text()
    except OSError:
        return 1
    match = re.search(r"numberOfSubdomains\s+(\d+)", text)
    return int(match.group(1)) if match else 1


def pids_in_dir(root: Path) -> list[int]:
    """PIDs of processes whose working directory is at or under ``root``.

    Used to find a case's running solver (blockMesh/mpirun/FGMFoam/...), which
    all run with their cwd inside the case directory. This targets one case's
    processes only — it never touches the batch orchestrator or other cases.
    """
    try:
        root = root.resolve()
    except OSError:
        return []
    proc = Path("/proc")
    if not proc.is_dir():
        return []
    self_pid = os.getpid()
    pids: list[int] = []
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid == self_pid:
            continue
        try:
            cwd = (entry / "cwd").resolve()
        except (OSError, RuntimeError):
            continue
        if cwd == root or root in cwd.parents:
            pids.append(pid)
    return pids


def kill_case(case_dir: Path) -> list[int]:
    """Stop a case: mark it cancelled, kill any of its processes, fix markers.

    Returns the PIDs that were signalled. Killing is by individual PID (not the
    process group), so a case running inside a batch is stopped without killing
    the orchestrator or sibling cases.
    """
    (case_dir / CANCEL_MARKER).touch()
    pids = pids_in_dir(case_dir)
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    if pids:
        time.sleep(2)
        for pid in pids_in_dir(case_dir):
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
    # The case is no longer running.
    (case_dir / ACTIVE_MARKER).unlink(missing_ok=True)
    (case_dir / DONE_MARKER).touch()
    return pids


def physical_core_count(override: Optional[int] = None) -> int:
    """Best-effort count of physical CPU cores (not hyperthreads).

    Uses ``override`` if given, else unique (physical id, core id) pairs from
    /proc/cpuinfo, then ``lscpu`` (sockets x cores-per-socket), then the logical
    count. Always returns at least 1.
    """
    if override and override > 0:
        return override
    try:
        pairs: set[tuple[Optional[str], str]] = set()
        phys: Optional[str] = None
        with open("/proc/cpuinfo", encoding="utf-8") as handle:
            for line in handle:
                key, _, value = line.partition(":")
                key, value = key.strip(), value.strip()
                if key == "physical id":
                    phys = value
                elif key == "core id":
                    pairs.add((phys, value))
                elif not line.strip():
                    phys = None
        if pairs:
            return len(pairs)
    except OSError:
        pass
    try:
        out = subprocess.run(
            ["lscpu"], capture_output=True, text=True, timeout=5
        ).stdout
        cps = sockets = None
        for line in out.splitlines():
            if "Core(s) per socket" in line:
                cps = int(line.split(":")[1].strip())
            elif line.strip().startswith("Socket(s)"):
                sockets = int(line.split(":")[1].strip())
        if cps and sockets:
            return cps * sockets
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    return max(1, os.cpu_count() or 1)


def plan_batch(
    cases: list[tuple[str, int]],
    requested_parallel: int,
    physical_cores: int,
) -> dict[str, Any]:
    """Decide how a batch of cases may run within the physical core budget.

    Args:
        cases: ``(case_name, n_procs)`` for each case to run.
        requested_parallel: how many to run at once (clamped to >= 1).
        physical_cores: physical cores available.

    Returns a dict with ``runnable`` / ``skipped`` case lists, the
    ``max_parallel`` that will actually be used, and human-readable ``warnings``.
    """
    requested_parallel = max(1, requested_parallel)
    runnable = [(n, p) for n, p in cases if p <= physical_cores]
    skipped = [(n, p) for n, p in cases if p > physical_cores]

    warnings: list[str] = []
    for name, nprocs in skipped:
        warnings.append(
            f"'{name}' requests {nprocs} cores, but this machine has only "
            f"{physical_cores} physical cores — it cannot run and is skipped."
        )

    if runnable:
        max_nprocs = max(p for _, p in runnable)
        fit = max(1, physical_cores // max_nprocs)
    else:
        max_nprocs, fit = 0, 1
    effective = min(requested_parallel, fit)

    if runnable and requested_parallel > effective:
        tail = " (sequentially)." if effective == 1 else "."
        warnings.append(
            f"Requested {requested_parallel} simulations in parallel, but only "
            f"{effective} fit within {physical_cores} physical cores "
            f"({max_nprocs} cores each); running {effective} at a time{tail}"
        )

    return {
        "runnable": runnable,
        "skipped": skipped,
        "max_parallel": effective,
        "physical_cores": physical_cores,
        "warnings": warnings,
    }


def _tail(path: Path, n: int = _TAIL_LINES) -> str:
    try:
        return "\n".join(path.read_text(errors="replace").splitlines()[-n:])
    except OSError:
        return ""


def _run_step(
    case_dir: Path, name: str, argv: list[str], timeout: Optional[int]
) -> dict[str, Any]:
    """Run one command in the case dir, logging to ``log.<name>``.

    ``timeout`` is in seconds, or ``None`` for no limit (used for the solver and
    reconstructPar, which legitimately run for a long time).
    """
    log_path = case_dir / f"log.{name}"
    try:
        with log_path.open("w") as log:
            proc = subprocess.run(
                argv, cwd=case_dir, stdout=log, stderr=subprocess.STDOUT,
                timeout=timeout, check=False,
            )
    except FileNotFoundError:
        return {
            "step": name, "command": " ".join(argv), "returncode": 127,
            "error": f"command not found: {argv[0]} "
                     "(is the OpenFOAM environment sourced and FGMFoam built?)",
            "log_file": None, "tail": "",
        }
    except subprocess.TimeoutExpired:
        return {
            "step": name, "command": " ".join(argv), "returncode": None,
            "error": f"timed out after {timeout}s (still running / incomplete)",
            "log_file": log_path.name, "tail": _tail(log_path),
        }
    return {
        "step": name, "command": " ".join(argv), "returncode": proc.returncode,
        "log_file": log_path.name, "tail": _tail(log_path),
    }


def _read_end_time(case_dir: Path) -> Optional[float]:
    try:
        text = (case_dir / "system" / "controlDict").read_text()
    except OSError:
        return None
    match = re.search(r"\bendTime\s+([0-9eE.+-]+)\s*;", text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _time_dirs(directory: Path) -> list[float]:
    """Numeric (time) subdirectory names directly under ``directory`` (e.g. 0.1)."""
    times: list[float] = []
    if not directory.is_dir():
        return times
    for sub in directory.iterdir():
        if sub.is_dir():
            try:
                times.append(float(sub.name))
            except ValueError:
                pass
    return times


def _reached_end_time(case_dir: Path, parallel: bool) -> bool:
    """True if the latest written time directory reaches the case endTime."""
    end_time = _read_end_time(case_dir)
    if end_time is None:
        return False
    times = _time_dirs(case_dir / "processor0" if parallel else case_dir)
    return bool(times) and max(times) >= end_time - 1e-12


def run_pipeline(
    case_dir: Path,
    n_procs: Optional[int] = None,
    timeout: int = DEFAULT_TIMEOUT_S,
    skip_checks: bool = False,
) -> dict[str, Any]:
    """Run the full solver pipeline for ``case_dir``.

    Returns a dict with the overall ``status`` ("completed" | "failed" |
    "blocked" | "error") and the per-step results. Stops at the first real
    failure, and refuses to launch the solver when the pre-run checks report an
    error (see :mod:`precheck`) unless ``skip_checks`` is set.
    """
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case directory not found: {case_dir}"}
    if n_procs is None:
        n_procs = read_n_procs(case_dir)
    parallel = bool(n_procs and n_procs > 1)
    pre, solver, post = _pipeline_steps(parallel, n_procs)

    results: list[dict[str, Any]] = []

    def _fail(step: str) -> dict[str, Any]:
        return {"status": "failed", "failed_step": step,
                "n_procs": n_procs, "steps": results}

    for name, argv in pre:
        result = _run_step(case_dir, name, argv, timeout)
        results.append(result)
        if result["returncode"] != 0:
            return _fail(name)

    # Everything above takes seconds; the solver below can take hours. This is
    # the last cheap moment to notice that the case cannot produce a meaningful
    # result — and the mesh now exists, so the Courant check can run.
    gate = precheck.summarise(precheck.check_case(case_dir))
    precheck_result = gate
    if gate["status"] == "errors" and not skip_checks:
        return {
            "status": "blocked",
            "n_procs": n_procs,
            "steps": results,
            "precheck": gate,
            "message": (
                "Not launching the solver: the case would not produce a "
                "meaningful result. Fix the errors below, or pass "
                "skip_checks=true to run anyway."
            ),
        }

    # The solver legitimately runs for a long time, so it is NOT time-limited
    # (timeout=None); ``timeout`` guards only the quick setup steps above against a
    # hang. Older FGMFoam builds abort during shutdown even after a successful run,
    # so a non-zero exit only counts as a failure if endTime was not reached.
    result = _run_step(case_dir, solver[0], solver[1], None)
    results.append(result)
    if result["returncode"] != 0:
        if _reached_end_time(case_dir, parallel):
            result["note"] = (
                f"FGMFoam exited with code {result['returncode']} but reached "
                "endTime; it aborts during shutdown, so the run is treated as "
                "complete."
            )
        else:
            return _fail("FGMFoam")

    # reconstructPar finishes the run (and can be slow for large cases), so it is
    # likewise not time-limited.
    for name, argv in post:
        result = _run_step(case_dir, name, argv, None)
        results.append(result)
        if result["returncode"] != 0:
            return _fail(name)

    # Clean up per-processor directories after a successful reconstruct.
    for proc_dir in case_dir.glob("processor*"):
        shutil.rmtree(proc_dir, ignore_errors=True)

    return {"status": "completed", "n_procs": n_procs, "steps": results,
            "precheck": precheck_result}


def run_case_tracked(case_dir: Path, timeout: int = DEFAULT_TIMEOUT_S) -> dict[str, Any]:
    """Run one case synchronously with the run markers set, so ``run_status``
    works for cases executed inside a batch."""
    (case_dir / DONE_MARKER).unlink(missing_ok=True)
    (case_dir / ACTIVE_MARKER).touch()
    try:
        return run_pipeline(case_dir, timeout=timeout)
    finally:
        (case_dir / ACTIVE_MARKER).unlink(missing_ok=True)
        (case_dir / DONE_MARKER).touch()


def _launch_background(case_dir: Path, solver, post) -> subprocess.Popen:
    """Launch the solver + post steps + cleanup detached, setting the run markers.

    Runs ``solver`` then each ``post`` step regardless of FGMFoam's exit code (it
    aborts on shutdown), then removes the processor*/ dirs, clears the active
    marker and writes the done marker. Returns the detached ``Popen``.
    """
    commands = [f"{shlex.join(solver[1])} > log.{solver[0]} 2>&1"]
    commands += [f"{shlex.join(argv)} > log.{name} 2>&1" for name, argv in post]
    commands += ["rm -rf processor*", f"rm -f {ACTIVE_MARKER}", f"touch {DONE_MARKER}"]
    script = "; ".join(commands)

    (case_dir / DONE_MARKER).unlink(missing_ok=True)
    (case_dir / ACTIVE_MARKER).touch()
    return subprocess.Popen(
        ["bash", "-c", script], cwd=case_dir,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        start_new_session=True,
    )


def run_pipeline_background(
    case_dir: Path,
    n_procs: Optional[int] = None,
    timeout: int = DEFAULT_TIMEOUT_S,
    skip_checks: bool = False,
) -> dict[str, Any]:
    """Run mesh/setup synchronously, then launch the solver in the background.

    blockMesh/setFields/decomposePar run synchronously (fast, so setup errors are
    reported at once). The solver, reconstructPar, and cleanup then run detached
    in the background, writing to ``log.FGMFoam`` etc. Returns immediately with
    the PID and the log to tail. Use :func:`run_status` to check progress.
    """
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case directory not found: {case_dir}"}
    if n_procs is None:
        n_procs = read_n_procs(case_dir)
    parallel = bool(n_procs and n_procs > 1)
    pre, solver, post = _pipeline_steps(parallel, n_procs)

    results: list[dict[str, Any]] = []
    for name, argv in pre:
        result = _run_step(case_dir, name, argv, timeout)
        results.append(result)
        if result["returncode"] != 0:
            return {"status": "failed", "failed_step": name,
                    "n_procs": n_procs, "steps": results}

    # Same gate as the synchronous path: the setup steps above are cheap, the
    # solver is not, and the mesh now exists so every check can run.
    gate = precheck.summarise(precheck.check_case(case_dir))
    if gate["status"] == "errors" and not skip_checks:
        return {
            "status": "blocked",
            "n_procs": n_procs,
            "steps": results,
            "precheck": gate,
            "message": (
                "Not launching the solver: the case would not produce a "
                "meaningful result. Fix the errors below, or pass "
                "skip_checks=true to run anyway."
            ),
        }

    proc = _launch_background(case_dir, solver, post)

    log_path = case_dir / f"log.{solver[0]}"
    return {
        "status": "launched",
        "precheck": gate,
        "pid": proc.pid,
        "n_procs": n_procs,
        "case_dir": str(case_dir),
        "log_file": log_path.name,
        "monitor": f"tail -f {log_path}",
        "prep_steps": [r["step"] for r in results],
        "message": (
            "Solver launched in the background. Monitor with the tail command; "
            "reconstructPar runs automatically when it finishes."
        ),
    }


def plan_restart(case_dir: Path, n_procs: Optional[int] = None) -> dict[str, Any]:
    """Decide whether/how a case can be restarted from its latest time (no I/O side effects).

    Returns a dict whose ``status`` is one of:
      - ``"error"`` — cannot restart (not meshed/run, or no time > 0 to continue from);
      - ``"nothing_to_do"`` — the latest time already reaches endTime (raise endTime first);
      - ``"ready"`` — restartable, with ``latest_time``, ``end_time``, ``n_procs``,
        ``parallel``, ``decomposed`` (processor*/ present) and ``need_decompose``
        (True when the latest reconstructed time must be decomposed because the
        decomposition is missing or older than that time).
    """
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case directory not found: {case_dir}"}
    if not (case_dir / "constant" / "polyMesh").is_dir():
        return {
            "status": "error",
            "error": "the case has no mesh (constant/polyMesh) yet — it has not been "
            "run; use run_case for the first run.",
        }
    if n_procs is None:
        n_procs = read_n_procs(case_dir)
    parallel = bool(n_procs and n_procs > 1)

    root_latest = max(_time_dirs(case_dir), default=None)
    proc_latest = max(_time_dirs(case_dir / "processor0"), default=None)
    decomposed = bool(list(case_dir.glob("processor*")))

    available = [t for t in (root_latest, proc_latest) if t is not None]
    latest = max(available) if available else None
    if latest is None or latest <= 1e-12:
        return {
            "status": "error",
            "error": "no completed time step to restart from (only time 0 present); "
            "run the case first with run_case.",
        }

    end_time = _read_end_time(case_dir)
    if end_time is not None and latest >= end_time - 1e-12:
        return {
            "status": "nothing_to_do",
            "latest_time": latest,
            "end_time": end_time,
            "message": (
                f"The latest available time ({latest}) already reaches endTime "
                f"({end_time}); increase endTime with change_settings before restarting."
            ),
        }

    # The decomposition is reusable only if the processors hold the latest time.
    proc_current = proc_latest is not None and abs(proc_latest - latest) <= 1e-12
    need_decompose = parallel and not proc_current
    return {
        "status": "ready",
        "n_procs": n_procs,
        "parallel": parallel,
        "latest_time": latest,
        "end_time": end_time,
        "decomposed": decomposed,
        "need_decompose": need_decompose,
    }


def restart_pipeline_background(
    case_dir: Path, timeout: int = DEFAULT_TIMEOUT_S
) -> dict[str, Any]:
    """Restart a case from its latest saved time, in the background.

    Reuses the existing mesh (no blockMesh/setFields). For a parallel case it runs
    ``decomposePar -force -latestTime`` only when :func:`plan_restart` reports the
    decomposition is missing or stale (e.g. the case was reconstructed); otherwise
    it reuses the existing processor*/ data. Then it launches FGMFoam (and
    reconstructPar) detached, exactly like :func:`run_pipeline_background`.

    The caller must have set controlDict ``startFrom latestTime`` (tools.py does).
    """
    plan = plan_restart(case_dir)
    if plan["status"] != "ready":
        return plan
    n_procs, parallel = plan["n_procs"], plan["parallel"]

    results: list[dict[str, Any]] = []
    if plan["need_decompose"]:
        result = _run_step(
            case_dir, "decomposePar", ["decomposePar", "-force", "-latestTime"], timeout
        )
        results.append(result)
        if result["returncode"] != 0:
            return {"status": "failed", "failed_step": "decomposePar",
                    "n_procs": n_procs, "steps": results}

    if parallel:
        solver = ("FGMFoam", ["mpirun", "-n", str(n_procs), "FGMFoam", "-parallel"])
        post = [("reconstructPar", ["reconstructPar"])]
    else:
        solver = ("FGMFoam", ["FGMFoam"])
        post = []

    proc = _launch_background(case_dir, solver, post)
    log_path = case_dir / f"log.{solver[0]}"
    return {
        "status": "launched",
        "pid": proc.pid,
        "n_procs": n_procs,
        "case_dir": str(case_dir),
        "restart_time": plan["latest_time"],
        "end_time": plan["end_time"],
        "decomposed_latest_time": plan["need_decompose"],
        "log_file": log_path.name,
        "monitor": f"tail -f {log_path}",
        "message": (
            f"Restarted from time {plan['latest_time']} toward endTime "
            f"{plan['end_time']} in the background"
            + (" (decomposed the latest time first)"
               if plan["need_decompose"] else " (reused the existing decomposition)")
            + ". Monitor with the tail command; reconstructPar runs when it finishes."
        ),
    }


def run_status(case_dir: Path) -> dict[str, Any]:
    """Report the state of a (possibly background) run for ``case_dir``."""
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case directory not found: {case_dir}"}

    active = (case_dir / ACTIVE_MARKER).exists()
    done = (case_dir / DONE_MARKER).exists()
    state = "running" if active else ("finished" if done else "not started")

    # Latest written time: processor0 while a parallel run is live, otherwise the
    # (reconstructed) case directory.
    search = case_dir / "processor0"
    if not search.is_dir():
        search = case_dir
    times = _time_dirs(search)
    latest = max(times) if times else None
    end_time = _read_end_time(case_dir)
    reached = (
        latest is not None and end_time is not None and latest >= end_time - 1e-12
    )

    log_path = case_dir / "log.FGMFoam"
    result: dict[str, Any] = {
        "state": state,
        "latest_time": latest,
        "end_time": end_time,
        "reached_end_time": reached,
        "log_file": "log.FGMFoam",
        "log_tail": _tail(log_path) if log_path.exists() else "",
    }
    if state == "finished":
        if reached:
            result["outcome"] = "succeeded"
            result["note"] = (
                "The run reached endTime and completed. FGMFoam aborts during "
                "shutdown even on success, so a SIGABRT / exit code 134 at the "
                "end of the log is expected and does NOT indicate failure."
            )
        else:
            result["outcome"] = "incomplete"
            result["note"] = (
                "The run finished but did not reach endTime; check log.FGMFoam "
                "for the cause."
            )
    return result
