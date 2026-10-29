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

"""Tests for the evaluation harness itself (agent/evals/): the sandbox keeps real
cases and campaigns untouched, no solver can start, and the checks score what the
tools really did. Runs offline — no model is called.

    agent/.venv/bin/python agent/tests/test_evals.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "evals"))

import run_eval  # noqa: E402  (puts agent/ on sys.path itself)
import tools  # noqa: E402
from scenarios import SUITES  # noqa: E402

SLOT_ARGS = {"case_type": "laminarSlotBurner", "case_name": "evaltest",
             "inlet_velocity": 0.4, "T_inlet": 300, "T_plate": 400}

run_eval.install_sandbox()
RECORDER = run_eval.Recorder()


def test_scenario_ids_are_unique_and_every_turn_has_a_check() -> None:
    ids = [s["id"] for suite in SUITES.values() for s in suite]
    assert len(ids) == len(set(ids))
    checks = {"tool", "no_status", "no_tool", "case", "reply_all", "reply_none", "grounded"}
    for suite in SUITES.values():
        for s in suite:
            for turn in s["turns"]:
                assert checks & set(turn), s["id"]


def test_builds_land_in_the_sandbox_not_testcases() -> None:
    result = tools.execute_tool("build_case", SLOT_ARGS)
    assert result["status"] == "built", result
    case_dir = Path(result["case_dir"])
    assert case_dir.parent == run_eval.SANDBOX / "testcases"
    assert not (run_eval.REPO_ROOT_TESTCASES / case_dir.name).exists()


def test_run_tools_are_stubbed() -> None:
    for name in run_eval._STUBBED_RUN_TOOLS:
        result = tools.execute_tool(name, {"case_name": "evaltest"})
        assert "sandbox" in result.get("note", ""), (name, result)


def test_prepared_campaign_goes_to_the_sandbox() -> None:
    result = tools.execute_tool("run_campaign", {
        "case_type": "laminarSlotBurner", "variable": "T_plate", "values": [350, 450],
        "inlet_velocity": 0.4, "T_inlet": 300, "T_plate": 400,
        "n_procs_per_case": 1, "n_concurrent_cases": 1,
        "campaign_name": "evaltest", "launch": False})
    assert result["status"] == "prepared", result
    assert (run_eval.SANDBOX / "campaigns" / "evaltest").is_dir()


def test_checks_score_what_the_tools_did() -> None:
    log = RECORDER.new_turn()
    tools.execute_tool("build_case", SLOT_ARGS)
    built = log.calls[-1].result
    spec = {"tool": {"name": "build_case", "args": {"U_inlet": 0.4} | {"T_plate": 400},
                     "status": "built"},
            "case": {"T_plate": 400.0},
            "no_tool": ["run_case"],
            "reply_all": [r"(?i)ready"],
            "reply_none": [r"(?i)running now"]}
    checks = run_eval.evaluate_turn(spec, log, "The case is ready.", built)
    by = {c["check"]: c["ok"] for c in checks}
    assert by["case parameters"] and by["no run_case"]
    assert by["reply ~ (?i)ready"] and by["reply !~ (?i)running now"]
    # U_inlet is not an alias, so the tool check must fail on it.
    assert not by["tool build_case"]
    spec["tool"]["args"] = {"inlet_velocity": 0.404, "case_type": "laminarslotburner"}
    assert run_eval.evaluate_turn(spec, log, "", built)[0]["ok"]
    assert not run_eval.evaluate_turn({"case": {"T_plate": 500}}, log, "", built)[0]["ok"]


def test_lists_compare_as_sets_of_numbers() -> None:
    assert run_eval._value_matches([0.002, 0.003], ["0.003", 0.002])
    assert not run_eval._value_matches([0.002, 0.003], [2, 3])


def test_forced_reads_are_told_apart_from_model_reads() -> None:
    log = RECORDER.new_turn()
    tools.execute_tool("read_case_file", {"path": "solver/src/lookUp/FGM/FGMlib.c"})
    RECORDER.on_event("implementation question -> reading the top source file")
    RECORDER.on_tool_call("read_case_file", {"path": "solver/src/lookUp/FGM/FGMlib.c"})
    checks = run_eval.evaluate_turn({"grounded": "FGMlib"}, log, "", None)
    assert checks[0]["ok"] and not checks[1]["ok"]
    tools.execute_tool("read_case_file", {"path": "solver/src/lookUp/Table/FGMTable.C"})
    checks = run_eval.evaluate_turn({"grounded": "FGMlib"}, log, "", None)
    assert checks[1]["ok"] and "FGMTable" in checks[1]["detail"]


def test_guard_detections_are_counted_with_the_first_one_marked() -> None:
    def g(name, acted=False):
        return {"kind": "block", "name": name, "message": "", "acted": acted}
    results = [
        {"model": "m", "turns": [{"guards": [g("dont_start")]},
                                 {"guards": [g("case_type_not_named")]}]},
        {"model": "m", "turns": [{"guards": [g("case_type_not_named")]}]},
        {"model": "m", "turns": [{"guards": []}]},
    ]
    table = "\n".join(run_eval._guard_table(results, ["m"]))
    assert "| block · case_type_not_named | 2 (first in 1) |" in table
    assert "| block · dont_start | 1 (first in 1) |" in table
    assert run_eval._guard_table([{"model": "m", "turns": [{}]}], ["m"]) == []


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    try:
        for test in tests:
            test()
            print(f"  ok  {test.__name__}")
    finally:
        RECORDER.restore()
        run_eval.shutil.rmtree(run_eval.SANDBOX, ignore_errors=True)
    print(f"\n{len(tests)} passed")
