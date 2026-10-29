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

"""Tests for campaigns: what may be swept, how the machine is shared, and the
results table and graph. Runs offline and starts no solver.

    agent/.venv/bin/python agent/tests/test_campaign.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import campaign as campaign_module  # noqa: E402
import tools  # noqa: E402


def test_equivalence_ratio_and_pressure_are_refused_with_the_reason() -> None:
    """Both are baked into the FGM table, so a sweep is impossible, not merely
    unsupported — and the message has to say why."""
    for name in ("equivalence_ratio", "p_0"):
        reason = campaign_module.check_variable("laminarFlamelet", name)
        assert reason and "table" in reason


def test_the_slot_burner_geometry_can_be_swept() -> None:
    """Its blockMeshDict is parameterised, so each case gets its own mesh."""
    for name in ("plate_thickness", "hole_pitch", "hole_size", "outlet_length",
                 "flame_thickness"):
        assert campaign_module.check_variable("laminarSlotBurner", name) is None
        assert name in campaign_module.GEOMETRIC


def test_unknown_geometry_names_point_at_the_supported_ones() -> None:
    for name in ("slot_width", "nozzle_diameter", "domain_height"):
        reason = campaign_module.check_variable("laminarSlotBurner", name)
        assert reason and "plate_thickness" in reason


def test_flamelet_geometry_is_still_refused() -> None:
    """Only the slot burner's mesh is parameterised."""
    for name in ("hole_size", "plate_thickness", "slot_width"):
        reason = campaign_module.check_variable("laminarFlamelet", name)
        assert reason
        assert "not parameterised" in reason or "not a sweepable" in reason


def test_sweepable_variables_are_per_case_type() -> None:
    assert campaign_module.check_variable("laminarFlamelet", "T_inlet") is None
    assert campaign_module.check_variable("laminarSlotBurner", "T_plate") is None
    assert campaign_module.check_variable("laminarSlotBurner", "inlet_velocity") is None
    # The flamelet has no plate, so sweeping its temperature is meaningless.
    reason = campaign_module.check_variable("laminarFlamelet", "T_plate")
    assert reason and "not a sweepable parameter" in reason


def test_cpu_options_never_oversubscribe_physical_cores() -> None:
    for cores in (2, 4, 8, 16):
        for n_cases in (2, 3, 5, 12):
            for option in campaign_module.cpu_options(n_cases, cores):
                assert option["cores_used"] <= cores
                assert option["n_concurrent_cases"] <= n_cases
                assert option["waves"] >= 1


def test_cpu_options_span_sequential_to_fully_parallel() -> None:
    options = campaign_module.cpu_options(4, 8)
    assert options[0]["n_procs_per_case"] == 8
    assert options[0]["n_concurrent_cases"] == 1      # one case, all cores
    assert options[-1]["n_procs_per_case"] == 1       # each case serial
    assert options[-1]["n_concurrent_cases"] == 4


def test_a_plan_that_exceeds_the_machine_is_rejected() -> None:
    reason = campaign_module.validate_cpu_plan(4, 4, 8)
    assert reason and "only 8 physical cores" in reason
    assert campaign_module.validate_cpu_plan(4, 2, 8) is None
    assert campaign_module.validate_cpu_plan(0, 2, 8) is not None


def test_the_cpu_plan_must_come_from_the_user() -> None:
    """Without a plan the tool must ask, not choose."""
    result = tools.execute_tool("run_campaign", {
        "case_type": "laminarFlamelet", "variable": "T_inlet",
        "values": [300, 400, 500],
    })
    assert result["status"] == "missing_cpu_plan"
    assert result["missing"] == ["n_procs_per_case", "n_concurrent_cases"]
    assert result["physical_cores"] >= 1
    assert result["options"]
    assert "do not choose for them" in result["message"].lower()


def test_a_campaign_needs_at_least_two_values() -> None:
    result = tools.execute_tool("run_campaign", {
        "case_type": "laminarFlamelet", "variable": "T_inlet", "values": [300],
        "n_procs_per_case": 1, "n_concurrent_cases": 1,
    })
    assert result["status"] == "error"
    assert "at least two" in result["error"]


def test_refused_variables_never_reach_the_cpu_question() -> None:
    result = tools.execute_tool("run_campaign", {
        "case_type": "laminarFlamelet", "variable": "equivalence_ratio",
        "values": [0.6, 0.7, 0.8],
    })
    assert result["status"] == "error"
    assert "table" in result["error"]


def test_plan_makes_one_case_per_value_with_the_chosen_core_count() -> None:
    plan = campaign_module.plan(
        "demo", "laminarFlamelet", "T_inlet", [300.0, 450.0],
        {"sim_time": 0.01}, n_procs_per_case=2, n_concurrent_cases=4,
    )
    assert len(plan.cases) == 2
    assert [c["value"] for c in plan.cases] == [300.0, 450.0]
    for entry in plan.cases:
        assert entry["case_name"].startswith("laminarFlamelet_demo_T_inlet_")
        assert entry["params"]["n_procs"] == 2
        assert entry["params"]["sim_time"] == 0.01
    # The swept variable is set per case, not inherited from the base.
    assert plan.cases[0]["params"]["T_inlet"] == 300.0


def test_slug_is_filesystem_safe() -> None:
    assert campaign_module._slug(300.0) == "300"
    assert campaign_module._slug(0.25) == "0p25"
    assert campaign_module._slug(-5.0) == "m5"


def _fake_campaign(tmp: Path, rows) -> tuple:
    """A campaign whose results are supplied directly, to test table and plot."""
    campaign_module.CAMPAIGN_DIR = tmp
    plan = campaign_module.plan(
        "synthetic", "laminarFlamelet", "T_inlet",
        [r[0] for r in rows], {}, 1, 1,
    )
    results = {
        "campaign": plan.name, "case_type": plan.case_type,
        "variable": "T_inlet", "metric": "laminar_flame_speed",
        "n_cases": len(rows), "n_measured": sum(1 for r in rows if r[1] is not None),
        "rows": [
            {"T_inlet": value, "case_name": f"c{i}",
             "status": "measured" if speed is not None else "blow_off",
             "laminar_flame_speed": speed}
            for i, (value, speed) in enumerate(rows)
        ],
    }
    return plan, results


def test_results_are_written_as_csv() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        plan, results = _fake_campaign(
            Path(tmp), [(300.0, 0.19), (400.0, 0.35), (500.0, 0.58)]
        )
        path = campaign_module.write_csv(plan, results)
        text = path.read_text()
    assert "T_inlet" in text and "laminar_flame_speed" in text
    assert "0.19" in text and "0.58" in text


def test_the_graph_is_produced() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        plan, results = _fake_campaign(
            Path(tmp), [(300.0, 0.19), (400.0, 0.35), (500.0, 0.58)]
        )
        path, error = campaign_module.plot(plan, results)
        assert error is None, error
        assert path is not None and path.is_file()
        assert path.stat().st_size > 5000, "a real PNG, not an empty canvas"


def test_failed_cases_are_kept_in_the_table_and_marked_on_the_graph() -> None:
    """A blow-off is a result, not a missing value."""
    with tempfile.TemporaryDirectory() as tmp:
        plan, results = _fake_campaign(
            Path(tmp), [(300.0, 0.19), (400.0, None), (500.0, 0.58)]
        )
        text = campaign_module.write_csv(plan, results).read_text()
        assert "blow_off" in text
        path, error = campaign_module.plot(plan, results)
        assert error is None and path.is_file()


def test_too_few_measurements_reports_instead_of_plotting_nonsense() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        plan, results = _fake_campaign(Path(tmp), [(300.0, 0.19), (400.0, None)])
        path, error = campaign_module.plot(plan, results)
    assert path is None
    assert "not enough to plot" in error


def test_domain_length_matches_the_mesh_dictionary() -> None:
    """The plugin's domain length must agree with l_in + l_plate + l_out."""
    import re
    import case_templates
    text = (case_templates.TESTCASES_DIR / "laminarSlotBurner" / "system.orig"
            / "blockMeshDict").read_text()
    l_in = float(re.search(r"^l_in\s+([0-9.eE+-]+);", text, re.M).group(1))
    assert re.search(r"^l_out\s+OUTLETLENGTH;", text, re.M)
    # The run script's value is the shipped default, so both routes agree.
    run = (case_templates.TESTCASES_DIR / "laminarSlotBurner" / "run").read_text()
    l_out = float(re.search(r"^outlet_length=([0-9.eE+-]+)", run, re.M).group(1))
    assert abs(case_templates.domain_length(
        "laminarSlotBurner", {"plate_thickness": 3.1139e-3})
        - (l_in + 3.1139e-3 + l_out)) < 1e-12
    assert abs(case_templates.domain_length(
        "laminarSlotBurner", {"plate_thickness": 3.1139e-3, "outlet_length": 0.035})
        - (l_in + 3.1139e-3 + 0.035)) < 1e-12


def test_flow_throughs_set_sim_time_per_case() -> None:
    """A shared sim_time would over-run fast cases and under-run slow ones."""
    plan = campaign_module.plan(
        "ft", "laminarSlotBurner", "inlet_velocity", [0.2, 0.4],
        {"T_inlet": 300.0, "T_plate": 300.0}, 1, 1, flow_throughs=2,
    )
    times = [c["params"]["sim_time"] for c in plan.cases]
    assert abs(times[0] - 2 * 0.0280222 / 0.2) < 1e-3
    assert abs(times[0] / times[1] - 2.0) < 0.01


def test_a_campaign_can_hold_a_longer_outlet_fixed() -> None:
    """Fast flames outgrow the default domain; the sweep keeps one outlet length
    for every case and the flow-through time follows it."""
    import tools
    props = next(t for t in tools.TOOLS if t["function"]["name"] == "run_campaign")
    assert "outlet_length" in props["function"]["parameters"]["properties"]
    plan = campaign_module.plan(
        "ft", "laminarSlotBurner", "inlet_velocity", [0.2, 0.4],
        {"T_inlet": 300.0, "T_plate": 300.0, "outlet_length": 0.035}, 1, 1,
        flow_throughs=2,
    )
    assert all(c["params"]["outlet_length"] == 0.035 for c in plan.cases)
    length = 6.2278e-3 + 3.1139e-3 + 0.035
    assert abs(plan.cases[0]["params"]["sim_time"] - 2 * length / 0.2) < 1e-3


def test_flow_throughs_are_refused_without_an_inlet() -> None:
    try:
        campaign_module.plan("ft", "laminarFlamelet", "T_inlet", [300.0, 400.0],
                             {}, 1, 1, flow_throughs=2)
    except ValueError as exc:
        assert "sim_time" in str(exc)
    else:
        raise AssertionError("expected a ValueError")


def test_a_prepared_campaign_is_built_but_not_started() -> None:
    import of_writer
    original = tools._run_run_cases
    launched = []
    tools._run_run_cases = lambda args: launched.append(args) or {"status": "launched"}
    try:
        result = tools.execute_tool("run_campaign", {
            "case_type": "laminarFlamelet", "variable": "T_inlet",
            "values": [300.0, 310.0], "sim_timestep": 2e-6, "sim_time": 1e-4,
            "n_procs_per_case": 1, "n_concurrent_cases": 1,
            "campaign_name": "pytestPrep", "launch": False,
        })
        assert result["status"] == "prepared" and not launched
        assert campaign_module.load("pytestPrep").status == "prepared"
        started = tools.execute_tool("start_campaign", {"campaign_name": "pytestPrep"})
        assert started["status"] == "launched" and len(launched) == 1
        again = tools.execute_tool("start_campaign", {"campaign_name": "pytestPrep"})
        assert again["status"] == "error", "a campaign must not start twice"
    finally:
        tools._run_run_cases = original
        for case in ("laminarFlamelet_pytestPrep_T_inlet_300",
                     "laminarFlamelet_pytestPrep_T_inlet_310"):
            of_writer.delete_case(case)
        import shutil
        shutil.rmtree(campaign_module.CAMPAIGN_DIR / "pytestPrep", ignore_errors=True)


def test_starting_a_campaign_is_a_run_step() -> None:
    import chat
    assert "start_campaign" in chat._RUN_TOOLS


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
