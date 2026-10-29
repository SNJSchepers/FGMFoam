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

"""Tests for the pre-run sanity gate.

These checks are what stands between a mistyped condition and a parallel job
that burns hours to produce nothing. They are pure Python and involve no model,
so they are tested directly against a real built case. Runs offline.

    agent/.venv/bin/python agent/tests/test_precheck.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import of_writer  # noqa: E402
import precheck  # noqa: E402
import tools  # noqa: E402

CASE = "laminarFlamelet_pytest_check"


def _fresh_case(**overrides) -> Path:
    args = {
        "case_type": "laminarFlamelet", "T_inlet": 300.0,
        "case_name": "pytest_check",
    }
    args.update(overrides)
    result = tools.execute_tool("build_case", args)
    assert result["status"] == "built", result
    return Path(result["case_dir"])


def _codes(case_dir: Path, severity: str | None = None) -> set[str]:
    return {
        f.code for f in precheck.check_case(case_dir)
        if severity is None or f.severity == severity
    }


def _set_delta_t(case_dir: Path, value: float) -> None:
    """Rewrite deltaT whatever it currently is — it is estimated per case now."""
    import re

    path = case_dir / "system" / "controlDict"
    text = path.read_text()
    patched, n = re.subn(r"deltaT\s+\S+;", f"deltaT            {value:g};", text)
    assert n == 1, "expected exactly one deltaT entry"
    path.write_text(patched)


def _patch(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    assert old in text, f"{path.name}: {old!r} not found"
    path.write_text(text.replace(old, new))


def test_a_freshly_built_case_has_no_errors() -> None:
    case_dir = _fresh_case()
    try:
        assert _codes(case_dir, "error") == set()
        assert precheck.summarise(precheck.check_case(case_dir))["n_errors"] == 0
    finally:
        of_writer.delete_case(CASE)


def test_pressure_must_match_the_table() -> None:
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "0" / "p", "uniform 101325.0", "uniform 500000.0")
        assert "pressure_mismatch" in _codes(case_dir, "error")
    finally:
        of_writer.delete_case(CASE)


def test_enthalpy_outside_the_table_is_an_error() -> None:
    """Every thermochemical quantity is looked up; outside the table there is
    nothing to look up."""
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "0" / "ht", "uniform -180661.27016190128", "uniform -9.9e6")
        assert "enthalpy_outside_table" in _codes(case_dir, "error")
    finally:
        of_writer.delete_case(CASE)


def test_progress_variable_outside_the_table_is_an_error() -> None:
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "system" / "setFieldsDict",
               "volScalarFieldValue Yc  0.002445377711595195",
               "volScalarFieldValue Yc  0.5")
        assert "progress_variable_outside_table" in _codes(case_dir, "error")
    finally:
        of_writer.delete_case(CASE)


def test_write_interval_below_the_timestep_is_a_warning() -> None:
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "system" / "controlDict",
               "writeInterval     0.001", "writeInterval     1e-09")
        assert "write_interval_below_timestep" in _codes(case_dir, "warning")
    finally:
        of_writer.delete_case(CASE)


def test_mesh_facts_are_reported() -> None:
    case_dir = _fresh_case()
    try:
        findings = {f.code: f for f in precheck.check_case(case_dir)}
        mesh = findings["mesh"]
        assert mesh.details["n_cells"] > 0
        # The flamelet mesh is strongly graded towards the flame; the smallest
        # spacing must reflect that, not an average.
        assert 0 < mesh.details["min_spacing"] < 1e-3
    finally:
        of_writer.delete_case(CASE)


def test_a_quiescent_case_says_courant_is_not_determinable() -> None:
    """Reporting Co = 0 for a flame that has not started would be misleading."""
    case_dir = _fresh_case()
    try:
        codes = _codes(case_dir)
        assert "courant_not_determinable" in codes
        assert "courant_ok" not in codes and "courant_high" not in codes
    finally:
        of_writer.delete_case(CASE)


def test_a_moving_case_gets_a_courant_estimate() -> None:
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "0" / "U", "uniform (0.0 0.0 0.0)", "uniform (5.0 0.0 0.0)")
        findings = {f.code: f for f in precheck.check_case(case_dir)}
        courant = findings.get("courant_high") or findings.get("courant_ok")
        assert courant is not None
        details = courant.details
        expected = details["velocity"] * details["delta_t"] / details["min_spacing"]
        assert abs(details["courant"] - expected) < 1e-9
    finally:
        of_writer.delete_case(CASE)


def test_a_too_large_step_is_flagged_with_a_suggestion() -> None:
    case_dir = _fresh_case()
    try:
        _patch(case_dir / "0" / "U", "uniform (0.0 0.0 0.0)", "uniform (5.0 0.0 0.0)")
        _set_delta_t(case_dir, 1e-3)
        findings = {f.code: f for f in precheck.check_case(case_dir)}
        assert "courant_high" in findings
        suggested = findings["courant_high"].details["suggested_timestep"]
        assert suggested < 1e-03
    finally:
        of_writer.delete_case(CASE)


def test_summarise_grades_by_worst_severity() -> None:
    case_dir = _fresh_case()
    try:
        assert precheck.summarise(precheck.check_case(case_dir))["status"] in (
            "ok", "warnings",
        )
        _patch(case_dir / "0" / "p", "uniform 101325.0", "uniform 500000.0")
        assert precheck.summarise(precheck.check_case(case_dir))["status"] == "errors"
    finally:
        of_writer.delete_case(CASE)


def test_check_case_is_exposed_as_a_read_only_tool() -> None:
    case_dir = _fresh_case()
    try:
        result = tools.execute_tool("check_case", {"case_name": CASE})
        assert result["status"] in ("ok", "warnings")
        assert result["case_name"] == CASE
        assert case_dir.is_dir(), "checking must not modify or remove the case"
    finally:
        of_writer.delete_case(CASE)


def _slot_case(sim_time: float) -> Path:
    result = tools.execute_tool("build_case", {
        "inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0,
        "sim_time": sim_time, "case_name": "pytest_check",
    })
    assert result["status"] == "built", result
    return Path(result["case_dir"])


def test_a_run_shorter_than_one_flow_through_is_flagged() -> None:
    """0.03 s on a 28 mm domain at 0.2 m/s is a fifth of one pass — the fresh
    gas has not crossed the domain, so the result is still the initial state."""
    case_dir = _slot_case(0.03)
    try:
        findings = {f.code: f for f in precheck.check_case(case_dir)}
        warning = findings["shorter_than_one_flow_through"]
        assert warning.severity == "warning"
        details = warning.details
        assert abs(details["flow_through_time"]
                   - details["domain_length"] / details["inlet_velocity"]) < 1e-12
        assert details["flow_throughs"] < 1.0
        assert details["recommended_end_time"] > details["flow_through_time"]
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_check")


def test_a_long_enough_run_passes_the_residence_check() -> None:
    case_dir = _slot_case(0.5)          # ~3.6 flow-throughs at 0.2 m/s
    try:
        codes = _codes(case_dir)
        assert "residence_time_ok" in codes
        assert "shorter_than_one_flow_through" not in codes
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_check")


def test_residence_time_does_not_apply_without_inflow() -> None:
    """The flamelet starts at rest; the flame sets its own time scale."""
    case_dir = _fresh_case()
    try:
        assert "residence_time_not_applicable" in _codes(case_dir)
    finally:
        of_writer.delete_case(CASE)


def test_missing_case_is_reported_not_raised() -> None:
    result = tools.execute_tool("check_case", {"case_name": "no_such_case"})
    assert result["status"] == "error"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
