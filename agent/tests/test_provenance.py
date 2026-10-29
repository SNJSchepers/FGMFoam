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

"""Tests for the per-case provenance record and the FGM table reader.

The record is what makes a generated case self-describing — which numbers were
the user's choice, which the agent filled in, and exactly which table and
mechanism produced the thermochemistry. Runs offline.

    agent/.venv/bin/python agent/tests/test_provenance.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import fgm_table  # noqa: E402
import of_writer  # noqa: E402
import provenance  # noqa: E402
import tools  # noqa: E402
from config import REPO_ROOT  # noqa: E402

TABLE = REPO_ROOT / "tables" / "database.fgm"


def _build(**overrides):
    args = {
        "case_type": "laminarFlamelet", "T_inlet": 300.0,
        "case_name": "pytest_prov",
    }
    args.update(overrides)
    return tools.execute_tool("build_case", args)


def test_table_header_is_parsed() -> None:
    info = fgm_table.read_table(TABLE)
    assert info.pressure == 101325.0
    assert info.dimension == 2
    assert info.datasize == (250, 100)
    assert info.n_rows == 250 * 100, "[DATASIZE] and [DATA] must agree"
    assert info.variables[:2] == ("CV1", "CV2")
    assert "TEMPERATURE" in info.variables


def test_table_control_variable_ranges_are_sane() -> None:
    info = fgm_table.read_table(TABLE)
    yc_lo, yc_hi = info.range_of(0)
    ht_lo, ht_hi = info.range_of(1)
    assert yc_lo <= 0.0 < yc_hi, "progress variable starts at the unburnt state"
    assert ht_lo < 0.0 < ht_hi, "enthalpy spans both signs"
    assert yc_hi < 1.0, "Yc is a specific mole number (kmol/kg), not normalised"


def test_table_reader_is_cached_per_file_identity() -> None:
    first = fgm_table.read_table(TABLE)
    second = fgm_table.read_table(TABLE)
    assert first is second


def test_every_built_case_gets_a_record() -> None:
    result = _build()
    try:
        assert result["status"] == "built"
        assert "provenance_error" not in result, result.get("provenance_error")
        record = provenance.read(Path(result["case_dir"]))
        assert record is not None
        assert record["schema"] == provenance.SCHEMA
        assert record["case"]["type"] == "laminarFlamelet"
        assert record["thermochemistry"]["mechanism"] == "gri30.yaml"
        assert record["fgm_table"]["sha256"] == fgm_table.read_table(TABLE).sha256
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_prov")


def test_record_separates_the_users_choices_from_the_agents() -> None:
    """The whole point: a reviewer must see which numbers were deliberate."""
    result = _build(sim_time=0.05)
    try:
        record = provenance.read(Path(result["case_dir"]))
        sources = record["parameter_sources"]
        assert sources["T_inlet"] == "user"
        assert sources["sim_time"] == "user"
        # Computed for these conditions rather than taken from the shipped
        # value — see test_an_estimated_step_is_recorded_as_estimated.
        assert sources["sim_timestep"] == "estimated"
        # Generic fallback, not specific to this case.
        assert sources["n_procs"] == "default"
        assert set(sources) == set(record["parameters"])
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_prov")


def test_record_carries_no_absolute_paths() -> None:
    """Records get shared and published; they must not leak a home directory."""
    result = _build()
    try:
        text = json.dumps(provenance.read(Path(result["case_dir"])))
        assert str(REPO_ROOT) not in text
        assert "/home/" not in text
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_prov")


def test_record_captures_the_routing_model() -> None:
    tools.set_session_model("anthropic", "claude-sonnet-5")
    try:
        result = _build()
        record = provenance.read(Path(result["case_dir"]))
        assert record["model"] == {
            "provider": "anthropic", "chat_model": "claude-sonnet-5",
        }
        assert isinstance(record["warnings"], list)
    finally:
        tools.set_session_model(None, None)
        of_writer.delete_case("laminarFlamelet_pytest_prov")


def test_warnings_are_carried_into_the_record() -> None:
    """Checked on the record itself, so it does not depend on which warnings a
    particular build happens to raise."""
    import case_templates
    from config import AGENT_VERSION

    template = case_templates.TEMPLATE_BY_NAME["laminarFlamelet"]
    warning = {"parameter": "sim_timestep", "suggested": 5e-7}
    record = provenance.build_record(
        case_dir=Path("testcases/laminarFlamelet_x"), template=template,
        params={"sim_time": 0.1, "sim_write_steps": 100, "sim_timestep": 2e-6},
        supplied={"sim_time"},
        thermo=type("T", (), {"Yc_u": 0.0, "Yc_b": 1.0, "ht": -1.0})(),
        table_path=TABLE, repo_root=REPO_ROOT, version=AGENT_VERSION,
        warnings=[warning],
    )
    assert record["warnings"] == [warning]
    assert record["parameter_sources"]["sim_time"] == "user"


def test_an_estimated_step_is_recorded_as_estimated() -> None:
    result = tools.execute_tool("build_case", {
        "inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0,
        "case_name": "pytest_prov",
    })
    try:
        record = provenance.read(Path(result["case_dir"]))
        assert record["parameter_sources"]["sim_timestep"] == "estimated"
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_prov")


def test_a_broken_record_never_fails_a_build() -> None:
    """A provenance failure must degrade, not lose the case."""
    original = provenance.build_record

    def boom(**_kwargs):
        raise RuntimeError("synthetic failure")

    provenance.build_record = boom
    try:
        result = _build()
        assert result["status"] == "built"
        assert "synthetic failure" in result["provenance_error"]
        assert Path(result["case_dir"]).is_dir()
    finally:
        provenance.build_record = original
        of_writer.delete_case("laminarFlamelet_pytest_prov")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
