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

"""Tests for automatic time-step selection.

The step is fixed for the whole run, so getting it wrong cannot be corrected
later. These check the formula against the values that real runs actually held.

    agent/.venv/bin/python agent/tests/test_timestep.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import timestep  # noqa: E402


@dataclass
class _Thermo:
    density_ratio: float = 6.13


def test_slot_burner_uses_the_burnt_gas_velocity() -> None:
    """Not the inlet velocity: the gas leaves the flame ~6x faster."""
    params = {"flame_thickness": 6.2278e-04, "inlet_velocity": 0.2}
    result = timestep.estimate("laminarSlotBurner", params, _Thermo())
    assert result is not None
    assert abs(result.velocity - 0.2 * 6.13) < 1e-6
    assert abs(result.cell_size - 6.2278e-05) < 1e-12
    # A real run at this step held a settled Courant of 0.11 against a 0.1 target.
    assert abs(result.timestep - 5.0e-06) < 1e-9


def test_the_estimate_matches_the_courant_actually_observed() -> None:
    """dt x U / dx must come back to the target, which is the whole point."""
    params = {"flame_thickness": 6.2278e-04, "inlet_velocity": 0.2}
    result = timestep.estimate("laminarSlotBurner", params, _Thermo())
    implied = result.timestep * result.velocity / result.cell_size
    assert abs(implied - timestep.TARGET_CFL) < 0.02


def test_a_faster_inlet_gives_a_smaller_step() -> None:
    slow = timestep.estimate(
        "laminarSlotBurner",
        {"flame_thickness": 6.2278e-04, "inlet_velocity": 0.2}, _Thermo())
    fast = timestep.estimate(
        "laminarSlotBurner",
        {"flame_thickness": 6.2278e-04, "inlet_velocity": 2.0}, _Thermo())
    assert fast.timestep < slow.timestep


def test_a_finer_mesh_gives_a_smaller_step() -> None:
    coarse = timestep.estimate(
        "laminarSlotBurner",
        {"flame_thickness": 6.2278e-04, "inlet_velocity": 0.2}, _Thermo())
    fine = timestep.estimate(
        "laminarSlotBurner",
        {"flame_thickness": 6.2278e-05, "inlet_velocity": 0.2}, _Thermo())
    assert fine.timestep < coarse.timestep


def test_the_step_is_snapped_to_divide_the_write_interval() -> None:
    """Otherwise the solver lands on 0.019998 instead of 0.02."""
    params = {
        "flame_thickness": 6.2278e-04, "inlet_velocity": 0.2,
        "sim_time": 0.05, "sim_write_steps": 50,
    }
    result = timestep.estimate("laminarSlotBurner", params, _Thermo())
    write_interval = 0.05 / 50
    ratio = write_interval / result.timestep
    assert abs(ratio - round(ratio)) < 1e-9


def test_no_estimate_without_a_velocity_scale() -> None:
    """A missing input must fall back, never guess."""
    assert timestep.estimate(
        "laminarSlotBurner", {"flame_thickness": 6.2278e-04}, _Thermo()) is None
    assert timestep.estimate(
        "laminarSlotBurner", {"inlet_velocity": 0.2}, _Thermo()) is None
    # A density ratio of 1 means the equilibrium calculation gave nothing usable.
    assert timestep.estimate(
        "laminarSlotBurner",
        {"flame_thickness": 6.2278e-04, "inlet_velocity": 0.2},
        _Thermo(density_ratio=1.0)) is None
    assert timestep.estimate("unknownCase", {}, _Thermo()) is None


def test_flamelet_cell_size_matches_its_mesh() -> None:
    """Its blockMeshDict refines to a fixed 10 um; the estimate must agree."""
    assert timestep.cell_size("laminarFlamelet", {}) == 1.0e-5
    assert timestep.cell_size("laminarSlotBurner", {"flame_thickness": 1e-3}) == 1e-4


def test_round_to_one_significant_figure() -> None:
    assert timestep._round_to_1_sigfig(2.4e-6) == 2e-6
    assert timestep._round_to_1_sigfig(3.7e-7) == 4e-7
    assert timestep._round_to_1_sigfig(0.0) == 0.0


def test_flamelet_estimate_uses_the_flame_speed() -> None:
    """No inlet drives a freely propagating flame, so S_L sets the scale.

    Uses the on-disk free-flame cache when present so the test stays offline;
    skipped rather than solving a 1D flame if the cache is cold.
    """
    import cantera_tool

    flame = None
    try:
        cached = __import__("json").loads(
            cantera_tool._FLAME_CACHE.read_text(encoding="utf-8"))
        key = cantera_tool._cache_key(0.7, 300.0, 101325.0)
        if key in cached:
            flame = cached[key]
    except Exception:  # noqa: BLE001
        pass
    if flame is None:
        print("    (skipped: free-flame cache is cold)")
        return

    result = timestep.estimate(
        "laminarFlamelet", {"T_inlet": 300.0}, _Thermo())
    assert result is not None
    expected = flame["flame_speed"] * flame["density_ratio"]
    assert abs(result.velocity - expected) < 1e-6
    assert "flame speed" in result.velocity_source


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
