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

"""Tests for post-run interpretation.

The measurements are checked against synthetic profiles with a *known* answer —
a front moving at a prescribed speed, a flame tip at a prescribed position — so
a wrong answer is detectable without a multi-hour solver run. The file-format
parsing is checked against a real `postProcess` output file where one is
available.

    agent/.venv/bin/python agent/tests/test_postrun.py
"""

from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import postrun  # noqa: E402
from postrun import Sample  # noqa: E402

T_UNBURNT = 300.0
T_BURNT = 1840.0          # adiabatic flame temperature at phi = 0.7
THICKNESS = 5.0e-4        # flame thickness used for the synthetic profiles


def _profile(x: list[float], front: float, rising: bool = True) -> list[float]:
    """A smooth temperature front centred on ``front``."""
    out = []
    for xi in x:
        s = (xi - front) / THICKNESS
        ramp = 0.5 * (1.0 + math.tanh(s if rising else -s))
        out.append(T_UNBURNT + (T_BURNT - T_UNBURNT) * ramp)
    return out


def _flamelet_samples(
    speed: float, u_unburnt: float, n_times: int = 40,
    transient: int = 10, dt: float = 2.0e-4,
) -> list[Sample]:
    """A front propagating towards -x at a known speed relative to the gas.

    The first ``transient`` samples move at a deliberately wrong speed, standing
    in for the step-function initial condition relaxing into a real flame.
    """
    x = [-0.05 + i * 1.0e-5 for i in range(10_000)]
    samples = []
    position = 0.015
    for i in range(n_times):
        time = i * dt
        if i > 0:
            # The front moves towards -x. In the lab frame it travels at the
            # flame speed *plus* the speed of the gas it is moving into, which
            # is what the analysis has to undo.
            rate = speed if i > transient else speed * 3.0
            position -= (rate + abs(u_unburnt)) * dt
        temperature = _profile(x, position, rising=True)
        velocity = [u_unburnt] * len(x)
        samples.append(Sample(time, list(x), {"T": temperature, "U": velocity}))
    return samples


def test_crossings_are_interpolated_between_points() -> None:
    x = [0.0, 1.0, 2.0]
    y = [0.0, 10.0, 0.0]
    assert postrun.crossings(x, y, 5.0) == [0.5, 1.5]
    assert postrun.crossings(x, y, 50.0) == []


def test_flame_speed_recovers_a_known_value() -> None:
    """The measured speed must match the one used to build the data."""
    speed, u_unburnt = 0.20, -0.8
    result = postrun.analyse_flamelet(_flamelet_samples(speed, u_unburnt))
    assert result["status"] == "measured", result
    # Front speed in the lab frame is the sum; removing the gas velocity must
    # recover the flame speed itself.
    assert abs(result["laminar_flame_speed"] - speed) < 0.02 * speed
    assert result["fit_r2"] > 0.999


def test_the_initialisation_transient_is_excluded() -> None:
    """Including the step-function start would inflate the speed threefold."""
    speed = 0.20
    samples = _flamelet_samples(speed, 0.0, transient=12)
    result = postrun.analyse_flamelet(samples)
    assert result["window"]["from_time"] > 0.0
    assert "transient" in result["transient_note"].lower()
    assert abs(result["laminar_flame_speed"] - speed) < 0.05 * speed


def test_too_few_snapshots_is_reported_not_guessed() -> None:
    result = postrun.analyse_flamelet(_flamelet_samples(0.2, 0.0, n_times=3))
    assert result["status"] == "insufficient_data"
    assert "sim_write_steps" in result["message"]


def _slot_samples(tip: float, upstream_hot: bool = False,
                  hottest: float = T_BURNT, wake_hot: bool = True) -> list[Sample]:
    """Both symmetry lines of a half-slot with a flame tip at a known position.

    The wake line (y = 0) has a plate-shaped gap and fills with burnt gas just
    behind it; the slot centreline (y = h) has no gap and crosses at ``tip``.
    """
    plate = (0.0062278, 0.0093417)     # l_in .. l_in + l_plate for lF = 6.2278e-4
    dx = 6.2278e-5
    x = [i * dx for i in range(int(0.028 / dx))]

    def profile(xs: list[float], front: float) -> list[float]:
        temperature = []
        for xi in xs:
            if xi < plate[0]:
                temperature.append(hottest if upstream_hot else T_UNBURNT)
            else:
                ramp = 0.5 * (1.0 + math.tanh((xi - front) / THICKNESS))
                temperature.append(T_UNBURNT + (hottest - T_UNBURNT) * ramp)
        return temperature

    wake_x = [xi for xi in x if not (plate[0] < xi < plate[1])]
    wake_front = plate[1] + 5e-4 if wake_hot else 1.0
    centre = Sample(0.1, x, {"T": profile(x, tip)})
    return [Sample(0.1, wake_x, {"T": profile(wake_x, wake_front)}, centre)]


def test_flame_length_is_measured_from_the_holder() -> None:
    tip = 0.015
    result = postrun.analyse_slot_burner(_slot_samples(tip))
    assert result["status"] == "measured", result
    # The gap edges must be found from the sampled points alone.
    assert abs(result["plate_trailing_edge"] - 0.0093417) < 1e-4
    assert abs(result["flame_tip_x"] - tip) < 2e-4
    assert abs(
        result["flame_length"] - (result["flame_tip_x"] - result["plate_trailing_edge"])
    ) < 1e-12


def test_the_tip_comes_from_the_slot_centreline_not_the_wake() -> None:
    """The wake is burnt just behind the plate at any velocity; reading the tip
    there gave the same ~0.5 mm for every inlet velocity of a campaign."""
    short = postrun.analyse_slot_burner(_slot_samples(0.012))
    long = postrun.analyse_slot_burner(_slot_samples(0.020))
    assert abs((long["flame_length"] - short["flame_length"]) - 0.008) < 2e-4
    assert long["flame_length"] > 0.01


def test_flashback_is_reported_instead_of_a_length() -> None:
    result = postrun.analyse_slot_burner(_slot_samples(0.015, upstream_hot=True))
    assert result["status"] == "flashback"
    assert "flame length" in result["message"]
    assert "flame_length" not in result


def test_blow_off_when_nothing_downstream_ignites() -> None:
    result = postrun.analyse_slot_burner(_slot_samples(0.015, hottest=600.0))
    assert result["status"] == "blow_off"
    assert result["max_axis_temperature"] < postrun.FLAME_TEMPERATURE


def test_a_tip_beyond_the_outlet_is_not_blow_off() -> None:
    """Burnt wake, cold centreline to the outlet: anchored, just too long."""
    result = postrun.analyse_slot_burner(_slot_samples(0.05))
    assert result["status"] == "longer_than_domain", result
    assert "flame_length" not in result
    assert "outlet_length" in result["message"]
    assert result["flame_length_lower_bound"] > 0.018


def test_a_cold_wake_and_centreline_is_blow_off() -> None:
    result = postrun.analyse_slot_burner(_slot_samples(0.05, wake_hot=False))
    assert result["status"] == "blow_off"


def test_a_tip_against_the_outlet_is_flagged() -> None:
    """A length set by the domain, not the flame, must not read as a result."""
    result = postrun.analyse_slot_burner(_slot_samples(0.0275))
    assert result["status"] == "longer_than_domain"
    assert "flame_length" not in result
    assert "bounded by the domain" in result["message"]


def test_plate_gap_detection_ignores_mere_grading() -> None:
    graded = [0.0]
    step = 1e-5
    for _ in range(200):
        graded.append(graded[-1] + step)
        step *= 1.02          # strong grading, but no obstacle
    assert postrun.axis_gap(graded) is None


def test_log_parsing_reads_courant_and_health() -> None:
    log = """
Time = 0.001
Courant Number mean: 0.01 max: 0.35
time step continuity errors : sum local = 1e-09, global = 2e-12, cumulative = 3e-11
ExecutionTime = 10.5 s  ClockTime = 11 s

Time = 0.002
Courant Number mean: 0.02 max: 1.42
time step continuity errors : sum local = 1e-09, global = 2e-12, cumulative = 4e-11
ExecutionTime = 21.0 s  ClockTime = 22 s
"""
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        (case / "log.FGMFoam").write_text(log)
        result = postrun.parse_log(case)
    assert result["n_timesteps"] == 2
    assert result["time_reached"] == 0.002
    assert result["courant_max"] == 1.42
    assert result["execution_time_s"] == 21.0
    assert result["cumulative_continuity_error"] == 4e-11
    assert not result["fatal_error"]
    # Above 1 the fixed step could not adapt, so it must be called out.
    assert "fixed time step" in result["courant_note"]


def test_the_ignition_spike_is_separated_from_settled_behaviour() -> None:
    """A run that sat at 0.11 for 6000 steps once reported "max Courant 4.11".

    The step-function ignition spikes the Courant number for a few steps. Judging
    the run on that number wrongly condemns a perfectly resolved simulation, so
    the settled value is reported alongside it.
    """
    lines = ["Courant Number mean: 0.01 max: 4.11"]
    lines += ["Courant Number mean: 0.01 max: 0.11"] * 200
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        (case / "log.FGMFoam").write_text("\n".join(lines))
        result = postrun.parse_log(case)
    assert result["courant_max"] == 4.11
    assert result["courant_max_settled"] == 0.11
    assert result["startup_steps_excluded"] == postrun.STARTUP_STEPS
    # The note must reassure rather than condemn.
    assert "start-up spike is" in result["courant_note"]
    assert "does not by itself" in result["courant_note"]


def test_a_genuinely_high_courant_is_still_flagged() -> None:
    """Settled above 1 is a real problem and must not be explained away."""
    lines = ["Courant Number mean: 0.5 max: 2.4"] * 200
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        (case / "log.FGMFoam").write_text("\n".join(lines))
        result = postrun.parse_log(case)
    assert result["courant_max_settled"] == 2.4
    assert "held 2.40" in result["courant_note"]
    assert "smaller sim_timestep" in result["courant_note"]


def test_log_parsing_flags_a_fatal_run() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        (case / "log.FGMFoam").write_text("--> FOAM FATAL ERROR:\nsomething broke\n")
        result = postrun.parse_log(case)
    assert result["fatal_error"] is True
    assert result["n_timesteps"] == 0


def test_written_times_ignores_zero_and_orig() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        for name in ("0", "0.orig", "0.01", "0.02", "constant", "system"):
            (case / name).mkdir()
        assert postrun.written_times(case) == [0.01, 0.02]


def test_incomplete_sampling_is_refused_not_measured() -> None:
    """A partially sampled case looks exactly like a flame that went out.

    This actually happened: a campaign analysed a case moments after it
    finished, postProcess saw only time 0, and the slot burner was reported as
    a confident "blow_off". Incomplete coverage must be an error, not a result.
    """
    samples = [Sample(0.0, [0.0, 1.0], {"T": [300.0, 300.0]})]

    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        for name in ("0", "0.01", "0.02", "0.03"):
            (case / name).mkdir()
        latest = postrun.written_times(case)[-1]
        # The guard compares the last sampled time against what is on disk.
        assert samples[-1].time < latest
        # And the analysis of such a sample set must not claim a measurement.
        assert postrun.analyse_flamelet(samples)["status"] == "insufficient_data"


def test_missing_log_is_reported() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        assert postrun.parse_log(Path(tmp))["status"] == "no_log"


def test_sampled_file_columns_match_openfoam_output() -> None:
    """Parse a real `postProcess` raw file when one is present in the repo."""
    from config import REPO_ROOT

    candidates = list(
        (REPO_ROOT / "testcases").glob("*/postProcessing/*/*/axis_*.xy")
    )
    if not candidates:
        print("    (skipped: no sampled output in the tree)")
        return
    xs, values = postrun._parse_xy(candidates[0], ["T", "Yc", "U"])
    assert len(xs) > 100
    assert set(values) == {"T", "Yc", "U"}
    assert all(len(v) == len(xs) for v in values.values())
    # Sanity: the unburnt state of these cases.
    assert min(values["T"]) >= 250.0
    assert min(values["Yc"]) >= 0.0


def test_a_front_still_accelerating_is_not_converged() -> None:
    # A 300 K flamelet stopped at 6 ms: 0.166 -> 0.185 m/s across its window
    # passed the 5 % scatter test and was reported 3 % low.
    times = [0.0005 * i for i in range(1, 13)]
    steady = [0.04 - 0.185 * t for t in times]
    start, drift = postrun.settled_start(times, steady, 2)
    assert drift is None
    rising = [0.04 - (0.15 * t + 3.0 * t * t) for t in times]   # speed 0.153 -> 0.186
    start, drift = postrun.settled_start(times, rising, 2)
    assert drift is not None and drift > postrun.TREND_TOLERANCE


def test_snapshots_draw_each_case_from_its_cell_fields() -> None:
    import tempfile
    from types import SimpleNamespace as NS
    import of_writer
    import snapshots

    def field(path, values):
        body = "\n".join(f"{v:g}" for v in values)
        path.write_text(f"FoamFile {{}}\ninternalField   nonuniform List<scalar> \n"
                        f"{len(values)}\n(\n{body}\n)\n;\n")
    with tempfile.TemporaryDirectory() as tmp:
        saved = of_writer.TESTCASES_DIR
        of_writer.TESTCASES_DIR = Path(tmp)
        try:
            rows = []
            for k, u in enumerate((0.2, 0.4)):
                time = Path(tmp) / f"laminarSlotBurner_snap_{k}" / "0.1"
                time.mkdir(parents=True)
                xs = [i * 1e-4 for i in range(30) for _ in range(10)]
                ys = [j * 1e-4 for _ in range(30) for j in range(10)]
                field(time / "Cx", xs); field(time / "Cy", ys)
                field(time / "T", [300 + 1500 * (x > 1e-3 * (1 + 5 * u)) for x in xs])
                rows.append({"case_name": time.parent.name, "inlet_velocity": u,
                             "flame_length": 0.005 + 0.02 * u})
            plan = NS(name="snap", case_type="laminarSlotBurner", variable="inlet_velocity",
                      directory=Path(tmp))
            spec = {"field": "T", "unit": "K", "range": (280.0, 1800.0), "mirror": "ymax",
                    "upright": True}
            path, error = snapshots.render(plan, {"rows": rows}, spec,
                                           {"key": "flame_length", "label": "Flame length [m]",
                                            "title": "Flame length"}, "Inlet velocity [m/s]", "m/s")
            assert error is None and path.is_file() and path.stat().st_size > 10000
        finally:
            of_writer.TESTCASES_DIR = saved


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
