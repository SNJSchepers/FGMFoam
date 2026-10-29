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

"""Interpret a finished run: did it behave, and what did the flame do?

Three things are reported, all measured rather than inferred by a model:

* **How the run behaved** — from `log.FGMFoam`: the Courant number actually
  reached (the pre-run gate only estimates it), the simulated time reached,
  execution time, and whether the solver diverged or was still healthy.
* **The laminar flame speed**, for the 1D flamelet: the front is tracked over
  the written times and its speed taken relative to the unburnt gas. The case is
  ignited as a *step function* in the progress variable, which is not a real
  flame structure, so the initial transient is detected and excluded rather than
  averaged in.
* **The flame length**, for the slot burner: the axial position on the slot
  centreline where the temperature first crosses a threshold, measured from the
  flame holder's trailing edge — together with the ways that measurement can be
  meaningless, **flashback** (flame upstream of the holder), **blow-off** (no
  flame left in the domain) and a tip beyond the outlet (**longer than the
  domain**).

Field data comes from OpenFOAM's own `postProcess -func` line sampling, not from
a reimplementation of its interpolation. Two lines are sampled, along the lower
and upper edges of the mesh bounding box; the flame holder is located from the
*gap* in the lower one, where the solid plate blocks it — so no `blockMeshDict`
is parsed.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import of_mesh
import provenance

# Temperature marking the flame front. Methane's autoignition temperature is
# about 810 K (537 C); 900 K sits just above it, comfortably inside the reaction
# zone and far from both the unburnt (300 K) and burnt (~1800 K at phi = 0.7)
# states, which makes the crossing sharp and insensitive to the exact value.
# What matters for comparing cases is using the *same* threshold throughout.
FLAME_TEMPERATURE: float = 900.0

# Name of the sampling function object written into the case.
SAMPLE_FUNC: str = "agentAxisSample"

# Fields written as three columns by the `raw` set format.
VECTOR_FIELDS: frozenset[str] = frozenset({"U"})

# Upper bound on sampling points, so a very fine mesh cannot make this slow.
MAX_SAMPLE_POINTS: int = 20000

# Fraction of the written times always discarded as initialisation, before the
# steady window is searched for.
MIN_TRANSIENT_FRACTION: float = 0.25

# Relative scatter below which the front speed counts as steady.
STEADY_TOLERANCE: float = 0.05

# Time steps treated as start-up when judging the Courant number. The ignition
# discontinuity relaxes within a handful of steps; by step 20 both shipped cases
# are an order of magnitude below their peak.
STARTUP_STEPS: int = 20


@dataclass(frozen=True)
class Sample:
    """One sampled line at one time: positions plus the requested fields.

    ``x``/``fields`` are the line along the domain's lower edge (y = y_min);
    ``upper`` is the same time sampled along its upper edge (y = y_max). On the
    slot burner these are the two symmetry planes: the wake behind the plate and
    the slot centreline, where the flame tip is.
    """

    time: float
    x: list[float]
    fields: dict[str, list[float]]
    upper: Optional["Sample"] = None


# --------------------------------------------------------------------------- #
# The solver log
# --------------------------------------------------------------------------- #


def parse_log(case_dir: Path, name: str = "FGMFoam") -> dict[str, Any]:
    """Extract how the run behaved from ``log.<name>``."""
    path = Path(case_dir) / f"log.{name}"
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return {"status": "no_log", "log": path.name}

    times = [float(m) for m in re.findall(r"^Time = ([0-9.eE+-]+)", text, re.M)]
    courant = [
        float(m) for m in re.findall(r"Courant Number mean: \S+ max: ([0-9.eE+-]+)", text)
    ]
    execution = re.findall(r"ExecutionTime = ([0-9.eE+-]+) s", text)
    # "time step continuity errors : sum local = ..., global = ..., cumulative = ..."
    # The label contains spaces, so match the value directly rather than the shape.
    continuity = re.findall(r"cumulative = ([-+0-9.eE]+)", text)
    # The ignition transient briefly spikes the Courant number: a step-function
    # initial condition is not a flame, and the first few steps relax it. The
    # overall maximum therefore says little about how the run behaved — a run
    # that sat at 0.11 for 6000 steps reported 4.1 because of step 3. Report
    # both, and judge the run on the settled value.
    settled = courant[STARTUP_STEPS:] if len(courant) > STARTUP_STEPS else courant
    fatal = "FOAM FATAL" in text
    # A diverging compressible run shows up as a runaway Courant number or as
    # non-finite values long before the solver gives up.
    nan = bool(re.search(r"\bnan\b|\binf\b", text, re.I))
    max_courant = max(courant) if courant else None

    result: dict[str, Any] = {
        "status": "read",
        "log": path.name,
        "n_timesteps": len(times),
        "time_reached": times[-1] if times else None,
        "execution_time_s": float(execution[-1]) if execution else None,
        "courant_max": max_courant,
        "courant_max_settled": max(settled) if settled else None,
        "courant_final": courant[-1] if courant else None,
        "startup_steps_excluded": STARTUP_STEPS if len(courant) > STARTUP_STEPS else 0,
        "cumulative_continuity_error": (
            float(continuity[-1]) if continuity else None
        ),
        "fatal_error": fatal,
        "non_finite_values": nan,
    }
    settled_max = max(settled) if settled else None
    if settled_max is not None and settled_max > 1.0:
        result["courant_note"] = (
            f"The Courant number held {settled_max:.2f} after the first "
            f"{STARTUP_STEPS} steps. The solver runs at a fixed time step, so it "
            "could not reduce it — treat the result with care and rerun with a "
            "smaller sim_timestep."
        )
    elif max_courant is not None and max_courant > 1.0:
        result["courant_note"] = (
            f"The Courant number peaked at {max_courant:.2f} during the first "
            f"{STARTUP_STEPS} steps, while the step-function initial condition "
            f"relaxed, then settled at {settled_max:.2f}. That start-up spike is "
            "expected and does not by itself make the run unreliable."
        )
    return result


# --------------------------------------------------------------------------- #
# Line sampling through OpenFOAM
# --------------------------------------------------------------------------- #


def _sample_dict_text(
    lines: dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]],
    fields: list[str], n_points: int,
) -> str:
    """Write the `sets` function-object dict used for sampling.

    ``uniform`` rather than ``midPoint``: a `midPoint` set walks through
    connected cells and *stops* at the first obstacle, so on the slot burner it
    ended at the plate's leading edge and never reached the axis downstream of
    it. A uniform set evaluates independent points and simply drops those
    outside the mesh, which both reaches the far side and leaves the solid body
    as a gap — the gap being how the flame holder is located.
    """
    sets = "".join(
        f"""    {name}
    {{
        type        uniform;
        axis        x;
        start       ({start[0]:.10g} {start[1]:.10g} {start[2]:.10g});
        end         ({end[0]:.10g} {end[1]:.10g} {end[2]:.10g});
        nPoints     {n_points};
    }}
"""
        for name, (start, end) in lines.items()
    )
    return f"""// Written by the FGMFoam agent for post-run analysis; safe to delete.
type            sets;
libs            (sampling);
writeControl    onEnd;
setFormat       raw;
interpolationScheme cell;
fields          ({' '.join(fields)});
sets
(
{sets});
"""


def _parse_xy(path: Path, fields: list[str]) -> tuple[list[float], dict[str, list[float]]]:
    """Parse a `raw`-format sampled line.

    Columns are the sampling coordinate followed by each field in order, with
    vectors expanded to three columns — only the first (axial) component of a
    vector is kept, which is the one every measurement here uses.
    """
    xs: list[float] = []
    columns: list[list[float]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [float(v) for v in line.split()]
        xs.append(parts[0])
        columns.append(parts[1:])
    values: dict[str, list[float]] = {}
    offset = 0
    for name in fields:
        # Vectors occupy three columns; only the axial component is ever used.
        span = 3 if name in VECTOR_FIELDS else 1
        values[name] = [row[offset] for row in columns]
        offset += span
    return xs, values


def written_times(case_dir: Path) -> list[float]:
    """Reconstructed time directories present in the case, excluding 0."""
    times: list[float] = []
    for entry in Path(case_dir).iterdir():
        if not entry.is_dir() or entry.name.endswith(".orig"):
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if value > 0.0:
            times.append(value)
    return sorted(times)


def sample_axis(
    case_dir: Path,
    fields: Optional[list[str]] = None,
    timeout: int = 900,
) -> tuple[list[Sample], Optional[str]]:
    """Sample ``fields`` along the domain's lower and upper x-wise edges.

    Both lines are sampled at every written time in one ``postProcess`` call; the
    upper one is attached to each ``Sample`` as ``upper``.

    Returns ``(samples, error)``; ``error`` is a message when OpenFOAM's
    `postProcess` could not be run (it needs an OpenFOAM-sourced shell).
    """
    case_dir = Path(case_dir)
    fields = fields or ["T", "Yc", "U"]
    box = of_mesh.bounds(case_dir)
    if box is None:
        return [], "no mesh found (constant/polyMesh); nothing to sample"

    (xmin, xmax), (ymin, ymax), (zmin, zmax) = box
    # Step just inside the domain so the end points are not on a boundary face.
    inset = (xmax - xmin) * 1e-4
    y_inset = (ymax - ymin) * 1e-4
    z_mid = 0.5 * (zmin + zmax)
    # Resolve at least the finest cell, so a graded mesh is not under-sampled
    # where it matters (both shipped meshes refine towards the flame).
    spacing = of_mesh.min_spacing(case_dir) or (xmax - xmin) / 1000.0
    n_points = int(min(MAX_SAMPLE_POINTS, max(200, round((xmax - xmin) / spacing))))
    lines = {
        name: ((xmin + inset, y, z_mid), (xmax - inset, y, z_mid))
        for name, y in (("axis", ymin + y_inset), ("upper", ymax - y_inset))
    }
    text = _sample_dict_text(lines, fields, n_points)
    dict_path = case_dir / "system" / SAMPLE_FUNC

    # Sampling output left from an earlier call can be *behind* the case: a
    # campaign may analyse a case moments after it finished, before every time
    # directory was reconstructed, and postProcess will not revisit times it
    # already wrote. Stale output once turned an unfinished case into a
    # confident "blow-off", so discard it and sample again. Output written for
    # other lines or fields is just as stale, so a changed dict discards it too.
    out_dir = case_dir / "postProcessing" / SAMPLE_FUNC
    on_disk = written_times(case_dir)
    if out_dir.is_dir():
        sampled = {d.name for d in out_dir.iterdir() if d.is_dir()}
        unchanged = dict_path.is_file() and dict_path.read_text() == text
        if not unchanged or len({f"{t:g}" for t in on_disk} - sampled) > 0:
            shutil.rmtree(out_dir, ignore_errors=True)
    dict_path.write_text(text)

    try:
        proc = subprocess.run(
            ["postProcess", "-func", SAMPLE_FUNC],
            cwd=case_dir, capture_output=True, text=True, timeout=timeout, check=False,
        )
    except FileNotFoundError:
        return [], (
            "postProcess not found — the analysis needs an OpenFOAM-sourced shell"
        )
    except subprocess.SubprocessError as exc:
        return [], f"postProcess failed: {exc}"
    if proc.returncode != 0:
        tail = (proc.stdout or proc.stderr or "").strip().splitlines()[-3:]
        return [], "postProcess failed: " + " | ".join(tail)

    samples: list[Sample] = []
    for time_dir in sorted(
        (d for d in out_dir.iterdir() if d.is_dir()),
        key=lambda d: float(d.name),
    ):
        files = list(time_dir.glob("axis_*.xy"))
        if not files:
            continue
        xs, values = _parse_xy(files[0], fields)
        upper = None
        upper_files = list(time_dir.glob("upper_*.xy"))
        if upper_files:
            ux, uvalues = _parse_xy(upper_files[0], fields)
            if ux:
                upper = Sample(float(time_dir.name), ux, uvalues)
        if xs:
            samples.append(Sample(float(time_dir.name), xs, values, upper))
    if not samples:
        return [], "postProcess produced no sampled data"
    # Refuse to measure a flame from data that does not reach the end of the
    # run: a partially sampled case looks exactly like a flame that went out.
    if on_disk and samples[-1].time < on_disk[-1] * (1.0 - 1e-6):
        return samples, (
            f"sampling only reached t = {samples[-1].time:g} s but the case has "
            f"data up to t = {on_disk[-1]:g} s. The case may still have been "
            "being written; no flame measurement was made from incomplete data"
        )
    return samples, None


# --------------------------------------------------------------------------- #
# Geometry of a profile
# --------------------------------------------------------------------------- #


def crossings(x: list[float], y: list[float], level: float) -> list[float]:
    """Positions where ``y`` crosses ``level``, linearly interpolated."""
    found: list[float] = []
    for i in range(len(x) - 1):
        a, b = y[i] - level, y[i + 1] - level
        if a == 0.0:
            found.append(x[i])
        elif a * b < 0.0:
            found.append(x[i] + (x[i + 1] - x[i]) * a / (a - b))
    return found


def axis_gap(x: list[float]) -> Optional[tuple[float, float]]:
    """Largest gap in the sampled positions — the solid body blocking the axis.

    The slot burner's plate is not meshed, so no sample points exist under it;
    the gap's edges are its leading and trailing faces. Returns None when no gap
    stands out against the local spacing.
    """
    if len(x) < 3:
        return None
    gaps = [(x[i + 1] - x[i], x[i], x[i + 1]) for i in range(len(x) - 1)]
    widest, lo, hi = max(gaps, key=lambda g: g[0])
    typical = sorted(g[0] for g in gaps)[len(gaps) // 2]
    # A real obstacle spans many cells; grading alone never does.
    if typical > 0 and widest > 20 * typical:
        return lo, hi
    return None


def _linear_fit(t: list[float], x: list[float]) -> tuple[float, float, float]:
    """Least-squares fit ``x = slope * t + intercept``; returns (slope, intercept, r2)."""
    n = len(t)
    mean_t = sum(t) / n
    mean_x = sum(x) / n
    s_tt = sum((ti - mean_t) ** 2 for ti in t)
    if s_tt == 0.0:
        return 0.0, mean_x, 0.0
    slope = sum((ti - mean_t) * (xi - mean_x) for ti, xi in zip(t, x)) / s_tt
    intercept = mean_x - slope * mean_t
    ss_tot = sum((xi - mean_x) ** 2 for xi in x)
    ss_res = sum((xi - (slope * ti + intercept)) ** 2 for ti, xi in zip(t, x))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return slope, intercept, r2


def steady_window(
    times: list[float], speeds: list[float],
    min_fraction: float = MIN_TRANSIENT_FRACTION,
    tolerance: float = STEADY_TOLERANCE,
) -> tuple[int, str]:
    """Index from which the front speed has settled.

    The flamelet is ignited as a step in the progress variable — a
    discontinuity, not a flame — so the first instants describe the profile
    relaxing into a real flame structure rather than the flame speed. The fixed
    fraction below is a floor; beyond it the longest trailing window whose
    scatter is within ``tolerance`` is taken, so a slow initialisation is
    excluded too.
    """
    n = len(speeds)
    floor = max(1, int(math.ceil(min_fraction * n)))
    if n - floor < 3:
        return floor, "too few written times to detect a steady window"
    best = None
    for start in range(floor, n - 2):
        window = speeds[start:]
        mean = sum(window) / len(window)
        if mean == 0:
            continue
        spread = math.sqrt(sum((v - mean) ** 2 for v in window) / len(window))
        if spread / abs(mean) <= tolerance:
            best = start
            break
    if best is None:
        return max(floor, n // 2), (
            "the front speed never settled to within "
            f"{tolerance * 100:.0f}%; the trailing half was used and the result "
            "should be treated as provisional"
        )
    return best, (
        f"initialisation transient excluded: the first {times[best]:.4g} s "
        f"({best} of {n} written times)"
    )


def analyse_flamelet(
    samples: list[Sample], threshold: float = FLAME_TEMPERATURE,
    field: str = "T", velocity_field: str = "U",
) -> dict[str, Any]:
    """Laminar flame speed from the tracked front of a 1D freely propagating flame.

    The front position is where the temperature crosses ``threshold``. Its speed
    in the lab frame is the slope of position against time; the flame speed is
    that relative to the **unburnt gas**, which is itself moving because the
    burnt gas expands and pushes it out of the open end. Ignoring that would
    overstate the flame speed by the expansion velocity.
    """
    positions: list[float] = []
    times: list[float] = []
    unburnt_velocity: list[float] = []

    for sample in samples:
        temperature = sample.fields.get(field)
        if temperature is None:
            continue
        found = crossings(sample.x, temperature, threshold)
        if not found:
            continue
        # The unburnt side is whichever end is cold; the front is the crossing
        # closest to it.
        cold_at_start = temperature[0] < temperature[-1]
        position = min(found) if cold_at_start else max(found)
        positions.append(position)
        times.append(sample.time)

        velocity = sample.fields.get(velocity_field)
        if velocity is not None:
            # Average the axial velocity over fresh gas only, away from the
            # front so the expansion through the flame is not included.
            margin = 0.1 * abs(sample.x[-1] - sample.x[0])
            fresh = [
                v for xi, v, ti in zip(sample.x, velocity, temperature)
                if ti < 0.5 * (threshold + min(temperature))
                and (xi < position - margin if cold_at_start else xi > position + margin)
            ]
            if fresh:
                unburnt_velocity.append(sum(fresh) / len(fresh))

    if len(positions) < 4:
        return {
            "status": "insufficient_data",
            "message": (
                "The flame front was found at fewer than four written times, "
                "which is not enough to measure a speed. Run longer, or write "
                "more snapshots (sim_write_steps)."
            ),
            "n_fronts": len(positions),
        }

    instantaneous = [
        (positions[i + 1] - positions[i]) / (times[i + 1] - times[i])
        for i in range(len(positions) - 1)
        if times[i + 1] > times[i]
    ]
    start, note = steady_window(times[1:], [abs(v) for v in instantaneous])
    keep = slice(start, len(positions))
    slope, _, r2 = _linear_fit(times[keep], positions[keep])

    u_unburnt = (
        sum(unburnt_velocity[keep]) / len(unburnt_velocity[keep])
        if len(unburnt_velocity) == len(positions) and unburnt_velocity[keep]
        else 0.0
    )
    flame_speed = abs(slope - u_unburnt)

    result = {
        "status": "measured",
        "laminar_flame_speed": flame_speed,
        "units": "m/s",
        "front_speed_lab_frame": slope,
        "unburnt_gas_velocity": u_unburnt,
        "threshold_K": threshold,
        "fit_r2": r2,
        "window": {
            "from_time": times[start], "to_time": times[-1],
            "n_points": len(positions) - start,
        },
        "transient_note": note,
        "front_positions": [
            {"time": t, "x": x} for t, x in zip(times, positions)
        ],
    }
    if r2 < 0.99:
        result["quality_note"] = (
            f"The front position is only fit to R^2 = {r2:.3f}, so the flame was "
            "not propagating steadily over the window used. Treat the value as "
            "indicative and check the profiles."
        )
    return result


def analyse_slot_burner(
    samples: list[Sample], threshold: float = FLAME_TEMPERATURE,
    field: str = "T",
) -> dict[str, Any]:
    """Flame length on the slot centreline, or flashback / blow-off.

    The domain is one half-slot between two symmetry planes: its lower edge
    (y = 0) runs through the plate and the wake behind it, its upper edge
    (y = h) is the slot centreline. The flame anchors on the plate and closes on
    the slot centreline, so the tip is taken from the upper line, and the plate
    is located from the gap it leaves in the lower one. The tip must not be read
    from the wake: that fills with burnt gas just behind the trailing edge at any
    inlet velocity, which once gave the same ~0.5 mm for every case.

    The length is the distance from the plate's trailing edge to where the slot
    centreline first reaches ``threshold``. Three outcomes make that number
    meaningless and are reported instead: the flame sitting upstream of the plate
    (**flashback**), no flame left downstream of it (**blow-off**), or a tip
    beyond or against the outlet (**longer than the domain** — the flame is
    anchored, but the domain sets the number).
    """
    if not samples:
        return {"status": "insufficient_data", "message": "no sampled data"}
    sample = samples[-1]
    centre = sample.upper
    wake_t = sample.fields.get(field)
    if wake_t is None or centre is None or field not in centre.fields:
        return {
            "status": "insufficient_data",
            "message": (
                f"{field} was not sampled on both the wake line (y = 0) and the "
                "slot centreline (upper edge of the domain)"
            ),
        }
    centre_t = centre.fields[field]

    x_min, x_max = centre.x[0], centre.x[-1]
    gap = axis_gap(sample.x)
    if gap is None:
        return {
            "status": "geometry_not_recognised",
            "message": (
                "No flame holder was found on the sampled axis (the plate should "
                "appear as a gap in the sampled points). The flame length is "
                "defined relative to its trailing edge, so it was not measured."
            ),
        }
    leading_edge, trailing_edge = gap

    hot_upstream = [
        xi
        for line_x, line_t in ((sample.x, wake_t), (centre.x, centre_t))
        for xi, ti in zip(line_x, line_t)
        if xi < leading_edge and ti >= threshold
    ]
    hottest_downstream = max(
        (ti
         for line_x, line_t in ((sample.x, wake_t), (centre.x, centre_t))
         for xi, ti in zip(line_x, line_t) if xi > trailing_edge),
        default=0.0,
    )

    base = {
        "time": sample.time,
        "threshold_K": threshold,
        "plate_leading_edge": leading_edge,
        "plate_trailing_edge": trailing_edge,
        "max_axis_temperature": max(centre_t),
        "max_wake_temperature": max(wake_t),
    }

    if hot_upstream:
        return {
            **base,
            "status": "flashback",
            "message": (
                f"Flashback: the gas reaches {threshold:g} K upstream of the "
                f"flame holder (as far back as x = {min(hot_upstream):.4g} m). "
                "The flame has propagated into the inlet channel, so there is no "
                "flame length to report."
            ),
            "upstream_extent": min(hot_upstream),
        }

    if hottest_downstream < threshold:
        return {
            **base,
            "status": "blow_off",
            "message": (
                f"Blow-off: nowhere downstream of the holder reaches "
                f"{threshold:g} K (hottest {hottest_downstream:.0f} K). The flame "
                "is not anchored — it has been blown out of the domain or "
                "extinguished."
            ),
        }

    longer = {
        **base,
        "status": "longer_than_domain",
        "flame_length_lower_bound": x_max - trailing_edge,
    }
    advice = (
        "The flame is anchored, but its length is set by the domain, not the "
        "flame; rebuild with a larger outlet_length to measure it."
    )
    found = crossings(centre.x, centre_t, threshold)
    if not found:
        # Burnt gas behind the plate but a cold centreline all the way out: the
        # flame closes beyond the outlet.
        longer["message"] = (
            f"The slot centreline stays below {threshold:g} K up to the outlet "
            f"(x = {x_max:.4g} m) while the wake behind the plate is burnt: the "
            f"flame tip lies beyond the domain. {advice}"
        )
        return longer
    tip = min(found)
    # A tip pressed against the outlet is bounded by the domain just the same.
    if (x_max - tip) < 0.05 * (x_max - x_min):
        longer["flame_tip_x"] = tip
        longer["message"] = (
            f"The flame tip is at x = {tip:.4g} m, within 5% of the outlet "
            f"(x = {x_max:.4g} m), so the length is bounded by the domain. {advice}"
        )
        return longer
    return {
        **base,
        "status": "measured",
        "flame_length": tip - trailing_edge,
        "flame_tip_x": tip,
        "units": "m",
        "definition": (
            f"axial distance from the flame holder's trailing edge "
            f"(x = {trailing_edge:.4g} m) to the first point on the slot "
            f"centreline reaching {threshold:g} K"
        ),
    }


def case_type_of(case_dir: Path) -> Optional[str]:
    """Which template a case came from.

    Read from the provenance record the build wrote, falling back to the
    directory prefix for a case built before provenance existed.
    """
    record = provenance.read(Path(case_dir))
    if record:
        case_type = record.get("case", {}).get("type")
        if case_type:
            return str(case_type)
    import case_templates

    template = case_templates.template_for_case_name(Path(case_dir).name)
    return template.name if template else None


def analyse(
    case_dir: Path,
    threshold: float = FLAME_TEMPERATURE,
    sample: bool = True,
) -> dict[str, Any]:
    """Report how a finished run behaved and what its flame did.

    ``sample`` can be turned off to read only the log, which needs no OpenFOAM
    environment.
    """
    case_dir = Path(case_dir)
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case not found: {case_dir}"}

    case_type = case_type_of(case_dir)
    result: dict[str, Any] = {
        "case_name": case_dir.name,
        "case_type": case_type,
        "run": parse_log(case_dir),
    }

    if not sample:
        result["status"] = "log_only"
        return result

    import case_templates

    template = case_templates.TEMPLATE_BY_NAME.get(case_type or "")
    fields = list(template.sample_fields) if template else None
    samples, error = sample_axis(case_dir, fields=fields)
    if error:
        result["status"] = "log_only"
        result["sampling_error"] = error
        return result

    result["n_sampled_times"] = len(samples)
    if template is not None and template.analyse is not None:
        # The measurement is declared by the case plugin (agent/cases/<name>).
        result["flame"] = template.analyse(samples, threshold)
    else:
        result["flame"] = {
            "status": "no_measurement_defined",
            "message": (
                f"Case type {case_type!r} declares no measurement (its case.py "
                "sets no `analyse`), so only the run log is reported."
            ),
        }
    result["status"] = "analysed"
    return result
