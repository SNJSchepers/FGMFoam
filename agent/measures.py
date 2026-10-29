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

"""Ready-made measurements for a case's ``analyse`` field.

Each function returns ``analyse(samples, threshold) -> dict``, where ``samples``
are fields sampled along the domain's lower x-wise edge at every written time
(by OpenFOAM's own ``postProcess``), each with the upper edge attached as
``sample.upper``. Most flame measurements are the same primitive
— *where a field crosses a threshold along a line* — seen differently:

=========================  ================================================
``front_speed``            crossing position tracked over time -> its speed
                           relative to the unburnt gas (laminar flame speed)
``front_behind_body``      crossing on the upper edge at the last time,
                           measured from the far edge of a solid body on the
                           lower edge (flame length), or flashback / blow-off
``crossing_position``      crossing at the last time, measured from x = 0
``field_extremum``         max or min of a field at the last time, and where
=========================  ================================================

For anything else, write your own ``analyse(samples, threshold)`` in the case's
``case.py`` — it receives the same samples.
"""

from __future__ import annotations

from typing import Any, Callable

Analyse = Callable[[list, float], dict]


def front_speed(field: str = "T", velocity_field: str = "U") -> Analyse:
    """Propagation speed of the threshold crossing, relative to the unburnt gas.

    The initial transient (e.g. a step-function ignition relaxing into a real
    flame) is detected and excluded. Result key: ``laminar_flame_speed``.
    """
    def analyse(samples: list, threshold: float) -> dict[str, Any]:
        import postrun
        return postrun.analyse_flamelet(samples, threshold, field, velocity_field)
    return analyse


def front_behind_body(field: str = "T") -> Analyse:
    """Distance from a solid body's trailing edge to the first crossing.

    The body is found as the gap it leaves in the lower line, so geometry changes
    need no code; the crossing is taken on the upper line (the slot centreline,
    where the flame closes). Reports ``flashback`` (crossing upstream of the
    body), ``blow_off`` and ``longer_than_domain`` instead of a length where one
    would be meaningless. Result key: ``flame_length``.
    """
    def analyse(samples: list, threshold: float) -> dict[str, Any]:
        import postrun
        return postrun.analyse_slot_burner(samples, threshold, field)
    return analyse


def crossing_position(field: str = "T") -> Analyse:
    """First position along x where ``field`` crosses the threshold, last time.

    Result key: ``crossing_position``.
    """
    def analyse(samples: list, threshold: float) -> dict[str, Any]:
        import postrun
        if not samples or field not in samples[-1].fields:
            return {"status": "insufficient_data", "message": f"{field} not sampled"}
        last = samples[-1]
        found = postrun.crossings(last.x, last.fields[field], threshold)
        if not found:
            return {"status": "no_crossing", "time": last.time, "threshold": threshold,
                    "message": f"{field} never crosses {threshold:g} along the axis"}
        return {"status": "measured", "time": last.time, "threshold": threshold,
                "crossing_position": found[0], "units": "m",
                "all_crossings": found}
    return analyse


def field_extremum(field: str = "T", kind: str = "max") -> Analyse:
    """Largest (or smallest) value of ``field`` along the axis at the last time.

    Result keys: ``<kind>_<field>`` and ``<kind>_<field>_position``.
    """
    if kind not in ("max", "min"):
        raise ValueError("kind must be 'max' or 'min'")

    def analyse(samples: list, threshold: float) -> dict[str, Any]:
        if not samples or field not in samples[-1].fields:
            return {"status": "insufficient_data", "message": f"{field} not sampled"}
        last = samples[-1]
        values = last.fields[field]
        pick = max if kind == "max" else min
        index = values.index(pick(values))
        return {"status": "measured", "time": last.time,
                f"{kind}_{field}": values[index],
                f"{kind}_{field}_position": last.x[index]}
    return analyse
