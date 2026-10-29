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

"""Choose the fixed time step from the mesh and the flow, before the run.

The solver runs at a fixed step (`adjustTimeStep no`, and `maxCo` in controlDict
is inert), so nothing corrects a bad choice at run time. The step therefore has
to be right up front:

    dt = CFL x dx / U

where ``dx`` is the smallest cell the flame will cross and ``U`` is the fastest
velocity that matters. The velocity is **not** the inlet velocity: gas leaves a
flame roughly ``rho_u/rho_b`` times faster than it enters (about 6x for
methane-air at phi = 0.7), and it is that burnt-gas velocity that sets the
Courant number.

Each case declares both in its plugin (``cell_size`` and ``velocity_scale`` in
``agent/cases/<name>/case.py``). The two shipped cases:

* **laminarSlotBurner** — the flow is driven by the inlet, so
  ``U = inlet_velocity x rho_u/rho_b``. The density ratio comes free from the
  equilibrium calculation already done for the progress variable, so this
  estimate costs nothing.
* **laminarFlamelet** — there is no inlet; the flame itself sets the velocity
  scale, so ``U = S_L x rho_u/rho_b``. That needs a 1D free-flame solve
  (~10-20 s, cached on disk).

Checked against real runs: for the slot burner this predicts dt = 5.07e-6 s at
the shipped conditions, and a run at 5e-6 s held a sustained Courant of 0.11
against the 0.1 target. For the flamelet at 2e-6 s it predicts 0.21 and the run
held 0.19.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Optional


# Target Courant number. 0.1 is deliberately well below the nominal stability
# limit: the step cannot adapt, the estimate is only as good as its velocity
# scale, and an ignition transient briefly exceeds it in any case.
TARGET_CFL: float = 0.1


@dataclass(frozen=True)
class Estimate:
    """A recommended time step and the reasoning behind it."""

    timestep: float
    cell_size: float
    velocity: float
    velocity_source: str
    cfl: float
    explanation: str


def _round_to_1_sigfig(value: float) -> float:
    """1 significant figure, so the reported step is a round number."""
    if value <= 0:
        return value
    magnitude = 10.0 ** math.floor(math.log10(value))
    return round(value / magnitude) * magnitude


def _snap_to_write_interval(timestep: float, write_interval: float) -> float:
    """Shrink the step so it divides the write interval a whole number of times.

    Otherwise the solver lands on 0.019998 instead of 0.02 and the written times
    drift away from the requested snapshots.
    """
    if write_interval <= 0 or timestep <= 0:
        return timestep
    steps_per_write = max(1, round(write_interval / timestep))
    return write_interval / steps_per_write


def _template(case_type: str):
    import case_templates
    return case_templates.TEMPLATE_BY_NAME.get(case_type)


def cell_size(case_type: str, params: dict[str, Any]) -> Optional[float]:
    """Smallest cell the flame crosses — from the case plugin's ``cell_size``."""
    template = _template(case_type)
    if template is None or template.cell_size is None:
        return None
    return template.cell_size(params)


def estimate(
    case_type: str,
    params: dict[str, Any],
    thermo: Any,
    cfl: float = TARGET_CFL,
) -> Optional[Estimate]:
    """Recommend a time step, or None when the case gives no velocity scale.

    The mesh and the velocity scale come from the case plugin's ``cell_size``
    and ``velocity_scale`` hooks. Never raises: a missing estimate falls back
    to the case's shipped value.
    """
    template = _template(case_type)
    if template is None or template.velocity_scale is None:
        return None
    dx = cell_size(case_type, params)
    if not dx or dx <= 0:
        return None
    try:
        scale = template.velocity_scale(params, thermo)
    except Exception:  # noqa: BLE001 - an estimate must never fail a build
        return None
    if not scale:
        return None
    velocity, source = scale
    if velocity <= 0:
        return None

    raw = _round_to_1_sigfig(cfl * dx / velocity)
    sim_time = params.get("sim_time")
    write_steps = params.get("sim_write_steps")
    timestep = raw
    snapped = ""
    if sim_time and write_steps:
        write_interval = float(sim_time) / float(write_steps)
        timestep = _snap_to_write_interval(raw, write_interval)
        if timestep != raw:
            snapped = (
                f", then snapped to {write_interval:.6g}/"
                f"{round(write_interval / timestep)} so it divides the write "
                "interval exactly"
            )

    return Estimate(
        timestep=timestep,
        cell_size=dx,
        velocity=velocity,
        velocity_source=source,
        cfl=cfl,
        explanation=(
            f"dt = CFL {cfl:g} x dx {dx:.3e} m / U {velocity:.3f} m/s "
            f"= {raw:.2e} s{snapped}. U is the {source}. The solver runs at a "
            "fixed step, so this is chosen before the run rather than adapted "
            "during it."
        ),
    )
