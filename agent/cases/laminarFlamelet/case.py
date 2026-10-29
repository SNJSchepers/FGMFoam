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

"""The 1D freely propagating flamelet: the minimal example of a case.

No parameters beyond the common ones, a fixed mesh, and a standard measurement
from measures.py. Start here when adding your own case.
"""

import cantera_tool
from case_api import CaseTemplate
from measures import front_speed

# Its blockMeshDict refines the flame region to a fixed 10 um (`nCellsRef`).
CELL_SIZE = 1.0e-5


def velocity_scale(params, thermo):
    """No inlet: only the flame moves the gas, so S_L sets the velocity scale.

    Needs a 1D free-flame solve (~10-20 s, cached on disk).
    """
    if (getattr(thermo, "density_ratio", 0.0) or 0.0) <= 1.0:
        return None
    flame = cantera_tool.free_flame(
        float(params.get("equivalence_ratio", 0.7)),
        float(params["T_inlet"]),
        float(params.get("p_0", 101325.0)),
    )
    if flame is None:
        return None
    return flame.flame_speed * flame.density_ratio, (
        f"laminar flame speed {flame.flame_speed:.4g} m/s x rho_u/rho_b "
        f"{flame.density_ratio:.2f} (burnt gas driven out by expansion), "
        "from a 1D free-flame solve"
    )


CASE = CaseTemplate(
    name="laminarFlamelet",
    summary=(
        "1D freely propagating laminar flame — a quiescent premixed mixture "
        "ignited in a box region; no inlet and no walls, so only the initial "
        "conditions and the run controls apply"
    ),
    required=("T_inlet",),
    optional=("equivalence_ratio", "n_procs", "sim_write_steps", "p_0",
              "sim_time", "sim_timestep"),
    defaults={"sim_time": 0.1, "sim_timestep": 2.0e-6},
    extra_placeholders={"INITIALTEMPERATURE": "T_inlet"},
    sweepable=("T_inlet",),
    axis_labels={"T_inlet": "Unburnt temperature  [K]"},
    analyse=front_speed("T", "U"),
    metric={"key": "laminar_flame_speed", "label": "Laminar flame speed  [m/s]",
            "title": "Laminar flame speed"},
    cell_size=lambda params: CELL_SIZE,
    velocity_scale=velocity_scale,
    order=1,
    keywords=("flamelet", "1d", "one-dimensional", "free flame", "freely propagating",
              "flame speed", "burning velocity"),
)
