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

"""The 2D laminar slot burner: an inlet-driven flame anchored on a plate.

Shows the fuller set of capabilities: extra parameters (the inlet and the
burner geometry), a parameterised mesh, a time-step velocity scale set by the
inlet, a domain length for the flow-through check, and a measurement relative
to a solid body (flame length, or flashback / blow-off).
"""

from case_api import CaseTemplate, ParamSpec
from measures import front_behind_body

# Streamwise lengths upstream (fixed) and downstream (`outlet_length`, this is
# its default) of the plate — `l_in` and `l_out` in
# testcases/laminarSlotBurner/system.orig/blockMeshDict.
L_IN = 6.2278e-03
L_OUT = 1.868340e-02
# dx = lF/10 in the same blockMeshDict.
CELLS_PER_FLAME = 10


def cell_size(params):
    thickness = params.get("flame_thickness")
    return float(thickness) / CELLS_PER_FLAME if thickness else None


def velocity_scale(params, thermo):
    """The burnt-gas velocity: gas leaves the flame rho_u/rho_b times faster."""
    inlet = params.get("inlet_velocity")
    ratio = getattr(thermo, "density_ratio", 0.0) or 0.0
    if not inlet or ratio <= 1.0:
        return None
    return float(inlet) * ratio, (
        f"inlet_velocity {float(inlet):g} m/s x rho_u/rho_b {ratio:.2f} "
        "(gas leaves the flame this much faster than it enters)"
    )


def domain_length(params):
    plate = params.get("plate_thickness")
    outlet = params.get("outlet_length") or L_OUT
    return None if plate is None else L_IN + float(plate) + float(outlet)


GEOMETRY = ("plate_thickness", "hole_pitch", "hole_size", "outlet_length",
            "flame_thickness")

def geometry_problem(params: dict) -> str | None:
    """Why this geometry cannot be meshed, or None.

    The plate's half-height is ``h_plate = pitch/2 - size/2`` (blockMeshDict), so
    a slot at least as wide as the pitch leaves no plate at all.
    """
    pitch, size = params.get("hole_pitch"), params.get("hole_size")
    lengths = {k: params.get(k) for k in GEOMETRY}
    negative = [k for k, v in lengths.items() if isinstance(v, (int, float)) and v <= 0]
    if negative:
        return "geometry must be positive: " + ", ".join(f"{k}={lengths[k]:g}" for k in negative)
    if isinstance(pitch, (int, float)) and isinstance(size, (int, float)) and size >= pitch:
        return (f"hole_size ({size:g} m) must be smaller than hole_pitch ({pitch:g} m): "
                "the plate between two slots is pitch - size wide, so this leaves "
                "no plate to stabilise the flame on")
    return None


CASE = CaseTemplate(
    name="laminarSlotBurner",
    summary=(
        "2D laminar slot burner — a premixed jet issuing between cooled "
        "plates, forming a steady Bunsen-type flame"
    ),
    required=("inlet_velocity", "T_inlet", "T_plate"),
    optional=("equivalence_ratio", "n_procs", "sim_write_steps", "p_0",
              "sim_time", "sim_timestep") + GEOMETRY,
    defaults={"sim_time": 0.1, "sim_timestep": 1.0e-5},
    extra_placeholders={
        "INLETTEMPERATURE": "T_inlet",
        "PLATETEMPERATURE": "T_plate",
        "INLETVELOCITY": "inlet_velocity",
        "FLAMETHICKNESS": "flame_thickness",
        "PLATETHICKNESS": "plate_thickness",
        "HOLEPITCH": "hole_pitch",
        "HOLESIZE": "hole_size",
        "OUTLETLENGTH": "outlet_length",
    },
    param_specs=(
        # Usual ranges (typical=): outside them the user is asked to confirm.
        # The laminar flame speed is ~0.19 m/s at the shipped table, and the
        # default outlet holds the flame tip only up to ~0.45 m/s.
        ParamSpec("inlet_velocity", float, True, None, "m/s", "inlet bulk velocity",
                  ("bulk velocity", "velocity"), typical=(0.05, 3.0)),
        ParamSpec("T_plate", float, True, None, "K", "plate wall temperature",
                  ("plate temperature", "wall temperature", "plate")),
        # Geometry: defaults reproduce the original 43000-cell mesh exactly.
        ParamSpec("flame_thickness", float, False, 6.2278e-04, "m",
                  "flame thickness lF; sets the cell size (dx = lF/10) and nothing "
                  "else, so it refines the mesh without changing the geometry",
                  ("flame thickness",), typical=(1.5e-4, 2.5e-3)),
        ParamSpec("plate_thickness", float, False, 3.1139e-03, "m",
                  "streamwise thickness of the plate (the flame holder)",
                  ("plate thickness", "plate", "thickness"), typical=(5e-4, 1.5e-2)),
        ParamSpec("hole_pitch", float, False, 1.245560e-02, "m",
                  "centre-to-centre distance between slots; the domain is one "
                  "half-slot of height hole_pitch/2", ("hole pitch", "slot pitch", "pitch"),
                  typical=(3e-3, 5e-2)),
        ParamSpec("hole_size", float, False, 7.473360e-03, "m",
                  "open height of the slot between plates",
                  ("hole size", "slot width", "slot size", "hole", "slot"),
                  typical=(1e-3, 4e-2)),
        ParamSpec("outlet_length", float, False, L_OUT, "m",
                  "length of the domain downstream of the plate; the flame tip "
                  "must fall inside it, so longer flames (higher inlet velocity) "
                  "need a longer outlet",
                  ("outlet length", "domain length", "outlet"), typical=(5e-3, 0.1)),
    ),
    sweepable=("T_inlet", "T_plate", "inlet_velocity") + GEOMETRY,
    geometric=GEOMETRY,
    axis_labels={
        "T_inlet": "Unburnt temperature  [K]",
        "T_plate": "Plate temperature  [K]",
        "inlet_velocity": "Inlet velocity  [m/s]",
        "plate_thickness": "Plate thickness  [m]",
        "hole_pitch": "Hole pitch  [m]",
        "hole_size": "Hole size  [m]",
        "outlet_length": "Outlet length  [m]",
        "flame_thickness": "Flame thickness (cell size)  [m]",
    },
    analyse=front_behind_body("T"),
    metric={"key": "flame_length", "label": "Flame length  [m]", "title": "Flame length"},
    cell_size=cell_size,
    velocity_scale=velocity_scale,
    domain_length=domain_length,
    check=geometry_problem,
    order=0,
    keywords=("slot", "burner", "2d", "two-dimensional", "flame length",
              "flashback", "blow-off", "blow off"),
)
