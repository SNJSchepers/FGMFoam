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

"""Thermochemical calculations for the laminar slot burner setup.

Given an equivalence ratio and an unburnt inlet temperature, Cantera is used to
derive the quantities that the OpenFOAM/FGM case needs:

  * the progress variable in the fresh (unburnt) mixture,
  * the progress variable in the fully burnt (chemical equilibrium) mixture, and
  * the specific enthalpy of the unburnt mixture.

The progress variable here is the specific mole number of CO2,
Yc = Y_CO2 / M_CO2 (mass fraction divided by molar mass), consistent with the
FGM table tabulation used by the solver.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import cantera as ct

# Reaction mechanism and mixture definition for the methane/air flame.
MECHANISM: str = "gri30.yaml"
FUEL: str = "CH4"
OXIDIZER: str = "O2:1, N2:3.76"

# Species whose mass fraction defines the FGM progress variable Yc.
PROGRESS_SPECIES: str = "CO2"


@dataclass
class Thermochemistry:
    """Derived thermochemical state used to fill in the OpenFOAM fields."""

    Yc_u: float  # progress variable Y_CO2/M_CO2 of the unburnt mixture [kmol/kg]
    Yc_b: float  # progress variable Y_CO2/M_CO2 at equilibrium [kmol/kg]
    ht: float    # specific enthalpy of the unburnt mixture [J/kg]
    # Thermal expansion across the flame, rho_u/rho_b. The gas leaves the flame
    # this many times faster than it enters, which is what sets the Courant
    # number — so the time-step estimate needs it. Free from the equilibrium
    # calculation already done for Yc_b.
    density_ratio: float = 1.0


def compute_thermochemistry(
    equivalence_ratio: float,
    t_inlet: float,
    p_0: float,
) -> Thermochemistry:
    """Compute progress-variable bounds and unburnt enthalpy with Cantera.

    Args:
        equivalence_ratio: Fuel/air equivalence ratio (phi) [-].
        t_inlet: Temperature of the unburnt inlet mixture [K].
        p_0: Operating/initial pressure [Pa].

    Returns:
        A populated :class:`Thermochemistry` instance.
    """
    gas = ct.Solution(MECHANISM)

    # Set the unburnt composition and thermodynamic state.
    gas.set_equivalence_ratio(equivalence_ratio, fuel=FUEL, oxidizer=OXIDIZER)
    gas.TP = t_inlet, p_0

    species_idx = gas.species_index(PROGRESS_SPECIES)
    # Molar mass of the progress-variable species [kg/kmol]; constant.
    molar_mass = gas.molecular_weights[species_idx]

    # Unburnt progress variable and specific enthalpy (captured before burning).
    yc_u = float(gas.Y[species_idx] / molar_mass)
    ht = float(gas.enthalpy_mass)

    density_unburnt = float(gas.density)

    # Burn to chemical equilibrium at constant enthalpy and pressure (adiabatic).
    gas.equilibrate("HP")
    yc_b = float(gas.Y[species_idx] / molar_mass)
    density_ratio = density_unburnt / float(gas.density)

    return Thermochemistry(
        Yc_u=yc_u, Yc_b=yc_b, ht=ht, density_ratio=density_ratio
    )


# --------------------------------------------------------------------------- #
# Freely propagating flame (only needed where the flame sets the velocity)
# --------------------------------------------------------------------------- #

# Solving a 1D flame takes ~10-20 s, so results are cached on disk: a campaign
# rebuilt or re-analysed repeatedly pays for each condition only once.
_FLAME_CACHE: Path = Path(__file__).resolve().parent / ".cache" / "free_flames.json"


@dataclass
class FreeFlame:
    """Laminar flame speed and thickness of a freely propagating flame."""

    flame_speed: float      # S_L, relative to the unburnt gas [m/s]
    flame_thickness: float  # thermal thickness (T_b - T_u)/max(dT/dx) [m]
    density_ratio: float    # rho_u/rho_b [-]


# Transport for the 1D flame. Unity Lewis and no thermal diffusion (Soret), to
# match the FGM model the solver uses (FGM_PM_UL_NA: unity Lewis). A reference
# with mixture-averaged transport would include differential diffusion that the
# table deliberately leaves out, and the comparison would measure that instead.
FLAME_TRANSPORT: str = "unity-Lewis-number"


def _cache_key(equivalence_ratio: float, t_unburnt: float, p_0: float) -> str:
    # Transport is part of the key: results for another model are not reusable.
    return (f"{MECHANISM}|{FLAME_TRANSPORT}|nosoret|{equivalence_ratio:.6g}|"
            f"{t_unburnt:.6g}|{p_0:.6g}")


def free_flame(
    equivalence_ratio: float, t_unburnt: float, p_0: float,
) -> Optional[FreeFlame]:
    """Solve a 1D freely propagating flame (cached), or None if it fails.

    Needed for the flamelet case, where nothing but the flame itself sets the
    velocity scale. Returns None rather than raising: a time-step estimate is a
    convenience, and failing to produce one must never fail a build.
    """
    key = _cache_key(equivalence_ratio, t_unburnt, p_0)
    try:
        cached = json.loads(_FLAME_CACHE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        cached = {}
    if key in cached:
        entry = cached[key]
        return FreeFlame(**entry)

    try:
        import numpy as np

        gas = ct.Solution(MECHANISM)
        gas.set_equivalence_ratio(equivalence_ratio, fuel=FUEL, oxidizer=OXIDIZER)
        gas.TP = t_unburnt, p_0
        density_unburnt = float(gas.density)
        flame = ct.FreeFlame(gas, width=0.03)
        flame.transport_model = FLAME_TRANSPORT
        flame.soret_enabled = False
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.2)
        flame.solve(loglevel=0, auto=True)
        thickness = float(
            (flame.T[-1] - flame.T[0]) / np.max(np.gradient(flame.T, flame.grid))
        )
        result = FreeFlame(
            flame_speed=float(flame.velocity[0]),
            flame_thickness=thickness,
            density_ratio=density_unburnt / float(flame.density[-1]),
        )
    except Exception:  # noqa: BLE001 - never fail a build over an estimate
        return None

    cached[key] = {
        "flame_speed": result.flame_speed,
        "flame_thickness": result.flame_thickness,
        "density_ratio": result.density_ratio,
    }
    try:
        _FLAME_CACHE.parent.mkdir(parents=True, exist_ok=True)
        _FLAME_CACHE.write_text(json.dumps(cached, indent=2), encoding="utf-8")
    except OSError:
        pass
    return result
