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

"""The plugin API: everything needed to teach the agent a new simulation case.

A case is a folder under ``agent/cases/`` (or on ``FGM_CASE_PATH``) holding a
``case.py`` that defines ``CASE = CaseTemplate(...)``. The folder may also hold
the OpenFOAM case itself and optional ``*.md`` notes that are indexed into the
agent's knowledge. Nothing else in the agent needs editing: provenance, the
pre-run checks, time-step selection, post-run measurement, campaigns, CPU
planning and graphs all work from this declaration.

See ``agent/cases/laminarFlamelet/case.py`` (the minimal example) and
``agent/cases/laminarSlotBurner/case.py`` (geometry, custom measurement).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

REPO_ROOT: Path = Path(__file__).resolve().parent.parent
TESTCASES_DIR: Path = REPO_ROOT / "testcases"

# Values available to every template's placeholders without declaring them:
# the tabulated thermochemistry (from Cantera) and the run controls.
DERIVED_VALUES: frozenset[str] = frozenset({"Yc_u", "Yc_b", "ht", "write_time"})

COMMON_PLACEHOLDERS: dict[str, str] = {
    "INITIALPROGRESSVARIABLE": "Yc_u",
    "FINALPROGRESSVARIABLE": "Yc_b",
    "INITIALENTHALPY": "ht",
    "INITIALPRESSURE": "p_0",
    "SIMULATIONTIME": "sim_time",
    "SIMULATIONDELTATIME": "sim_timestep",
    "SIMULATIONWRITETIME": "write_time",
    "NPROCS": "n_procs",
}


@dataclass(frozen=True)
class ParamSpec:
    """One simulation parameter: its name, type, default, unit and meaning."""

    name: str               # internal/JSON key
    py_type: type           # float or int
    required: bool          # informational; templates decide what is required
    default: Optional[Any]  # used when the parameter is optional and not given
    unit: str               # physical unit, for prompts and the summary
    description: str        # shown to the model and in the parameter table
    # How people say it ("plate thickness", "slot width"). With a number and a
    # matching unit next to it in the user's message ("a plate thickness of
    # 5 mm"), the chat loop makes sure the tool call carries that value — small
    # models drop or mangle stated values. The name with spaces always counts.
    phrases: tuple[str, ...] = ()
    # The range a value of this parameter usually has in this case. A value
    # outside it is not refused but must be confirmed by the user first, with a
    # hint when it looks like a unit slip (40 m/s that should be 0.40, cm/s).
    typical: Optional[tuple[float, float]] = None


@dataclass(frozen=True)
class CaseTemplate:
    """A simulation case the agent can build, run, check, measure and sweep.

    Only ``name``, ``summary``, ``required``, ``optional`` and
    ``extra_placeholders`` are needed for a case that can be built and run.
    Each further field switches on one more capability.
    """

    name: str                        # template directory name; no underscores
    summary: str                     # one line, shown to the user and the model
    required: tuple[str, ...]        # parameters the user must give
    optional: tuple[str, ...]        # parameters that fall back to a default
    # Placeholder token in the case files -> parameter (or derived value).
    extra_placeholders: dict[str, str]
    # Recommended values overriding the generic defaults, e.g. run controls.
    defaults: dict[str, float] = field(default_factory=dict)
    reference_temperature: float = 300.0
    # Case-relative files that contain placeholders; missing ones are skipped.
    target_files: tuple[str, ...] = (
        "0/T", "0/U", "0/Yc", "0/ht", "0/p",
        "system/controlDict", "system/decomposeParDict", "system/setFieldsDict",
        "system/blockMeshDict",
    )
    # Where the OpenFOAM case lives; default testcases/<name>/.
    template_path: Optional[Path] = None
    # Parameters this case introduces beyond the common ones.
    param_specs: tuple[ParamSpec, ...] = ()
    # --- campaigns ------------------------------------------------------------
    sweepable: tuple[str, ...] = ()
    geometric: tuple[str, ...] = ()          # sweeping these rebuilds the mesh
    axis_labels: dict[str, str] = field(default_factory=dict)
    # --- post-run measurement -------------------------------------------------
    # analyse(samples, threshold) -> dict; build one from measures.py.
    analyse: Optional[Callable[[list, float], dict]] = None
    # {"key": result field to plot, "label": y-axis label, "title": graph title}
    metric: Optional[dict[str, str]] = None
    sample_fields: tuple[str, ...] = ("T", "Yc", "U")
    # --- time step: dt = CFL * cell_size(params) / velocity_scale(...) ---------
    cell_size: Optional[Callable[[dict], Optional[float]]] = None
    # velocity_scale(params, thermo) -> (velocity [m/s], explanation) or None
    velocity_scale: Optional[Callable[[dict, Any], Optional[tuple[float, str]]]] = None
    # --- pre-run checks -------------------------------------------------------
    domain_length: Optional[Callable[[dict], Optional[float]]] = None
    # check(params) -> why these values cannot make a valid case, or None. Hard:
    # runs on the complete parameter set (defaults applied) before anything is
    # written, e.g. a slot wider than the pitch leaves no plate to mesh.
    check: Optional[Callable[[dict], Optional[str]]] = None
    # Listing order; the lowest is the default case.
    order: int = 100
    # Words that show the USER asked for this case ("flamelet", "slot burner").
    # A small model sometimes picks a case type the user never named; the chat
    # loop only accepts a case type the conversation mentions (by its name, one of
    # these words, or a parameter unique to it), and otherwise asks.
    keywords: tuple[str, ...] = ()

    @property
    def template_dir(self) -> Path:
        return Path(self.template_path) if self.template_path else TESTCASES_DIR / self.name

    @property
    def prefix(self) -> str:
        """Generated cases are ``<template>_<descriptor>``, never the bare name."""
        return f"{self.name}_"

    @property
    def default_case_name(self) -> str:
        return f"{self.name}_agent"

    def mentioned_in(self, text: str) -> bool:
        """True when ``text`` names this case: its template name or a keyword."""
        lowered = text.lower()
        return any(word.lower() in lowered for word in (self.name, *self.keywords))

    @property
    def placeholders(self) -> dict[str, str]:
        return {**COMMON_PLACEHOLDERS, **self.extra_placeholders}

    @property
    def parameters(self) -> tuple[str, ...]:
        return self.required + self.optional


# Parameters every case has. A case declares only what it adds (param_specs).
COMMON_PARAM_SPECS: tuple[ParamSpec, ...] = (
    ParamSpec("equivalence_ratio", float, False, 0.7, "-",
              "fuel/air equivalence ratio (phi); fixed at 0.7 by the current FGM "
              "table — a different value needs a regenerated table"),
    ParamSpec("T_inlet", float, True, None, "K",
              "unburnt-gas temperature (inlet, or initial for a quiescent mixture)",
              ("inlet temperature", "unburnt temperature", "unburnt gas temperature",
               "unburnt-gas temperature", "initial temperature")),
    ParamSpec("n_procs", int, False, 4, "-",
              "number of MPI processes / mesh subdomains"),
    ParamSpec("sim_time", float, True, None, "s",
              "total simulated physical time (endTime)",
              ("end time", "simulation time", "simulated time", "endtime")),
    ParamSpec("sim_timestep", float, True, None, "s",
              "fixed solver time step (deltaT)", ("time step", "timestep", "deltat")),
    ParamSpec("sim_write_steps", int, False, 100, "-",
              "number of write snapshots over the run"),
    ParamSpec("p_0", float, False, 101325.0, "Pa",
              "operating/initial pressure; must match the pressure at which the "
              "FGM table was generated"),
)
COMMON_PARAMETERS: frozenset[str] = frozenset(s.name for s in COMMON_PARAM_SPECS)
