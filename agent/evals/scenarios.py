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

"""Evaluation scenarios for the FGMFoam agent: simulation manager + solver Q&A.

Each scenario is one fresh conversation. A turn is a user message plus the checks
its outcome must satisfy. Checks are deterministic — they look at the tools the
model actually executed (and what those tools returned), at the case files the
tools wrote into the sandbox, and at the reply text. Nothing here is graded by
another LLM.

Check vocabulary (all keys optional):

  tool       {"name", "args"?, "status"?}  — at least one EXECUTED call of that
             tool (or of any tool, when name is a list) whose args match (numbers within 1 %, strings case-insensitive,
             parameter aliases resolved) and whose result status matches.
  no_status  [(tool, status), ...]         — no executed call of that tool may
             return that status (e.g. no build_case -> "built").
  no_tool    [name, ...]                   — tool must not be executed. A call
             blocked by one of chat.py's guards counts as not executed.
  case       {param: value}                — the case most recently built in this
             scenario carries these parameters in agent-provenance.json.
  reply_all  [regex | [regex, ...], ...]   — every item must match the reply; a
             list item is an any-of group.
  reply_none [regex, ...]                  — none may match (hallucination traps).
  grounded   regex                         — a file whose path matches was read
             this turn (by the model, or by the deterministic forced read).
  model_read regex                         — as grounded, but only reads the MODEL
             chose count. Reported, and scored as its own column.

Ground truth for the Q&A checks was taken from the solver source (file:line in
each ``why``) — re-check them after changing the solver.
"""

from __future__ import annotations

SLOT_BURNER = "laminarSlotBurner"
FLAMELET = "laminarFlamelet"

# --------------------------------------------------------------------------- #
# A. Simulation manager                                                       #
# --------------------------------------------------------------------------- #

SIMULATION: list[dict] = [
    {
        "id": "S01_flamelet_build",
        "title": "Build a flamelet from one condition",
        "turns": [{
            "user": "Set up a laminar flamelet case with an unburnt gas temperature of 400 K.",
            "tool": {"name": "build_case", "args": {"case_type": FLAMELET, "T_inlet": 400},
                     "status": "built"},
            "case": {"T_inlet": 400},
            "no_tool": ["run_case", "run_cases"],
        }],
    },
    {
        "id": "S02_slotburner_build",
        "title": "Build a slot burner from three conditions",
        "turns": [{
            "user": ("Set up a slot burner case with an inlet velocity of 0.4 m/s, "
                     "an inlet temperature of 300 K and a plate temperature of 400 K."),
            "tool": {"name": "build_case",
                     "args": {"inlet_velocity": 0.4, "T_inlet": 300, "T_plate": 400},
                     "status": "built"},
            "case": {"inlet_velocity": 0.4, "T_inlet": 300, "T_plate": 400},
            "no_tool": ["run_case", "run_cases"],
        }],
    },
    {
        "id": "S03_ambiguous_case_type",
        "title": "Ask which case when the type is unclear",
        "turns": [{
            "user": "Set up a simulation with an unburnt temperature of 350 K.",
            "no_status": [("build_case", "built")],
            "reply_all": [r"(?i)flamelet|\b1D\b|freely propagating",
                          r"(?i)slot|\b2D\b|burner"],
        }],
    },
    {
        "id": "S04_missing_parameters",
        "title": "Ask for missing slot burner parameters instead of inventing them",
        "turns": [{
            "user": "Set up a slot burner case with an inlet velocity of 0.5 m/s.",
            "no_status": [("build_case", "built")],
            "reply_all": [r"(?i)temperature|T_inlet|T_plate"],
        }],
    },
    {
        "id": "S05_geometry_change",
        "title": "Change the burner geometry of a built case",
        "turns": [
            {
                "user": ("Set up a slot burner case named geo with inlet velocity 0.4 m/s, "
                         "inlet temperature 300 K and plate temperature 400 K."),
                "tool": {"name": "build_case", "status": "built"},
            },
            {
                "user": "Now rebuild that case with a plate thickness of 5 mm.",
                "tool": {"name": ["build_case", "change_settings"], "args": {"plate_thickness": 0.005},
                         "status": "built"},
                # The conditions must survive the geometry change.
                "case": {"plate_thickness": 0.005, "inlet_velocity": 0.4,
                         "T_inlet": 300, "T_plate": 400},
            },
            {
                "user": "And make the slot (hole size) 6 mm wide, keeping everything else.",
                "tool": {"name": ["build_case", "change_settings"], "args": {"hole_size": 0.006},
                         "status": "built"},
                "case": {"hole_size": 0.006, "plate_thickness": 0.005,
                         "inlet_velocity": 0.4},
            },
        ],
    },
    {
        "id": "S06_condition_change",
        "title": "Change an operating condition and a run setting",
        "turns": [
            {
                "user": "Set up a laminar flamelet case named cond at 300 K.",
                "tool": {"name": "build_case", "status": "built"},
            },
            {
                "user": "Rebuild it with an unburnt temperature of 450 K instead.",
                "tool": {"name": ["build_case", "change_settings"], "args": {"T_inlet": 450}, "status": "built"},
                "case": {"T_inlet": 450},
            },
            {
                "user": "Change the simulation end time of that case to 0.02 s.",
                "tool": {"name": "change_settings", "args": {"sim_time": 0.02}},
                "no_status": [("change_settings", "error")],
            },
        ],
    },
    {
        "id": "S07_build_then_run_same_turn",
        "title": "Build and run requested together: build only, run is a separate step",
        "turns": [{
            "user": "Build a laminar flamelet case at 300 K and run it.",
            "tool": {"name": "build_case", "status": "built"},
            "no_tool": ["run_case", "run_cases"],
        }],
    },
    {
        "id": "S08_campaign_cpu_plan",
        "title": "Campaign: ask the CPU plan, then prepare without launching",
        "turns": [
            {
                "user": ("Run a campaign of laminar flamelets at unburnt temperatures "
                         "of 300, 400 and 500 K."),
                "tool": {"name": "run_campaign",
                         # "variable" may be left out: T_inlet is the flamelet's
                         # only sweepable parameter, so the tool infers it.
                         "args": {"case_type": FLAMELET, "values": [300, 400, 500]},
                         "status": "missing_cpu_plan"},
                "no_status": [("run_campaign", "prepared"), ("run_campaign", "started")],
                "reply_all": [r"(?i)core|processor|cpu|proc"],
            },
            {
                "user": ("Use 2 processors per case and run 2 cases at the same time. "
                         "Only prepare it, do not start it yet."),
                "tool": {"name": "run_campaign",
                         "args": {"n_procs_per_case": 2, "n_concurrent_cases": 2,
                                  "launch": False},
                         "status": "prepared"},
                "no_tool": ["start_campaign", "run_cases"],
            },
        ],
    },
    {
        "id": "S09_geometry_campaign",
        "title": "Campaign over a geometry parameter (plate thickness)",
        "turns": [{
            "user": ("Prepare, but do not start, a slot burner campaign that sweeps the plate "
                     "thickness over 2, 3 and 4 mm, at inlet velocity 0.4 m/s, inlet "
                     "temperature 300 K and plate temperature 400 K. Use 4 processors per "
                     "case and 2 cases concurrently."),
            "tool": {"name": "run_campaign",
                     "args": {"case_type": SLOT_BURNER, "variable": "plate_thickness",
                              "values": [0.002, 0.003, 0.004], "launch": False},
                     "status": "prepared"},
            "no_tool": ["start_campaign", "run_cases"],
        }],
    },
    {
        "id": "S10_refused_sweep",
        "title": "Refuse a sweep the table cannot support (equivalence ratio)",
        "turns": [{
            "user": ("Run a flamelet campaign over equivalence ratios 0.6, 0.7 and 0.8, "
                     "with 2 processors per case and 2 cases concurrently."),
            "no_status": [("run_campaign", "prepared"), ("run_campaign", "started"),
                          ("build_case", "built")],
            "reply_all": [r"(?i)table|fixed|0\.7"],
        }],
    },
    {
        "id": "S11_feasibility_pressure",
        "title": "Answer a feasibility question instead of acting (5 bar)",
        "turns": [{
            "user": "Is it possible to run the slot burner at 5 bar?",
            "no_status": [("build_case", "built")],
            "reply_all": [r"(?i)table", [r"(?i)101325|1 ?atm|atmospheric|1 bar|1\.01"]],
            "reply_none": [r"(?i)^yes\b"],
        }],
    },
]

# --------------------------------------------------------------------------- #
# B. Solver information                                                       #
# --------------------------------------------------------------------------- #

QA: list[dict] = [
    {
        "id": "Q01_lookup_enthalpy_bc",
        "title": "lookUpEnthalpy boundary condition",
        "why": ("lookUpEnthalpyFvPatchScalarField.C:339-352: iterates ht += 0.1*cp*(T0-T1) "
                "until |T0-T1| < 1 K, at most 1000 iterations"),
        "turns": [{
            "user": "How does the lookUpEnthalpy boundary condition work?",
            "grounded": r"lookUpEnthalpy",
            "reply_all": [r"(?i)temperature", r"(?i)enthalpy", r"(?i)iterat",
                          [r"1000", r"(?i)1 ?K\b|Tresidual|residual"]],
            "reply_none": [r"(?i)neural", r"(?i)from the mesh file"],
        }],
    },
    {
        "id": "Q02_table_sharing_mpi",
        "title": "How the table is shared in parallel",
        "why": "FGMTable.C:89-148: MPI_Comm_split_type + MPI_Win_allocate_shared, one copy per node",
        "turns": [{
            "user": "How is the FGM table shared between MPI ranks when running in parallel?",
            "grounded": r"FGMTable",
            "reply_all": [r"(?i)shared[- ]memory|MPI_Win_allocate_shared",
                          r"(?i)node|MPI_COMM_TYPE_SHARED|same (machine|host)"],
            "reply_none": [r"MPI_Bcast", r"(?i)each (rank|process) (reads|loads) (its own|the (whole )?(file|table))"],
        }],
    },
    {
        "id": "Q03_out_of_range_lookup",
        "title": "Lookup outside the table range (trap: it extrapolates)",
        "why": ("FGMlib.c:220-253: cell index clamped to the edge cell but the weights are "
                "not, so the result is linear extrapolation; always returns EXIT_SUCCESS"),
        "turns": [{
            "user": ("Explain what the table lookup does when Yc or ht falls outside "
                     "the range of the table."),
            "grounded": r"FGMlib",
            "reply_all": [r"(?i)extrapolat"],
            # Claiming it aborts is wrong; "no error, warning, or abort" is right.
            "reply_none": [r"FatalError",
                           r"(?i)\b(?:aborts|will abort|is aborted|abort(?:s|ing)? (?:the )?(?:run|simulation|solver))",
                           r"(?i)(ranges?|spans?|goes|varies) from 0 to 1"],
        }],
    },
    {
        "id": "Q04_combustion_properties",
        "title": "Valid combustionProperties options",
        "why": "FGM.C:63-90: only premixed / unityLewis / no / none, anything else FatalError",
        "turns": [{
            "user": ("Describe which options I can set in constant/combustionProperties "
                     "for the FGM model."),
            "grounded": r"FGM\.[CH]|combustionProperties",
            "reply_all": [r"premixed", r"unityLewis", r"(?i)\bnone\b", r"(?i)adiabatic"],
            "reply_none": [r"\bEDC\b", r"PaSR", r"(?i)nonPremixed is supported"],
        }],
    },
    {
        "id": "Q05_update_fields",
        "title": "What update() writes back into the thermo",
        "why": ("FGM_PM_UL_NA.C:366-460: T, thermo:psi, mu, alpha, species, HeatRelease; "
                "negative sourceYc clipped to 0"),
        "turns": [{
            "user": ("Walk me through what the FGM model's update step writes into the "
                     "flow fields, and what it does with negative source terms."),
            "grounded": r"FGM_PM_UL_NA",
            "reply_all": [r"\bT\b|(?i:temperature)", r"(?i)psi|compressibility",
                          r"(?i)\bmu\b|viscosity", [r"(?i)\balpha\b|conductivity|diffusivity"],
                          [r"(?i)(set|clip|clamp|limit)\w* to (zero|0)|to 0(\.0)?\b|= ?0"]],
        }],
    },
    {
        "id": "Q06_progress_variable",
        "title": "Progress variable definition (trap: not normalised)",
        "why": "knowledge/solver_notes.md:359: Yc = Y_CO2 / M_CO2, not normalised to 0-1",
        "turns": [{
            "user": "How is the progress variable Yc defined in FGMFoam?",
            "reply_all": [r"CO2|CO_2|CO₂", r"(?i)molar mass|M_CO2|mole|molecular weight"],
            "reply_none": [r"(?i)normali[sz]ed (between|from|to) 0",
                           r"(?i)(ranges?|goes|varies) from 0 to 1"],
        }],
    },
    {
        "id": "Q07_equivalence_ratio",
        "title": "Other equivalence ratio (table limitation)",
        "why": "solver_notes.md: shipped table is phi = 0.7, p = 101325 Pa",
        "turns": [{
            "user": "Can I run a simulation at an equivalence ratio of 1.0 with the shipped table?",
            "no_status": [("build_case", "built")],
            "reply_all": [r"0\.7", r"(?i)regenerat|new table|another table|different table|generate"],
            "reply_none": [r"(?i)^yes\b"],
        }],
    },
    {
        "id": "Q08_transport_equations",
        "title": "How Yc and ht are transported",
        "why": ("FGM_PM_UL_NA.C solve(): Yc and ht equations with a shared multivariate "
                "convection scheme; sourceYc from the table"),
        "turns": [{
            "user": "How does FGMFoam solve the transport equations for Yc and ht?",
            "grounded": r"FGM_PM_UL_NA|FGMFoam",
            "reply_all": [r"(?i)source", r"(?i)table|look ?up"],
            "reply_none": [r"(?i)neural", r"\bLTE\b", r"(?i)reaction mechanism is integrated"],
        }],
    },
    {
        "id": "Q09_table_column_names",
        "title": "Table columns are matched by name",
        "why": "FGM_PM_UL_NA.C:240-280: std::strcmp against SOURCE_CV1, TEMPERATURE, CP, ...",
        "turns": [{
            "user": ("Explain how the solver finds the right columns in the FGM table, "
                     "and what happens if I rename a column."),
            "grounded": r"FGM_PM_UL_NA|FGMlib",
            "reply_all": [r"(?i)name", [r"SOURCE_CV1", r"TEMPERATURE", r"strcmp"]],
        }],
    },
    {
        "id": "Q10_table_path",
        "title": "Where the table is read from (trap: hardcoded)",
        "why": "FGM_PM_UL_NA.C: FGMTable::getFGMInstance(\"constant/lookUp/database.fgm\")",
        "turns": [{
            "user": ("Describe where the solver reads the FGM table from. Can I point it to a "
                     "different file in combustionProperties?"),
            "reply_all": [r"constant/lookUp/database\.fgm",
                          r"(?i)hard[- ]?coded|fixed|cannot|can't|not configurable|no\b"],
            "reply_none": [r"(?i)(tableFile|tablePath|tableName)\s"],
        }],
    },
]

SUITES: dict[str, list[dict]] = {"sim": SIMULATION, "qa": QA}
