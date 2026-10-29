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

"""Tools the conversational agent can call.

Each tool wraps the existing deterministic pipeline (`validator`,
`cantera_tool`, `of_writer`). The LLM decides *when* to call them; the actual
work and all validation stay in Python, so case building remains reproducible
and the model cannot bypass the required-parameter checks.

`TOOLS` holds the Mistral function-calling schemas; `execute_tool` dispatches a
call by name and always returns a JSON-serialisable dict (including on error, so
the model can react conversationally).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

import cantera_tool
import campaign as campaign_module
import case_templates
import config
import fgm_table
import of_writer
import postrun
import precheck
import providers
import provenance
import runner
import timestep
import validator
from config import REPO_ROOT
from rag import corpus

# The shipped FGM table bakes in a fixed equivalence ratio and pressure; a
# different value is not supported without regenerating the table.
TABLE_EQUIVALENCE_RATIO: float = validator.SPEC_BY_NAME["equivalence_ratio"].default
TABLE_PRESSURE: float = validator.SPEC_BY_NAME["p_0"].default

# read_case_file scope: only files under these repo roots (and these exact files)
# may be read. The FGM table returns its header only; geometry/binary files are
# refused. This keeps the tool to solver/case/doc text and out of agent internals
# (e.g. agent/.env, config.py).
_READABLE_DIR_ROOTS: tuple[str, ...] = (
    "solver/", "testcases/", "tables/", "agent/knowledge/",
)
_READABLE_EXACT_FILES: frozenset[str] = frozenset({"README.md", "agent/README.md"})
_UNREADABLE_NAMES: frozenset[str] = frozenset({"points", "faces", "owner", "neighbour"})
_UNREADABLE_SUFFIXES: frozenset[str] = frozenset(
    {".o", ".dep", ".so", ".png", ".jpg", ".pdf", ".swp"}
)
_MAX_READ_LINES: int = 800


def _json_type(py_type: type) -> str:
    return "integer" if py_type is int else "number"


def _build_case_schema() -> dict[str, Any]:
    """Derive the build_case parameter schema from the validator's spec.

    Which parameters are *required* depends on the case type, which JSON Schema
    cannot express conditionally in a way these models follow reliably. So the
    schema marks only the parameters common to every template as required and
    describes the rest; the per-template requirement is enforced in Python by
    ``validator.find_missing``, which returns a structured "still need X" result
    the model relays as a follow-up question.
    """
    properties: dict[str, Any] = {}
    # Required in the schema only where every template needs it; the rest is
    # enforced per template in code.
    always_required = set.intersection(
        *(set(t.required) for t in case_templates.TEMPLATES)
    )
    required: list[str] = []
    for spec in validator.PARAM_SPECS:
        desc = f"{spec.description} [{spec.unit}]"
        per_template = {
            t.name: t.defaults[spec.name]
            for t in case_templates.TEMPLATES if spec.name in t.defaults
        }
        if per_template:
            recommended = ", ".join(f"{v:g} for {k}" for k, v in per_template.items())
            desc += (
                f" (optional; recommended {recommended}, the value each case "
                "ships with at 300 K — say so when using it rather than "
                "presenting it as the user's choice)"
            )
        elif not spec.required and spec.default is not None:
            desc += f" (optional, default {spec.default})"
        applies_to = [
            t.name for t in case_templates.TEMPLATES if spec.name in t.parameters
        ]
        if len(applies_to) < len(case_templates.TEMPLATES):
            desc += f" (only for case_type: {', '.join(applies_to)})"
        properties[spec.name] = {"type": _json_type(spec.py_type), "description": desc}
        if spec.name in always_required:
            required.append(spec.name)

    properties["case_type"] = {
        "type": "string",
        "enum": [t.name for t in case_templates.TEMPLATES],
        "description": "which template case to build. "
        + "; ".join(
            f"{t.name}: {t.summary} — needs {', '.join(t.required)}"
            for t in case_templates.TEMPLATES
        )
        + ". REQUIRED whenever the user has not made the case obvious: if they just "
        "say 'set up a case' without naming one, ask them which case they want "
        "instead of guessing or defaulting.",
    }
    # Optional case directory name (not a physical parameter).
    properties["case_name"] = {
        "type": "string",
        "description": (
            "directory name for the generated case, e.g. "
            "'laminarSlotBurner_lowPhi' or 'laminarFlamelet_hot'. A bare "
            "descriptor is automatically prefixed with the case type, e.g. "
            "'lowPhi' -> 'laminarSlotBurner_lowPhi'. Defaults to "
            "'<case_type>_agent'."
        ),
    }
    properties["confirm_unusual"] = _CONFIRM_PROPERTY
    properties["n_steps"] = _N_STEPS_PROPERTY
    properties["use_shipped_values"] = _SHIPPED_PROPERTY
    return {"type": "object", "properties": properties, "required": required}


_SHIPPED_PROPERTY: dict[str, Any] = {
    "type": "boolean",
    "description": (
        "true ONLY when the user explicitly asks for the default / standard / "
        "shipped conditions ('use the default values'). Fills each required "
        "condition the user did not give with the value the case's own run script "
        "ships with, and reports which. Never set it to avoid asking."
    ),
}


_N_STEPS_PROPERTY: dict[str, Any] = {
    "type": "integer",
    "description": (
        "run for this many solver time steps ('run only 200 timesteps'): sets "
        "sim_time = n_steps x the time step, computed per case after the time step "
        "is known. Use it instead of sim_time, never together with it."
    ),
}


_CONFIRM_PROPERTY: dict[str, Any] = {
    "type": "boolean",
    "description": (
        "true ONLY after the user has confirmed, in a later message, values the "
        "tool reported as unusual (status confirm_values). Never set it on a "
        "first attempt: the tool ignores it for values it has not asked about."
    ),
}


def _campaign_geometry_properties() -> dict[str, Any]:
    """Geometry a campaign may hold fixed for every case (e.g. a longer outlet).

    Generated from the case plugins, like the build_case schema, so a new
    geometric parameter needs no schema edit here.
    """
    properties: dict[str, Any] = {}
    for spec in validator.PARAM_SPECS:
        templates = [t.name for t in case_templates.TEMPLATES if spec.name in t.geometric]
        if templates:
            properties[spec.name] = {
                "type": _json_type(spec.py_type),
                "description": (
                    f"held fixed for every case unless swept: {spec.description} "
                    f"[{spec.unit}] (only for case_type: {', '.join(templates)})"
                ),
            }
    return properties


TOOLS: list[dict[str, Any]] = [
    {
        "type": "function",
        "function": {
            "name": "build_case",
            "description": (
                "Build one of the FGMFoam template cases with the given "
                "parameters: compute the thermochemistry with Cantera and write "
                "the patched case under testcases/<case_name>/. ONLY call this when "
                "the user has EXPLICITLY asked to set up / build / create / "
                "configure a case AND has provided the required parameter values. "
                "Do NOT call it to explain the solver, answer a question, or when a "
                "request is vague — answer in text instead."
            ),
            "parameters": _build_case_schema(),
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_case",
            "description": (
                "Run the OpenFOAM solver pipeline for a case already built with "
                "build_case: blockMesh, setFields, decomposePar, the parallel "
                "FGMFoam solver, and reconstructPar. By default it launches the "
                "solver in the BACKGROUND and returns immediately with a log "
                "file to tail (use run_status to check progress). Requires an "
                "OpenFOAM-sourced shell with FGMFoam built. Only run when the "
                "user explicitly asks to."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "case to run (default 'laminarSlotBurner_agent')",
                    },
                    "skip_checks": {
                        "type": "boolean",
                        "description": "run even if the pre-run checks report "
                        "errors (default false). Only set this when the user has "
                        "seen the errors and explicitly asked to run anyway.",
                    },
                    "background": {
                        "type": "boolean",
                        "description": "launch in the background and return at once "
                        "(default true); set false to wait for completion",
                    },
                    "timeout_s": {
                        "type": "integer",
                        "description": "max seconds for each quick SETUP step "
                        "(blockMesh/setFields/decomposePar), as a hang guard "
                        "(default 900); the FGMFoam solver and reconstructPar are "
                        "never time-limited and run to completion",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_status",
            "description": (
                "Check the status of a (background) solver run for a case: "
                "whether it is still running or has finished, the latest "
                "simulated time, and the tail of log.FGMFoam."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "case to check (default 'laminarSlotBurner_agent')",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_cases",
            "description": (
                "Run several already-built cases as a batch. By DEFAULT they run "
                "SEQUENTIALLY (one after another). Only set max_parallel > 1 if the "
                "user explicitly asks to run several at the same time. The tool "
                "enforces the physical-core budget: it refuses any case that needs "
                "more cores than the machine has, and caps how many run at once so "
                "the concurrently-running cases never exceed the physical cores "
                "(running the rest sequentially). Launches in the background."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "the built cases to run",
                    },
                    "max_parallel": {
                        "type": "integer",
                        "description": "how many to run at once (default 1 = "
                        "sequential); set > 1 only if the user explicitly asks",
                    },
                },
                "required": ["case_names"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "kill_jobs",
            "description": (
                "Stop and discard one or more simulation jobs by case name. For a "
                "running case it kills the solver processes; for a case still "
                "queued in a batch it removes it from the queue (it will be "
                "skipped). Use when the user asks to kill, stop, cancel, discard, "
                "or remove jobs."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "the cases/jobs to stop or remove",
                    },
                },
                "required": ["case_names"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "delete_cases",
            "description": (
                "Delete (remove from disk) agent-generated cases set up earlier — "
                "e.g. after a setup mistake. Pass case_names for specific ones, or "
                "all=true to delete every generated case. It ONLY removes "
                "agent-generated cases (named <template>_<name>, e.g. "
                "laminarSlotBurner_hot or laminarFlamelet_T400) and can NEVER "
                "delete the original template cases (laminarSlotBurner, "
                "laminarFlamelet). Any running case is stopped first."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_names": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "the generated cases to delete",
                    },
                    "all": {
                        "type": "boolean",
                        "description": "delete every generated case (ignores case_names)",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "compute_thermochemistry",
            "description": (
                "Compute the FGM progress-variable bounds (Yc_u unburnt, Yc_b "
                "burnt) and the unburnt specific enthalpy ht for an equivalence "
                "ratio and inlet temperature, without building a case."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "equivalence_ratio": {
                        "type": "number",
                        "description": "fuel/air equivalence ratio phi [-]",
                    },
                    "T_inlet": {
                        "type": "number",
                        "description": "inlet (unburnt) temperature [K] (default 300)",
                    },
                    "p_0": {
                        "type": "number",
                        "description": "pressure [Pa] (default 101325)",
                    },
                },
                "required": ["equivalence_ratio"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "read_case_file",
            "description": (
                "Return the exact, verbatim contents of a solver or case file in "
                "the repository (e.g. constant/combustionProperties, 0.orig/T, "
                "system/controlDict, or a source file like "
                "solver/src/combustionModel/FGM/FGM/FGM.C). Use this whenever a "
                "question is about what a specific file contains, and answer from "
                "the returned text rather than from memory or retrieved fragments. "
                "A bare case-relative path resolves against the laminarSlotBurner "
                "template case."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "repository-relative or case-relative path to the file",
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "change_settings",
            "description": (
                "Change run-control settings of an ALREADY-BUILT case WITHOUT "
                "rebuilding it: the total simulated time (sim_time -> endTime), the "
                "fixed solver time step (sim_timestep -> deltaT), and/or the number "
                "of saved snapshots (sim_write_steps -> writeInterval = endTime / "
                "sim_write_steps). It edits only system/controlDict and does NOT "
                "change the mesh, the thermochemistry, or the initial/boundary "
                "fields (use build_case to change those). The new settings take "
                "effect on the next run_case or restart_case. Give at least one of "
                "sim_time, sim_timestep, sim_write_steps. (If you pass a condition "
                "or geometry parameter here anyway, the case is rebuilt with it.)"
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "the built case to modify (default 'laminarSlotBurner_agent')",
                    },
                    "sim_time": {
                        "type": "number",
                        "description": "new total simulated time / endTime [s]",
                    },
                    "sim_timestep": {
                        "type": "number",
                        "description": "new fixed solver time step / deltaT [s]",
                    },
                    "sim_write_steps": {
                        "type": "integer",
                        "description": "new number of write snapshots over the run "
                        "(sets writeInterval = endTime / sim_write_steps)",
                    },
                    "n_steps": _N_STEPS_PROPERTY,
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "restart_case",
            "description": (
                "Restart / continue an already-run case from its latest saved time "
                "WITHOUT re-meshing or re-initialising (no blockMesh/setFields). It "
                "reuses the existing mesh, decomposes the latest reconstructed time "
                "into processors only if needed (i.e. if the case was reconstructed "
                "or the decomposition is stale), sets startFrom latestTime, and "
                "launches FGMFoam in the background. Use this to extend or continue "
                "a run; use run_case only for the first run. If the run already "
                "reached endTime, first increase endTime with change_settings. Only "
                "call when the user explicitly asks to restart, continue, or extend."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "the case to restart (default 'laminarSlotBurner_agent')",
                    },
                    "timeout_s": {
                        "type": "integer",
                        "description": "max seconds for the synchronous decomposePar "
                        "step, if one is needed (default 900)",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_campaign",
            "description": (
                "Run a CAMPAIGN: one variable swept over several values, built, "
                "run and analysed together (e.g. flame speed against unburnt "
                "temperature, or flame length against plate temperature). Use "
                "this whenever the user wants several cases varying one "
                "quantity, not one case. "
                "Call it FIRST, WITHOUT n_procs_per_case and n_concurrent_cases, "
                "as soon as the case type and the values are clear — do not ask "
                "about cores beforehand. That call builds and starts nothing: it only "
                "reports the machine's core count "
                "and the sensible splits, which you present to the user so they "
                "choose how many cores each case gets and how many cases run at "
                "once. Never pick the split yourself."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_type": {
                        "type": "string",
                        "enum": [t.name for t in case_templates.TEMPLATES],
                        "description": "which template case to sweep",
                    },
                    "variable": {
                        "type": "string",
                        "description": (
                            "the parameter to sweep. laminarFlamelet: T_inlet. "
                            "laminarSlotBurner: T_inlet, T_plate or "
                            "inlet_velocity. Equivalence ratio and pressure "
                            "cannot be swept (fixed by the FGM table), nor can "
                            "geometry (the mesh comes from the template)."
                        ),
                    },
                    "values": {
                        "type": "array",
                        "items": {"type": "number"},
                        "description": "the values to sweep over (at least two)",
                    },
                    "n_procs_per_case": {
                        "type": "integer",
                        "description": "MPI processes per case — from the user's "
                        "choice, never invented",
                    },
                    "n_concurrent_cases": {
                        "type": "integer",
                        "description": "how many cases run at the same time — "
                        "from the user's choice, never invented",
                    },
                    "campaign_name": {
                        "type": "string",
                        "description": "short name for the campaign, used for "
                        "the case names and the output directory",
                    },
                    "T_inlet": {"type": "number", "description": "held fixed unless swept [K]"},
                    "T_plate": {"type": "number", "description": "slot burner only [K]"},
                    "inlet_velocity": {"type": "number", "description": "slot burner only [m/s]"},
                    "sim_time": {"type": "number", "description": "per case [s]"},
                    "sim_timestep": {"type": "number", "description": "per case [s]"},
                    "sim_write_steps": {"type": "integer", "description": "snapshots per case"},
                    **_campaign_geometry_properties(),
                    "flow_throughs": {
                        "type": "number",
                        "description": "run each case for this many flow-through "
                        "times (domain length / inlet velocity), so sim_time is "
                        "set per case. Slot burner only. Prefer it over sim_time "
                        "when inlet velocity or plate thickness is swept.",
                    },
                    "confirm_unusual": _CONFIRM_PROPERTY,
                    "n_steps": _N_STEPS_PROPERTY,
                    "use_shipped_values": _SHIPPED_PROPERTY,
                    "launch": {
                        "type": "boolean",
                        "description": "false = build and save the campaign but do "
                        "NOT start it (the user starts it later with "
                        "start_campaign). Default true. Use false whenever the "
                        "user says to set up, prepare, or wait before running.",
                    },
                },
                "required": ["case_type", "variable", "values"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "start_campaign",
            "description": (
                "Start a campaign that was prepared earlier with "
                "run_campaign(launch=false), using the CPU plan saved with it. "
                "Only when the user explicitly asks to start or run it."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "campaign_name": {"type": "string", "description": "which campaign"},
                },
                "required": ["campaign_name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "analyse_campaign",
            "description": (
                "Collect a finished (or partly finished) campaign into a results "
                "table and a graph: flame speed against the swept variable for "
                "laminarFlamelet, flame length for laminarSlotBurner. Cases that "
                "flashed back or blew off are kept in the table as outcomes and "
                "marked on the graph. Writes results.csv and a .png next to the "
                "campaign. Use it when the user asks for the results, the trend "
                "or a plot."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "campaign_name": {"type": "string", "description": "which campaign"},
                    "temperature_threshold": {
                        "type": "number",
                        "description": "flame-front temperature [K] (default 900)",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "analyse_run",
            "description": (
                "Interpret a finished (or running) case: how the solver behaved "
                "— the Courant number actually reached, simulated time, "
                "execution time, whether it diverged — and what the flame did. "
                "For laminarFlamelet that is the LAMINAR FLAME SPEED, measured "
                "from the tracked front with the step-function initialisation "
                "transient excluded. For laminarSlotBurner it is the FLAME "
                "LENGTH on the symmetry axis, or a report of FLASHBACK (flame "
                "upstream of the flame holder) or BLOW-OFF (no flame left in "
                "the domain). Use it when the user asks how a run went, what "
                "the flame speed or flame length was, or whether the result "
                "looks trustworthy. Read-only; it samples the written fields "
                "with OpenFOAM's postProcess and changes nothing."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "the case to analyse",
                    },
                    "temperature_threshold": {
                        "type": "number",
                        "description": "temperature [K] marking the flame front "
                        "(default 900, just above methane's ~810 K autoignition "
                        "temperature). Only change it if the user asks.",
                    },
                    "log_only": {
                        "type": "boolean",
                        "description": "report only the solver log and skip "
                        "field sampling (default false). Use when OpenFOAM is "
                        "not available in the shell.",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "check_case",
            "description": (
                "Check a built case before running it: whether its pressure, "
                "initial enthalpy and progress variable lie inside the FGM "
                "table's tabulated range, whether the time step is sane against "
                "the mesh (estimated Courant number), whether the decomposition "
                "and write interval make sense. Fast and read-only — it changes "
                "nothing and starts nothing. Use it when the user asks whether a "
                "case looks right, or before a long run. run_case performs the "
                "same checks automatically and refuses to launch on errors."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "case_name": {
                        "type": "string",
                        "description": "the built case to check",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "switch_model",
            "description": (
                "Switch the model answering this conversation, or report which "
                "models are available. Call it WITH 'model' only when the user "
                "asks to use, switch to, or try a specific model or provider "
                "('use Claude', 'switch to the local model'). Call it WITHOUT "
                "'model' when they are only asking what exists or whether "
                "something is available ('which models can I use?', 'is Claude "
                "supported?') — that lists the options and switches nothing. "
                "A switch applies from the next reply onwards, lasts for this "
                "session only, and keeps the conversation; it changes nothing "
                "about the solver or any case, and is undone by switching back."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "model": {
                        "type": "string",
                        "description": (
                            "which model to switch to: a preset name "
                            f"({', '.join(p.name for p in providers.PRESETS)}), a "
                            "'preset:model-id' pair such as "
                            "'claude:claude-haiku-4-5', or a bare model id to keep "
                            "the current provider. Omit to just list the available "
                            "models without switching."
                        ),
                    },
                },
                "required": [],
            },
        },
    },
]

# Tools the stateless tool layer cannot execute because they mutate the live
# session (the provider client). ChatAgent intercepts these; the entry below
# exists so a direct call still fails loudly instead of silently doing nothing.
SESSION_TOOLS: frozenset[str] = frozenset({"switch_model"})

# Which model is currently answering. Set by ChatAgent (including after a
# /model switch) purely so a generated case can record what routed the request;
# nothing in the deterministic pipeline reads it.
_SESSION: dict[str, Any] = {"provider": None, "chat_model": None}


def set_session_model(provider: Optional[str], chat_model: Optional[str]) -> None:
    """Record the active provider/model for provenance."""
    _SESSION["provider"] = provider
    _SESSION["chat_model"] = chat_model


# The case built most recently in this conversation. Tools that act on one case
# default to it rather than to the fixed "<template>_agent" name: a small model
# often calls check_case/run_case/change_settings with no case_name right after a
# build, and the fixed default then reports "not found" for a case that exists.
# Also the base for a rebuild (see _inherit_from_last_build). ChatAgent resets it
# per conversation.
_LAST_BUILD: dict[str, Any] = {}

# Tools that act on a single existing case and take an optional case_name.
_SINGLE_CASE_TOOLS: frozenset[str] = frozenset({
    "run_case", "run_status", "change_settings", "restart_case", "check_case",
    "analyse_run",
})


# A campaign request that stopped at "how should the cores be split?". When the
# user answers, the chat loop completes this exact call (pending_campaign()).
_PENDING_CAMPAIGN: dict[str, Any] = {}


def reset_session_cases() -> None:
    """Forget the last built case — a new conversation starts clean."""
    _LAST_BUILD.clear()
    _PENDING_CAMPAIGN.clear()


def pending_campaign() -> dict[str, Any]:
    """The arguments of a campaign waiting for its CPU plan ({} if none)."""
    return dict(_PENDING_CAMPAIGN)


def _default_case_name(args: dict[str, Any]) -> str:
    """The case a single-case tool acts on: the one named, else the last built."""
    return (args.get("case_name") or _LAST_BUILD.get("case_name")
            or of_writer.DEFAULT_CASE_NAME)


def _no_case_named_error() -> Optional[dict[str, Any]]:
    """Error for a single-case tool called with no name and nothing to default to."""
    if _LAST_BUILD or of_writer.case_path(of_writer.DEFAULT_CASE_NAME).is_dir():
        return None
    return {
        "status": "error",
        "error": ("no case_name given and no case has been built in this "
                  "conversation. Ask the user which case they mean."),
        "existing_cases": of_writer.generated_cases(),
    }


def missing_case_type_result() -> dict[str, Any]:
    """The reply that makes the model ask which case to set up."""
    return {
        "status": "missing_case_type",
        "missing": ["case_type"],
        "choices": [
            {
                "case_type": t.name,
                "description": t.summary,
                "needs": list(t.required),
            }
            for t in case_templates.TEMPLATES
        ],
        "message": (
            "Which case should be set up? Ask the user to choose — do "
            "NOT pick one for them. Options: "
            + "; ".join(
                f"{t.name} ({t.summary})" for t in case_templates.TEMPLATES
            )
        ),
    }


# change_settings edits these in place; any other parameter means a rebuild.
_RUN_CONTROLS: frozenset[str] = frozenset({"sim_time", "sim_timestep", "sim_write_steps"})


# Lengths are in metres. Every geometric dimension of an FGMFoam case is far below
# this, so a larger value is a unit mistake (a small model passes "5" for 5 mm)
# that would otherwise mesh a burner a thousand times too big.
MAX_LENGTH_M: float = 0.1


# Temperatures are in kelvin. Anything outside this band is a wrong unit (degC)
# or a value meant for another parameter (an equivalence ratio of 0.7 is not a
# temperature) — never a premixed methane-air condition.
TEMPERATURE_RANGE_K: tuple[float, float] = (200.0, 3000.0)


def _unit_errors(params: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Reject lengths clearly in mm (not metres) and temperatures not in kelvin."""
    lengths, temperatures = {}, {}
    low, high = TEMPERATURE_RANGE_K
    for spec in validator.PARAM_SPECS:
        value = params.get(spec.name)
        if value is None or spec.unit not in ("m", "K"):
            continue
        values = value if isinstance(value, list) else [value]
        try:
            numbers = [float(v) for v in values]
        except (TypeError, ValueError):
            continue
        if spec.unit == "m" and any(abs(n) > MAX_LENGTH_M for n in numbers):
            lengths[spec.name] = value
        if spec.unit == "K" and any(not low <= n <= high for n in numbers):
            temperatures[spec.name] = value
    if not lengths and not temperatures:
        return None
    messages = []
    if lengths:
        messages.append(
            "Lengths are in METRES. "
            + ", ".join(f"{k}={v}" for k, v in lengths.items())
            + f" is larger than {MAX_LENGTH_M} m, which is almost certainly "
            "millimetres. Convert (e.g. 5 mm = 0.005) and call the tool again."
        )
    if temperatures:
        messages.append(
            "Temperatures are in KELVIN. "
            + ", ".join(f"{k}={v}" for k, v in temperatures.items())
            + f" is outside {low:g}-{high:g} K: a wrong unit, or a value that "
            "belongs to another parameter. Do not retry with a guess — ask the user."
        )
    return {
        "status": "error",
        "error": "value out of physical range",
        "parameters": {**lengths, **temperatures},
        "message": " ".join(messages),
    }


def _chosen_parameters(case_dir: Any) -> tuple[Optional[str], dict[str, Any]]:
    """(case type, parameters the user chose) from an existing case's provenance."""
    record = provenance.read(case_dir) if case_dir.is_dir() else None
    if not record:
        return None, {}
    sources = record.get("parameter_sources", {})
    chosen = {k: v for k, v in record.get("parameters", {}).items()
              if sources.get(k) in ("user", "inherited")}
    return record.get("case", {}).get("type"), chosen


def _rebuild_base(args: dict[str, Any], template: Any) -> tuple[Optional[str], dict[str, Any]]:
    """The case a rebuild starts from: the one named, if it exists, else the last
    built. Returns (its name, the parameters the user chose for it)."""
    if args.get("case_name"):
        try:
            named = of_writer.case_path(str(args["case_name"]), template)
        except ValueError:
            named = None
        if named is not None and named.is_dir():
            case_type, chosen = _chosen_parameters(named)
            return (named.name, chosen) if case_type == template.name else (None, {})
    if _LAST_BUILD.get("case_type") == template.name:
        return _LAST_BUILD["case_name"], dict(_LAST_BUILD.get("chosen", {}))
    return None, {}


def _inherit(params: dict[str, Any], template: Any,
             args: dict[str, Any]) -> tuple[Optional[str], dict[str, Any]]:
    """Parameters the user chose for the case being rebuilt and did not change
    now. A rebuild ("now with a 5 mm plate") restates only what changes; without
    this a small model drops the other conditions, or invents them. Only
    user-chosen values carry over — never defaults or the estimated step."""
    base, chosen = _rebuild_base(args, template)
    return base, {k: v for k, v in chosen.items() if k not in params}


def _table_constraint_error(
    args: dict[str, Any], key: str, table_value: float, label: str
) -> Optional[dict[str, Any]]:
    """Return an error dict if ``args[key]`` differs from the table's fixed value.

    The FGM table bakes in a single equivalence ratio and pressure, so building a
    case at a different value would be inconsistent with the tabulated chemistry.
    """
    requested = args.get(key)
    if requested is None:
        return None
    try:
        if abs(float(requested) - table_value) <= 1e-9:
            return None
    except (TypeError, ValueError):
        return None
    return {
        "status": "table_unsupported",
        "parameter": key,
        "requested": requested,
        "table_value": table_value,
        "message": (
            f"The FGM table is generated at a fixed {label} of {table_value}, "
            f"which is baked into the table (not a free dimension). "
            f"{key}={requested} is not supported by the current table — regenerate "
            f"the FGM table at that value, or use {table_value}."
        ),
    }


# Unusual values (outside a parameter's usual range) are confirmed by the user,
# not by the model: confirm_unusual counts only for values the tool already asked
# about in an EARLIER turn, so a model cannot flag and confirm in one go.
_TURN: list[int] = [0]
_ASKED: dict[tuple[str, float], int] = {}


def begin_turn() -> None:
    """Called by the chat loop at the start of every user message."""
    _TURN[0] += 1


def _unconfirmed(unusual: list[dict[str, Any]], confirmed: bool) -> list[dict[str, Any]]:
    """The unusual values the user has not confirmed yet; marks them as asked."""
    pending = []
    for entry in unusual:
        key = (entry["parameter"], round(entry["value"], 12))
        asked = _ASKED.get(key)
        if confirmed and asked is not None and asked < _TURN[0]:
            continue
        _ASKED.setdefault(key, _TURN[0])
        pending.append(entry)
    return pending


def _confirm_result(pending: list[dict[str, Any]], template: Any,
                    missing: Optional[list[str]] = None) -> dict[str, Any]:
    lines = []
    for e in pending:
        low, high = e["usual_range"]
        line = (f"{e['parameter']} = {e['value']:g} {e['unit']} (usual for "
                f"{template.name}: {low:g} to {high:g} {e['unit']})")
        if e.get("maybe_meant"):
            line += f"; did the user mean {e['maybe_meant']}?"
        lines.append(line)
    message = ("Nothing was built. These values are unusual: " + "; ".join(lines)
               + ". Ask the user whether they are sure (and whether the unit is "
               "right). Only if they confirm in their next message, call again with "
               "confirm_unusual=true.")
    if missing:
        message += " Also still needed: " + ", ".join(missing) + "."
    return {"status": "confirm_values", "case_type": template.name,
            "unusual": pending, "message": message}


def _shipped_fill(params: dict[str, Any], template: Any, args: dict[str, Any],
                  swept: str = "") -> dict[str, Any]:
    """The shipped values for required conditions the user did not give, if asked."""
    if not args.get("use_shipped_values"):
        return {}
    return {k: v for k, v in template.shipped.items()
            if k in template.required and params.get(k) is None and k != swept}


def _n_steps(args: dict[str, Any]) -> tuple[Optional[int], Optional[dict[str, Any]]]:
    """(n_steps, error): the requested number of time steps, validated."""
    value = args.get("n_steps")
    if value in (None, ""):
        return None, None
    try:
        n = int(float(value))
    except (TypeError, ValueError):
        return None, {"status": "error", "error": f"n_steps must be a whole number, got {value!r}"}
    if n <= 0 or n != float(value):
        return None, {"status": "error", "error": f"n_steps must be a positive whole number, got {value!r}"}
    if args.get("sim_time") is not None:
        return None, {"status": "error", "error": (
            "give n_steps or sim_time, not both: n_steps sets sim_time = n_steps x "
            "the time step")}
    return n, None


def _steps_to_controls(n_steps: int, timestep: float, write_steps: int) -> tuple[float, int]:
    """(sim_time, write snapshots) for ``n_steps`` steps of ``timestep``.

    The snapshot count is lowered to a divisor of ``n_steps``, so every write
    falls on a time step (the write interval is a whole number of steps).
    """
    write_steps = max(1, min(int(write_steps), n_steps))
    while n_steps % write_steps:
        write_steps -= 1
    return float(f"{n_steps * timestep:.12g}"), write_steps


def unburnt_temperature_outside_table(t_inlet: float, phi: float = 0.7,
                                      p_0: Optional[float] = None) -> Optional[str]:
    """Why an unburnt mixture at ``t_inlet`` K is outside the shipped table, or None."""
    problems = _table_temperature_problems(
        {"T_inlet": t_inlet, "equivalence_ratio": phi, "p_0": p_0}, None)
    return " ".join(problems) or None


def plate_temperature_outside_table(t_plate: float) -> Optional[str]:
    """Why a plate at ``t_plate`` K is outside the shipped table, or None."""
    return " ".join(_table_temperature_problems({"T_plate": t_plate}, None)) or None


def _table_temperature_problems(params: dict[str, Any], template: Any) -> list[str]:
    """Temperatures the FGM table cannot represent, each explained (hard limits).

    ``T_inlet``: the unburnt enthalpy must lie in the table's enthalpy range, and
    the temperature not below the coldest tabulated one. ``T_plate``: a wall at a
    fixed temperature is inverted to an enthalpy by ``lookUpEnthalpy``, which can
    only find temperatures the table returns. Lists (campaign values) are checked
    value by value.
    """
    info = fgm_table.read_table_safe(of_writer.TABLE_SRC)
    if info is None or info.dimension < 2:
        return []
    defaults = template.defaults if template is not None else {}
    phi = float(params.get("equivalence_ratio") or defaults.get("equivalence_ratio", 0.7))
    p_0 = float(params.get("p_0") or defaults.get("p_0", info.pressure))
    h_lo, h_hi = info.range_of(1)
    t_lo, t_hi = info.temperature_range or (0.0, float("inf"))

    def numbers(value: Any) -> list[float]:
        values = value if isinstance(value, list) else [value]
        return [float(v) for v in values if isinstance(v, (int, float))]

    problems = []
    for t_inlet in numbers(params.get("T_inlet")):
        try:
            ht = cantera_tool.unburnt_enthalpy(phi, t_inlet, p_0)
        except Exception:  # noqa: BLE001 - leave it to the later checks
            continue
        if h_lo <= ht <= h_hi and t_inlet >= t_lo:
            continue
        t_max = cantera_tool.max_unburnt_temperature(phi, p_0, h_hi)
        problems.append(
            f"An unburnt-gas temperature of {t_inlet:g} K is outside the FGM table: "
            f"at phi = {phi:g} and {p_0:g} Pa the table covers unburnt mixtures from "
            f"{t_lo:.0f} K to about {t_max:.0f} K (its enthalpy range is "
            f"[{h_lo:.4g}, {h_hi:.4g}] J/kg; {t_inlet:g} K gives {ht:.4g} J/kg).")
    for t_plate in numbers(params.get("T_plate")):
        if not t_lo <= t_plate <= t_hi:
            problems.append(
                f"A plate temperature of {t_plate:g} K is outside the FGM table, whose "
                f"temperatures run from {t_lo:.0f} to {t_hi:.0f} K: the plate's "
                "enthalpy boundary condition can only find temperatures the table holds.")
    return problems


def _outside_table(params: dict[str, Any], template: Any) -> Optional[dict[str, Any]]:
    """Refuse temperatures the FGM table does not cover, at build time.

    ``precheck`` finds the same problem, but only when a case is checked or run;
    by then the user has been asked for every other value and a case has been
    written for conditions that cannot be simulated.
    """
    problems = _table_temperature_problems(params, template)
    if not problems:
        return None
    info = fgm_table.read_table_safe(of_writer.TABLE_SRC)
    defaults = template.defaults if template is not None else {}
    phi = float(params.get("equivalence_ratio") or defaults.get("equivalence_ratio", 0.7))
    p_0 = float(params.get("p_0") or defaults.get("p_0", info.pressure))
    return {
        "status": "outside_table",
        "case_type": template.name if template is not None else None,
        "max_T_inlet": round(cantera_tool.max_unburnt_temperature(phi, p_0, info.range_of(1)[1])),
        "table_temperature_range_K": [round(t) for t in info.temperature_range or ()],
        "message": " ".join(problems) + (
            " Nothing was built and no other values are needed: tell the user, and "
            "ask for a temperature inside the table (or explain that this needs a "
            "regenerated table)."),
    }


def _run_build_case(args: dict[str, Any], session: bool = True) -> dict[str, Any]:
    """Build one case. ``session=False`` (a campaign's own builds) neither
    inherits from nor becomes the conversation's last built case: a campaign
    spells out every parameter, and its cases are not what "that case" means."""
    # case_type and case_name are not physical parameters: resolve them before
    # validating, because the case type decides which parameters are required.
    requested_type = args.get("case_type")
    if requested_type:
        try:
            template = case_templates.find_template(requested_type)
        except ValueError as exc:
            return {"status": "error", "error": str(exc)}
    else:
        # No case type given. Infer it when the supplied parameters can only
        # mean one case (an inlet velocity implies the slot burner); otherwise
        # ASK rather than silently defaulting — quietly building the wrong
        # geometry is worse than one extra question, and enforcing it here
        # means it holds for every model.
        template = case_templates.infer_template(
            validator.canonical_name(name) for name in args
        )
        if template is None and args.get("case_name"):
            # "laminarFlamelet_cond" names its template in the prefix.
            template = case_templates.template_for_case_name(str(args["case_name"]))
        if template is None:
            return missing_case_type_result()
    case_name = args.get("case_name") or template.default_case_name

    # The FGM table bakes in a fixed equivalence ratio and pressure; reject other
    # values rather than build a case inconsistent with the tabulated chemistry.
    for key, table_value, label in (
        ("equivalence_ratio", TABLE_EQUIVALENCE_RATIO, "equivalence ratio"),
        ("p_0", TABLE_PRESSURE, "pressure (Pa)"),
    ):
        error = _table_constraint_error(args, key, table_value, label)
        if error is not None:
            return error

    params = validator.coerce_types(args)
    # coerce_types drops what it cannot read, so that a missing REQUIRED value is
    # asked about again — but an optional one (any geometry) would silently fall
    # back to its default. Refuse instead: the 3B passed plate_thickness=[2, 3, 4]
    # and got a "built" case with the default plate.
    unreadable = sorted(
        validator.canonical_name(k) for k, v in args.items()
        if v is not None and validator.canonical_name(k) in validator.SPEC_BY_NAME
        and validator.canonical_name(k) not in params
    )
    if unreadable:
        return {
            "status": "error",
            "error": "unreadable parameter values",
            "parameters": {k: args.get(k) for k in unreadable},
            "message": (
                "Each parameter takes ONE number: " + ", ".join(unreadable)
                + ". To try several values, use run_campaign with variable and values."
            ),
        }
    unit_error = _unit_errors(params)
    if unit_error is not None:
        return unit_error
    base, inherited = _inherit(params, template, args) if session else (None, {})
    params = {**inherited, **params}
    shipped = _shipped_fill(params, template, args)
    params = {**shipped, **params}
    # Checked before asking for anything else: a 3000 K unburnt mixture is not
    # worth collecting a velocity and a plate temperature for.
    outside = _outside_table(params, template)
    if outside:
        return outside
    missing = validator.find_missing(params, template)
    # Inherited values were confirmed (or usual) when they were first built.
    # A campaign asks once for all its cases (session=False builds are its own).
    pending = _unconfirmed(
        validator.unusual_values({k: v for k, v in params.items() if k not in inherited},
                                 template),
        bool(args.get("confirm_unusual"))) if session else []
    if pending:
        return _confirm_result(pending, template, missing)
    if missing:
        return {
            "status": "missing_parameters",
            "case_type": template.name,
            "missing": missing,
            "message": (
                f"Cannot build the {template.name} case yet; still need: "
                + ", ".join(missing)
            ),
        }
    n_steps, step_error = _n_steps({**args, "sim_time": args.get("sim_time")})
    if step_error:
        return step_error
    supplied = set(params) - set(inherited) - set(shipped)
    params = validator.apply_defaults(params, template)
    problem = template.check(params) if template.check else None
    if problem:
        return {"status": "error", "error": "impossible geometry", "message":
                problem + ". Nothing was built; ask the user for consistent values."}
    thermo = cantera_tool.compute_thermochemistry(
        params["equivalence_ratio"], params["T_inlet"], params["p_0"]
    )

    # The solver cannot adapt its step, so choose one from the mesh and the flow
    # rather than falling back to the value this case happens to ship with.
    estimated: set[str] = set()
    dt_estimate = None
    if "sim_timestep" not in supplied:
        # With n_steps the end time follows from the step, so the step is not
        # snapped to the default write interval (0.4 m/s gave 3.003e-6 s).
        dt_estimate = timestep.estimate(
            template.name, {**params, "sim_time": None} if n_steps else params, thermo)
        if dt_estimate is not None:
            params["sim_timestep"] = float(f"{dt_estimate.timestep:.12g}")
            estimated.add("sim_timestep")
    if n_steps is not None:
        # After the time step is known (estimated or given), per case.
        params["sim_time"], params["sim_write_steps"] = _steps_to_controls(
            n_steps, float(params["sim_timestep"]), int(params["sim_write_steps"]))
        estimated |= {"sim_time"}
    case_dir = of_writer.prepare_case(case_name, template)
    replacements = of_writer.build_replacements(params, thermo, template)
    of_writer.patch_placeholders(case_dir, replacements, template)
    result = {
        "status": "built",
        "case_type": template.name,
        "case_name": case_dir.name,
        "case_dir": str(case_dir),
        "parameters": params,
        "thermochemistry": {"Yc_u": thermo.Yc_u, "Yc_b": thermo.Yc_b, "ht": thermo.ht},
        "write_interval": params["sim_time"] / params["sim_write_steps"],
    }
    if shipped:
        result["shipped_values_used"] = {
            "values": shipped,
            "note": ("the user asked for the default conditions; these are the values "
                     f"{template.name}'s own run script ships with. Tell the user."),
        }
    if n_steps is not None:
        result["n_steps"] = {
            "n_steps": n_steps, "sim_time": params["sim_time"],
            "sim_timestep": params["sim_timestep"],
            "sim_write_steps": params["sim_write_steps"],
            "note": (f"sim_time = {n_steps} steps x {params['sim_timestep']:g} s = "
                     f"{params['sim_time']:g} s"),
        }
    # Say which run-control values were not chosen by the user, so the model
    # reports them as recommendations rather than as settings they asked for.
    defaulted = sorted(set(template.defaults) - supplied - estimated)
    if defaulted:
        result["used_recommended_values"] = {
            name: template.defaults[name] for name in defaulted
        }
    if dt_estimate is not None:
        result["sim_timestep_estimated"] = {
            "value": dt_estimate.timestep,
            "target_courant": dt_estimate.cfl,
            "cell_size": dt_estimate.cell_size,
            "velocity": dt_estimate.velocity,
            "explanation": dt_estimate.explanation,
        }
    # A recommended step chosen at 300 K may be too large for a hotter mixture,
    # and the solver runs at a fixed step, so it will not adapt on its own.
    warnings: list[dict[str, Any]] = []
    # Only warn when the shipped value was used unchanged: an estimated step
    # already accounts for the temperature through the density ratio and the
    # flame speed, so repeating the rule of thumb would be noise.
    if "sim_timestep" not in supplied and dt_estimate is None:
        advice = case_templates.timestep_advice(template, params)
        if advice is not None:
            result["timestep_warning"] = advice
            warnings.append(advice)

    # Make the case self-describing: what was chosen, what was filled in, which
    # table and mechanism produced the numbers. Never allowed to fail a build.
    try:
        record = provenance.build_record(
            case_dir=case_dir,
            template=template,
            params=params,
            supplied=supplied,
            estimated=estimated,
            inherited=set(inherited),
            thermo=thermo,
            table_path=of_writer.TABLE_SRC,
            repo_root=REPO_ROOT,
            version=config.AGENT_VERSION,
            model=dict(_SESSION),
            warnings=warnings,
            shipped=set(shipped),
        )
        result["provenance"] = str(provenance.write(case_dir, record).name)
    except Exception as exc:  # noqa: BLE001 - a record is not worth losing a case
        result["provenance_error"] = f"{type(exc).__name__}: {exc}"
    if inherited:
        result["inherited_from"] = {
            "case_name": base,
            "parameters": inherited,
            "note": "kept from the previous case because the user did not change them",
        }
    if not session:
        return result
    _LAST_BUILD.clear()
    _LAST_BUILD.update({
        "case_name": case_dir.name,
        "case_type": template.name,
        "chosen": {k: params[k] for k in (supplied | set(inherited))},
    })
    return result


def _run_run_campaign(args: dict[str, Any]) -> dict[str, Any]:
    case_type = args.get("case_type")
    if not case_type:
        # find_template(None) means "the default case"; a campaign must never
        # be built on a case type nobody chose.
        return missing_case_type_result()
    try:
        template = case_templates.find_template(case_type)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}

    # The variable is never guessed here: "equivalence ratios 0.6, 0.7, 0.8" with
    # the variable left out once became a T_inlet sweep at 0.6-0.8 K. The chat
    # loop fills it in only from the user's own words (chat._variable_from_user).
    variable = str(args.get("variable") or "")
    problem = campaign_module.check_variable(template.name, variable)
    if problem:
        return {"status": "error", "error": problem}

    values = args.get("values") or []
    try:
        values = [float(v) for v in values]
    except (TypeError, ValueError):
        return {"status": "error", "error": "values must be numbers"}
    if len(values) < 2:
        return {
            "status": "error",
            "error": "a campaign needs at least two values to sweep over",
        }
    swept = {**validator.coerce_types(args), validator.canonical_name(variable): values}
    unit_error = _unit_errors(swept)
    if unit_error is not None:
        return unit_error
    outside = _outside_table(swept, template)
    if outside is not None:
        return outside
    pending = _unconfirmed(validator.unusual_values(swept, template),
                           bool(args.get("confirm_unusual")))
    if pending:
        return _confirm_result(pending, template)
    if template.check:
        fixed = validator.apply_defaults(
            {k: v for k, v in validator.coerce_types(args).items()
             if k in template.parameters and k != validator.canonical_name(variable)},
            template)
        for v in values:
            problem = template.check({**fixed, validator.canonical_name(variable): v})
            if problem:
                return {"status": "error", "error": "impossible geometry",
                        "message": f"{variable} = {v:g}: {problem}. Nothing was built."}

    n_steps, step_error = _n_steps(args)
    if step_error:
        return step_error
    if n_steps is not None and args.get("flow_throughs"):
        return {"status": "error", "error": (
            "give n_steps or flow_throughs, not both: each sets the simulated time")}
    base: dict[str, Any] = {
        name: args[name] for name in
        ("T_inlet", "T_plate", "inlet_velocity", "sim_time", "sim_timestep",
         "sim_write_steps", "n_steps") + tuple(template.geometric)
        if args.get(name) is not None
    }
    base.pop(variable, None)
    shipped = _shipped_fill(base, template, args, swept=validator.canonical_name(variable))
    base.update(shipped)
    # Checked once, before any case is built, rather than case by case.
    missing = [k for k in validator.find_missing({**base, variable: values[0]}, template)]
    if missing:
        offer = {k: v for k, v in template.shipped.items() if k in missing}
        return {
            "status": "missing_parameters", "case_type": template.name,
            "missing": missing,
            "message": (
                "Nothing was built; the campaign still needs " + ", ".join(missing)
                + ". These conditions have no defaults: ask the user for them."
                + (" If the user wants the conditions the case ships with, that is "
                   + ", ".join(f"{k} = {v:g}" for k, v in offer.items())
                   + " (call again with use_shipped_values=true once they say so)."
                   if offer else "")),
        }

    cores = runner.physical_core_count(
        config.load_config(require_key=False).physical_cores
    )
    procs = args.get("n_procs_per_case")
    concurrent = args.get("n_concurrent_cases")

    # The CPU plan is the user's decision: it changes wall-clock time by a large
    # factor and depends on what else the machine is doing. Ask, never assume.
    if not procs or not concurrent:
        _PENDING_CAMPAIGN.clear()
        _PENDING_CAMPAIGN.update({k: v for k, v in args.items()
                                  if k not in ("n_procs_per_case", "n_concurrent_cases")})
        _PENDING_CAMPAIGN["variable"] = variable
        return {
            "status": "missing_cpu_plan",
            "missing": ["n_procs_per_case", "n_concurrent_cases"],
            "physical_cores": cores,
            "n_cases": len(values),
            "options": campaign_module.cpu_options(len(values), cores),
            "message": (
                f"This campaign is {len(values)} cases and the machine has "
                f"{cores} physical cores. ASK THE USER how to spend them — how "
                "many cores per case and how many cases at once — and present "
                "the options below. Do NOT choose for them."
            ),
        }

    procs, concurrent = int(procs), int(concurrent)
    problem = campaign_module.validate_cpu_plan(procs, concurrent, cores)
    if problem:
        return {
            "status": "error", "error": problem, "physical_cores": cores,
            "options": campaign_module.cpu_options(len(values), cores),
        }


    name = campaign_module.safe_name(
        args.get("campaign_name") or f"{variable}sweep"
    )
    flow_throughs = args.get("flow_throughs")
    try:
        plan = campaign_module.plan(
            name, template.name, variable, values, base, procs, concurrent,
            flow_throughs=float(flow_throughs) if flow_throughs else None,
        )
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}

    built, failed = [], []
    for entry in plan.cases:
        # Shipped values go in as shipped (not as the user's), so each case's
        # provenance says where they came from.
        result = _run_build_case({
            "case_type": template.name,
            "case_name": entry["case_name"],
            **{k: v for k, v in entry["params"].items() if k not in shipped},
            **({"use_shipped_values": True} if shipped else {}),
        }, session=False)
        if result.get("status") == "built":
            built.append(result["case_name"])
        else:
            failed.append({"case": entry["case_name"], "result": result})
    if failed:
        return {
            "status": "error",
            "error": "some cases could not be built; nothing was run",
            "failures": failed,
        }
    summary = {
        "campaign": plan.name, "case_type": template.name, "variable": variable,
        "values": values, "n_cases": len(built),
        **({"n_steps": n_steps} if n_steps is not None else {}),
        **({"shipped_values_used": shipped} if shipped else {}),
        "cpu_plan": {
            "n_procs_per_case": procs, "n_concurrent_cases": concurrent,
            "physical_cores": cores, "cores_used": procs * concurrent,
        },
        "cases": [
            {"case": e["case_name"], variable: e["value"],
             "sim_time": e["params"].get("sim_time")}
            for e in plan.cases
        ],
    }
    if args.get("launch") is False:
        _PENDING_CAMPAIGN.clear()
        plan.status = "prepared"
        campaign_module.save(plan)
        return {
            **summary, "status": "prepared",
            "next": f"Built and saved, not started. Start it with "
                    f"start_campaign('{plan.name}') when the user says so.",
        }

    _PENDING_CAMPAIGN.clear()
    plan.status = "launched"
    campaign_module.save(plan)
    launched = _run_run_cases({
        "case_names": built, "max_parallel": concurrent,
    })
    return {
        "status": launched.get("status", "unknown"),
        "campaign": plan.name,
        "case_type": template.name,
        "variable": variable,
        "values": values,
        "n_cases": len(built),
        "cpu_plan": {
            "n_procs_per_case": procs,
            "n_concurrent_cases": concurrent,
            "physical_cores": cores,
            "cores_used": procs * concurrent,
        },
        "cases": built,
        "batch": launched,
        "next": (
            "Runs are in the background. Check progress with run_status on any "
            f"case, then call analyse_campaign('{plan.name}') for the table and "
            "the graph."
        ),
    }


def _run_start_campaign(args: dict[str, Any]) -> dict[str, Any]:
    plan = campaign_module.load(str(args.get("campaign_name") or ""))
    if plan is None:
        return {"status": "error", "error": "no such campaign; available: "
                + (", ".join(campaign_module.list_campaigns()) or "none")}
    if plan.status == "launched":
        return {"status": "error", "error": f"campaign {plan.name!r} was already "
                "started; use run_status or analyse_campaign"}
    cores = runner.physical_core_count(
        config.load_config(require_key=False).physical_cores
    )
    # The machine may have changed since the plan was saved.
    problem = campaign_module.validate_cpu_plan(
        plan.n_procs_per_case, plan.n_concurrent_cases, cores
    )
    if problem:
        return {"status": "error", "error": problem}
    names = [e["case_name"] for e in plan.cases]
    missing = [n for n in names if not of_writer.case_path(n).is_dir()]
    if missing:
        return {"status": "error", "error": "cases no longer built: " + ", ".join(missing)}
    plan.status = "launched"
    campaign_module.save(plan)
    launched = _run_run_cases({"case_names": names, "max_parallel": plan.n_concurrent_cases})
    return {"status": launched.get("status", "unknown"), "campaign": plan.name,
            "n_cases": len(names), "batch": launched}


def _run_analyse_campaign(args: dict[str, Any]) -> dict[str, Any]:
    name = args.get("campaign_name")
    if not name:
        available = campaign_module.list_campaigns()
        if len(available) == 1:
            name = available[0]
        else:
            return {
                "status": "error",
                "error": "which campaign? " + (
                    "none have been run yet" if not available
                    else "available: " + ", ".join(available)
                ),
            }
    plan = campaign_module.load(str(name))
    if plan is None:
        return {
            "status": "error",
            "error": f"no campaign named {name!r}; "
                     f"available: {', '.join(campaign_module.list_campaigns()) or 'none'}",
        }
    threshold = args.get("temperature_threshold")
    results = campaign_module.collect(
        plan, threshold=float(threshold) if threshold else postrun.FLAME_TEMPERATURE
    )
    csv_path = campaign_module.write_csv(plan, results)
    png_path, plot_error = campaign_module.plot(plan, results)
    results["csv"] = str(csv_path.relative_to(REPO_ROOT))
    if png_path:
        results["plot"] = str(png_path.relative_to(REPO_ROOT))
    if plot_error:
        results["plot_error"] = plot_error
    results["status"] = "collected"
    return results


def _run_analyse_run(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    case_dir = of_writer.case_path(case_name)
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case '{case_dir.name}' not found"}
    threshold = args.get("temperature_threshold")
    return postrun.analyse(
        case_dir,
        threshold=float(threshold) if threshold else postrun.FLAME_TEMPERATURE,
        sample=not bool(args.get("log_only")),
    )


def _run_check_case(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    case_dir = of_writer.case_path(case_name)
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case '{case_dir.name}' not found"}
    result = precheck.summarise(precheck.check_case(case_dir))
    result["case_name"] = case_dir.name
    return result


def _run_run_case(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    case_dir = of_writer.case_path(case_name)
    if not case_dir.is_dir():
        return {
            "status": "error",
            "error": f"case '{case_dir.name}' not found; build it first with build_case",
        }
    background = args.get("background")
    background = True if background is None else bool(background)
    timeout = int(args.get("timeout_s") or runner.DEFAULT_TIMEOUT_S)
    skip_checks = bool(args.get("skip_checks"))
    if background:
        return runner.run_pipeline_background(
            case_dir, timeout=timeout, skip_checks=skip_checks
        )
    return runner.run_pipeline(case_dir, timeout=timeout, skip_checks=skip_checks)


def _run_run_status(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    case_dir = of_writer.case_path(case_name)
    if not case_dir.is_dir():
        return {"status": "error", "error": f"case '{case_dir.name}' not found"}
    return runner.run_status(case_dir)


def _run_change_settings(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    try:
        case_dir = of_writer.case_path(case_name)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}
    if not case_dir.is_dir():
        return {"status": "error",
                "error": f"case '{case_dir.name}' not found; build it first with build_case"}

    # Conditions and geometry are baked into the fields, thermochemistry and
    # mesh, so they cannot be edited in place: rebuild the case, keeping
    # everything else the user chose for it. Without this, a small model sent
    # "make the slot 6 mm" here, the value was silently ignored, and the tool
    # still answered "updated".
    physical = {k: v for k, v in validator.coerce_types(args).items()
                if k not in _RUN_CONTROLS}
    if physical:
        if runner.pids_in_dir(case_dir):
            return {"status": "error",
                    "error": f"case '{case_dir.name}' is running; stop it with "
                             "kill_jobs before changing its conditions or geometry"}
        case_type, _ = _chosen_parameters(case_dir)
        template = (case_templates.TEMPLATE_BY_NAME.get(case_type or "")
                    or case_templates.template_for_case_name(case_dir.name))
        if template is None:
            return {"status": "error",
                    "error": f"cannot tell which case type '{case_dir.name}' is"}
        result = _run_build_case({**args, "case_type": template.name,
                                  "case_name": case_dir.name})
        if result.get("status") == "built":
            result["note"] = (
                "rebuilt: " + ", ".join(physical) + " change the mesh, fields or "
                "thermochemistry, so the case was regenerated with the other "
                "settings kept; any previous results in it are gone"
            )
        return result

    sim_time = args.get("sim_time")
    sim_timestep = args.get("sim_timestep")
    sim_write_steps = args.get("sim_write_steps")
    n_steps, step_error = _n_steps(args)
    if step_error:
        return step_error
    if n_steps is not None:
        # The step is the new one if given, else the case's own.
        current = of_writer.read_control_entries(
            case_dir, ["deltaT", "endTime", "writeInterval"])
        try:
            dt = float(sim_timestep if sim_timestep is not None else current["deltaT"])
            writes = int(sim_write_steps) if sim_write_steps is not None else round(
                float(current["endTime"]) / float(current["writeInterval"]))
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            return {"status": "error",
                    "error": "cannot apply n_steps: deltaT/endTime/writeInterval not readable"}
        sim_time, sim_write_steps = _steps_to_controls(n_steps, dt, writes)
    if sim_time is None and sim_timestep is None and sim_write_steps is None:
        return {"status": "error",
                "error": "nothing to change; give sim_time, sim_timestep, sim_write_steps "
                         "and/or n_steps"}

    entries: dict[str, Any] = {}
    try:
        if sim_time is not None:
            sim_time = float(sim_time)
            if sim_time <= 0:
                return {"status": "error", "error": "sim_time must be > 0"}
            entries["endTime"] = sim_time
        if sim_timestep is not None:
            sim_timestep = float(sim_timestep)
            if sim_timestep <= 0:
                return {"status": "error", "error": "sim_timestep must be > 0"}
            entries["deltaT"] = sim_timestep
        if sim_write_steps is not None:
            sim_write_steps = int(sim_write_steps)
            if sim_write_steps <= 0:
                return {"status": "error", "error": "sim_write_steps must be a positive integer"}
            # writeInterval is a time, so it needs the (possibly new) endTime.
            end_for_write = sim_time
            if end_for_write is None:
                current = of_writer.read_control_entries(case_dir, ["endTime"]).get("endTime")
                end_for_write = float(current) if current is not None else None
            if end_for_write is None:
                return {"status": "error",
                        "error": "cannot set write steps: endTime not found in controlDict"}
            entries["writeInterval"] = end_for_write / sim_write_steps
    except (TypeError, ValueError) as exc:
        return {"status": "error", "error": f"invalid value: {exc}"}

    try:
        changes = of_writer.set_control_entries(case_dir, entries)
    except (KeyError, OSError) as exc:
        return {"status": "error", "error": str(exc)}
    return {
        "status": "updated",
        "case_name": case_dir.name,
        "changes": changes,
        "message": (
            "Updated " + ", ".join(f"{k} {v['old']} -> {v['new']}" for k, v in changes.items())
            + ". This takes effect on the next run or restart; the mesh, fields and "
            "thermochemistry are unchanged."
        ),
    }


def _run_restart_case(args: dict[str, Any]) -> dict[str, Any]:
    case_name = _default_case_name(args)
    try:
        case_dir = of_writer.case_path(case_name)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}
    if not case_dir.is_dir():
        return {"status": "error",
                "error": f"case '{case_dir.name}' not found; build it first with build_case"}
    if runner.pids_in_dir(case_dir):
        return {"status": "error",
                "error": f"case '{case_dir.name}' is currently running; stop it with "
                "kill_jobs before restarting"}

    # Validate restartability before touching controlDict.
    plan = runner.plan_restart(case_dir)
    if plan["status"] != "ready":
        return plan
    try:
        of_writer.set_control_entries(case_dir, {"startFrom": "latestTime"})
    except (KeyError, OSError) as exc:
        return {"status": "error", "error": str(exc)}
    timeout = int(args.get("timeout_s") or runner.DEFAULT_TIMEOUT_S)
    return runner.restart_pipeline_background(case_dir, timeout=timeout)


def _run_kill_jobs(args: dict[str, Any]) -> dict[str, Any]:
    names = args.get("case_names") or []
    if isinstance(names, str):
        names = [names]
    if not names:
        return {"status": "error", "error": "no case_names given"}

    results: list[dict[str, Any]] = []
    for name in names:
        try:
            case_dir = of_writer.case_path(name)
        except ValueError:
            results.append({"case": str(name), "status": "not_found"})
            continue
        if not case_dir.is_dir():
            results.append({"case": case_dir.name, "status": "not_found"})
            continue
        was_running = bool(runner.pids_in_dir(case_dir))
        killed = runner.kill_case(case_dir)
        results.append({
            "case": case_dir.name,
            "status": "killed" if was_running else "cancelled",
            "killed_pids": killed,
        })

    return {
        "status": "done",
        "results": results,
        "message": "; ".join(
            f"{r['case']}: {r['status']}"
            + (f" (pids {r['killed_pids']})" if r.get("killed_pids") else "")
            for r in results
        ),
    }


def _run_delete_cases(args: dict[str, Any]) -> dict[str, Any]:
    delete_all = bool(args.get("all"))
    names = args.get("case_names") or []
    if isinstance(names, str):
        names = [names]
    if delete_all:
        names = of_writer.generated_cases()
    if not names:
        return {
            "status": "noop",
            "message": "No generated cases to delete."
            if delete_all
            else "No case_names given.",
        }

    results: list[dict[str, Any]] = []
    for name in names:
        try:
            case_dir = of_writer.case_path(name)
        except ValueError:
            results.append({"case": str(name), "status": "refused"})
            continue
        if not case_dir.is_dir():
            results.append({"case": case_dir.name, "status": "not_found"})
            continue
        runner.kill_case(case_dir)  # stop it if running, then delete
        try:
            of_writer.delete_case(name)
            results.append({"case": case_dir.name, "status": "deleted"})
        except (ValueError, OSError) as exc:
            results.append({"case": case_dir.name, "status": "error", "reason": str(exc)})

    deleted = [r["case"] for r in results if r["status"] == "deleted"]
    return {
        "status": "done",
        "deleted_count": len(deleted),
        "results": results,
        "message": "; ".join(f"{r['case']}: {r['status']}" for r in results),
    }


def _run_run_cases(args: dict[str, Any]) -> dict[str, Any]:
    names = args.get("case_names") or []
    if isinstance(names, str):
        names = [names]
    if not names:
        return {"status": "error", "error": "no case_names given"}
    try:
        requested_parallel = int(args.get("max_parallel") or 1)
    except (TypeError, ValueError):
        requested_parallel = 1

    physical = runner.physical_core_count(
        config.load_config(require_key=False).physical_cores
    )

    cases: list[tuple[str, int]] = []
    not_found: list[str] = []
    for name in names:
        try:
            case_dir = of_writer.case_path(name)
        except ValueError:
            not_found.append(str(name))
            continue
        if not case_dir.is_dir():
            not_found.append(case_dir.name)
            continue
        cases.append((case_dir.name, runner.read_n_procs(case_dir)))

    plan = runner.plan_batch(cases, requested_parallel, physical)
    warnings = list(plan["warnings"])
    if not_found:
        warnings.append("not built (build them first; skipped): " + ", ".join(not_found))

    runnable_names = [name for name, _ in plan["runnable"]]
    skipped = [{"case": name, "n_procs": procs} for name, procs in plan["skipped"]]

    if not runnable_names:
        return {
            "status": "refused",
            "physical_cores": physical,
            "skipped": skipped,
            "not_found": not_found,
            "warnings": warnings,
            "message": "No cases could be run within the available physical cores.",
        }

    # Fresh start: clear any stale cancel markers from a previous batch.
    for name in runnable_names:
        (of_writer.case_path(name) / runner.CANCEL_MARKER).unlink(missing_ok=True)

    batch_log = REPO_ROOT / "agent" / "batch.log"
    batch_log.write_text(
        f"# FGMFoam batch: {len(runnable_names)} case(s), "
        f"{plan['max_parallel']} at a time on {physical} physical cores\n"
    )
    subprocess.Popen(
        [
            sys.executable, str(REPO_ROOT / "agent" / "run_batch.py"),
            "--max-parallel", str(plan["max_parallel"]),
            "--log", str(batch_log), *runnable_names,
        ],
        cwd=str(REPO_ROOT / "agent"),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )

    return {
        "status": "launched",
        "physical_cores": physical,
        "max_parallel": plan["max_parallel"],
        "running_sequentially": plan["max_parallel"] == 1,
        "running": runnable_names,
        "skipped": skipped,
        "not_found": not_found,
        "warnings": warnings,
        "batch_log": str(batch_log),
        "monitor": f"tail -f {batch_log}",
        "message": (
            f"Launched {len(runnable_names)} case(s) in the background, "
            f"{plan['max_parallel']} at a time"
            + (" (sequentially)." if plan["max_parallel"] == 1 else ".")
        ),
    }


def _run_compute_thermochemistry(args: dict[str, Any]) -> dict[str, Any]:
    phi = float(args["equivalence_ratio"])
    t_inlet = float(args.get("T_inlet", 300.0))
    p_0 = float(args.get("p_0", 101325.0))
    thermo = cantera_tool.compute_thermochemistry(phi, t_inlet, p_0)
    return {
        "equivalence_ratio": phi,
        "T_inlet": t_inlet,
        "p_0": p_0,
        "Yc_u": thermo.Yc_u,
        "Yc_b": thermo.Yc_b,
        "ht": thermo.ht,
    }


def _resolve_repo_file(path: str) -> Optional[Path]:
    """Resolve a requested file to a safe, readable path inside the repository."""
    path = (path or "").strip()
    if not path:
        return None
    candidates = [REPO_ROOT / path]
    # A bare case-relative path (e.g. "constant/combustionProperties") resolves
    # against each template case in turn, the default first.
    if not path.startswith(("solver/", "testcases/", "tables/", "agent/", "README")):
        candidates += [t.template_dir / path for t in case_templates.TEMPLATES]
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        if not resolved.is_relative_to(REPO_ROOT) or not resolved.is_file():
            continue
        rel = resolved.relative_to(REPO_ROOT).as_posix()
        if rel in _READABLE_EXACT_FILES or rel.startswith(_READABLE_DIR_ROOTS):
            return resolved
    return None


def _run_read_case_file(args: dict[str, Any]) -> dict[str, Any]:
    requested = args.get("path") or ""
    resolved = _resolve_repo_file(requested)
    if resolved is None:
        return {
            "status": "not_found",
            "message": (
                f"No readable repository file matches '{requested}'. Use a path "
                "like 'testcases/laminarSlotBurner/constant/combustionProperties', "
                "'0.orig/T', or 'solver/src/combustionModel/FGM/FGM/FGM.C'."
            ),
        }
    rel = resolved.relative_to(REPO_ROOT).as_posix()
    if resolved.name in _UNREADABLE_NAMES or resolved.suffix in _UNREADABLE_SUFFIXES:
        return {
            "status": "not_readable",
            "path": rel,
            "message": "That file is binary or mesh-geometry data and is not returned.",
        }
    if resolved.suffix == ".fgm":
        header = corpus._read_fgm_header(resolved) or ""
        return {
            "status": "ok",
            "path": rel,
            "content": header,
            "note": "FGM table header only; the numeric data after [DATA] is omitted.",
        }
    try:
        lines = resolved.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return {"status": "error", "path": rel, "message": str(exc)}
    truncated = len(lines) > _MAX_READ_LINES
    result: dict[str, Any] = {
        "status": "ok",
        "path": rel,
        "content": "\n".join(lines[:_MAX_READ_LINES]),
    }
    if truncated:
        result["note"] = f"truncated to the first {_MAX_READ_LINES} of {len(lines)} lines"
    return result


_DISPATCH = {
    "build_case": _run_build_case,
    "run_case": _run_run_case,
    "run_cases": _run_run_cases,
    "run_status": _run_run_status,
    "change_settings": _run_change_settings,
    "restart_case": _run_restart_case,
    "kill_jobs": _run_kill_jobs,
    "delete_cases": _run_delete_cases,
    "compute_thermochemistry": _run_compute_thermochemistry,
    "read_case_file": _run_read_case_file,
    "check_case": _run_check_case,
    "analyse_run": _run_analyse_run,
    "run_campaign": _run_run_campaign,
    "start_campaign": _run_start_campaign,
    "analyse_campaign": _run_analyse_campaign,
}


def execute_tool(name: str, args: dict[str, Any]) -> dict[str, Any]:
    """Run a tool by name, returning a dict (errors are returned, not raised)."""
    if name in SESSION_TOOLS:
        return {
            "error": (
                f"'{name}' changes the live session and must be handled by the "
                "chat agent, not the stateless tool layer"
            )
        }
    handler = _DISPATCH.get(name)
    if handler is None:
        return {"error": f"unknown tool '{name}'"}
    if name in _SINGLE_CASE_TOOLS and not args.get("case_name"):
        error = _no_case_named_error()
        if error is not None:
            return error
    try:
        return handler(args)
    except Exception as exc:  # noqa: BLE001 - surface failures to the model
        return {"error": f"{type(exc).__name__}: {exc}"}
