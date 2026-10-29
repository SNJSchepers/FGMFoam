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

"""Conversational, code-aware interface to the FGMFoam agent.

Each user turn retrieves relevant chunks from the RAG index, injects them into a
grounded system prompt, and calls the configured chat model (hosted Mistral or a
local OpenAI-compatible model, via `llm.make_client`) with the case-building tools
available. The model answers questions (grounded, with citations) or calls a tool
to compute thermochemistry / build the case. All validation and file writing stay
in Python (see `tools.py`).
"""

from __future__ import annotations

import json
import re
import time
from typing import Any, Callable, Optional

import case_templates
import tools as tools_module
import validator
import dataclasses

import config as config_module
import providers
from config import Config, ConfigError
from llm import make_client
from rag.index import Hit, VectorIndex, retrieve
from rag.index import _lexical_query_terms as query_terms

# Guard against a model that keeps calling tools without finishing a turn.
MAX_TOOL_ROUNDS: int = 5

# Retry of chat completions on transient API errors (rate limits, 5xx).
MAX_COMPLETION_ATTEMPTS: int = 6


def _param_lines() -> str:
    """Render the build_case parameter schema for the system prompt.

    Required means required by a template, as in the tool schema
    (``ParamSpec.required`` is only informational). Listing the run controls as
    required made careful models ask for values build_case fills in itself.
    """
    templates = case_templates.TEMPLATES
    required, optional = [], []
    for spec in validator.PARAM_SPECS:
        line = f"- {spec.name} [{spec.unit}]: {spec.description}"
        needed_by = [t.name for t in templates if spec.name in t.required]
        if needed_by:
            if len(needed_by) < len(templates):
                line += f" (only for {', '.join(needed_by)})"
            required.append(line)
        elif spec.name == "sim_timestep":
            optional.append(f"{line} (optional; computed from the mesh and flow "
                            "when omitted)")
        elif any(spec.name in t.defaults for t in templates):
            optional.append(f"{line} (optional; the case's recommended value "
                            "when omitted)")
        else:
            optional.append(f"{line} (optional, default {spec.default})")
    return "Required:\n" + "\n".join(required) + "\n\nOptional:\n" + "\n".join(optional)


SYSTEM_PROMPT = f"""You are a knowledgeable assistant for FGMFoam, an OpenFOAM \
solver for premixed laminar flames using Flamelet-Generated Manifold (FGM) \
chemistry. You help users understand the solver and set up its laminar slot \
burner case.

You can: (1) ANSWER questions about the solver, its libraries, and the example \
cases; (2) SET UP one of the template cases with build_case — laminarSlotBurner \
(2D slot burner) or laminarFlamelet (1D freely propagating flame); if the user \
has not made clear which one they want, ASK before building, never assume; (3) RUN and \
MANAGE cases — run_case (first run), restart_case (continue from the latest \
time), run_cases (a sweep, sequential unless the user asks for parallel), \
run_status, change_settings (edit controlDict without rebuilding), kill_jobs, \
delete_cases. Use compute_thermochemistry for thermochemistry numbers without \
building a case, and read_case_file to read a file's exact contents. Each tool's \
arguments and behaviour are described in its own schema — this prompt is the \
POLICY for when to use them and how to stay grounded.

ANSWERING vs. ACTING (read this first): If the user asks a QUESTION or for an \
EXPLANATION — e.g. "explain the solver", "what is ...", "how does ... work", \
"describe ...", "tell me about ...", "can you explain ..." — you MUST answer in \
text, grounded in the retrieved context (reading source files with read_case_file \
when needed). Do NOT call build_case, run_case, run_cases, restart_case, \
change_settings, kill_jobs, or delete_cases to answer a question. Those action \
tools are ONLY for explicit requests to set up, run, modify, or manage a case \
(e.g. "set up a case ...", "run it", "change the end time to ..."). In \
particular, NEVER call build_case unless the user explicitly asked to set up / \
build / create / configure a case AND has given parameter values — never call it \
to explain the solver or because a question is vague. When in doubt, answer in \
words.

Grounding rules:
- Base statements about THIS solver and its cases only on the retrieved context \
provided below, or on OpenFOAM facts you are certain of. When you use a fact from \
a source, name the file and line range in plain parentheses, e.g. \
(solver/applications/FGMFoam/createFields.H:80-85), using the labels shown in the \
context. Write every citation as ordinary text — never use bracketed reference \
markers or any structured citation format.
- Never invent or guess the contents of a file. When a question is about what a \
specific file contains (a dictionary such as combustionProperties, a 0/ field, a \
system/ dictionary, or a solver source file), call read_case_file to fetch its \
exact text and answer from that. Only describe a file's contents from text you \
have actually seen — the read_case_file result or the retrieved context — and if \
you do not have it, say so rather than filling it in from generic OpenFOAM \
knowledge.
- For how-it-works questions (an equation, boundary condition, class, the FGM \
model, the table lookup, the time loop), READ the source file(s) with \
read_case_file and describe from what you read — retrieved fragments are partial \
and you will get formulas/steps wrong. Key files: the PIMPLE time loop is \
solver/applications/FGMFoam/FGMFoam.C; the FGM model is FGM_PM_UL_NA.C; the table \
interpolation lookupFGM_2D is in solver/src/lookUp/FGM/FGMlib.c (called per cell \
in FGM_PM_UL_NA.C's update() and in the enthalpy BC \
solver/src/lookUpBoundaryCondition/lookUp/lookUpEnthalpyFvPatchScalarField.C). For \
"the lookup procedure" read FGMlib.c and the update() loop.
- When asked about "the algorithm", "how the solver works", or the "solution \
procedure" generally (NOT specifically the table lookup), read FGMFoam.C and LEAD \
with the per-time-step ORDER: (1) momentum predictor (UEqn.H); (2) \
combustion->correct() — the FGM model solves the two control-variable transport \
equations (progress variable Yc and enthalpy ht), then loops over cells to look up \
the FGM table and set temperature, compressibility psi (density), transport \
properties, species, heat release and the Yc source term; (3) the pressure \
corrector loop (pEqn.H / pcEqn.H), which sets rho = psi*p; (4) turbulence->correct() \
(a no-op for the laminar cases). Continuity is enforced THROUGH the pressure \
equation — do NOT list it as a separate equation, and treat rhoEqn.H as inactive. \
The table lookup (lookupFGM_2D) is ONE sub-step of combustion->correct(); focus on \
it only when explicitly asked (more in agent/knowledge/solver_notes.md).
- Never reproduce source code, equations, or dictionary blocks from memory or \
generic knowledge — you would get the details (coefficients, key names, values) \
wrong. Only show a verbatim snippet that is actually in your context (a \
read_case_file result or the retrieved context). When the user asks about a \
specific equation, boundary condition, or other code detail AND you have the file \
in context, DO show the short relevant snippet (e.g. the fvScalarMatrix that \
defines an equation) and then explain it — do not just describe it abstractly. \
For broader "how does it work" questions, prefer describing over dumping large \
blocks of code.
- OpenFOAM finite-volume operators, so you read the .C source correctly: \
fvm::ddt(rho, X) is the TRANSIENT term d(rho*X)/dt (implicit; fvc::ddt is its \
explicit form) — NOT convection. fvm::div(phi, X) or mvConvection().fvmDiv(phi, X) \
is CONVECTION, div(phi*X), where phi is the MASS FLUX (not the velocity), and \
mvConvection is a multivariate convection scheme (NOT a table lookup). \
fvm::laplacian(D, X) is DIFFUSION, div(D*grad(X)). A field on the RHS after == is \
an explicit SOURCE; fvm::Sp/fvm::SuSp are implicit sources. fvScalarMatrix \
assembles these into a linear system; .relax() under-relaxes and .solve() solves \
it. In FGMFoam the Yc/ht diffusivity is (alpha + turbulence().mut()/Sct_) = \
laminar (alpha) + turbulent (mut/Sct), and sourceYc_ is the tabulated reaction \
source (clamped to be non-negative).
- Answer the question that was asked and keep replies focused. Do not proactively \
list what a file or model does NOT contain (keys, parameters, equations) unless \
the user asks or you are correcting a misconception they actually stated.
- If the context does not cover something, say so plainly and point to the most \
relevant file rather than inventing behaviour.
- Distinguish what the SOLVER supports from how the example CASES are configured. \
The shipped cases (laminarSlotBurner, laminarFlamelet) are laminar CH4/air slot \
burners — their turbulenceProperties set simulationType laminar — but the solver \
is NOT laminar-only: it supports turbulent LES/DNS for a well-resolved flame (no \
sub-grid turbulence–chemistry-interaction model is implemented, so coarse LES and \
RANS are not accurate). Ground capability claims in the source/notes; do not claim \
turbulence is unsupported.
- Be honest about your own capabilities, especially regarding RESULTS. You do NOT \
read, analyse, post-process, plot, or interpret simulation results — flow fields, \
temperature/species/velocity profiles, convergence behaviour, or anything in the \
written time directories — because there is no tool that reads the numeric field \
data. For a run, what you CAN do is report its status with run_status: whether it \
is running or finished, the latest simulated time, and the tail of log.FGMFoam \
(useful for seeing whether it completed or errored). Do not offer to interpret, \
visualise, or analyse results, and do not list that as something you can do; if \
the user asks for it, say it is not currently supported and point them to ParaView \
or the standard OpenFOAM post-processing utilities. When you point to ParaView, \
keep it to one or two sentences (open the case with paraFoam, or create a .foam \
file and open it) plus the standard postProcess utility; do NOT write a \
step-by-step ParaView GUI/menu walkthrough or name specific ParaView \
buttons/menus/"widgets" or files to click — you will fabricate them.

Key facts about this project (always true; never contradict these, even if a \
file name or your general OpenFOAM knowledge suggests otherwise):
- The FGM lookup table is always constant/lookUp/database.fgm. No file configures \
a table path.
- constant/combustionProperties sets combustionModel FGM and exactly four \
FGMCoeffs switches — flameType, transportModel, adiabatic, subGridModel — that \
select the FGM model variant. It does NOT contain a table path, progress-variable \
bounds, or Prandtl/Schmidt coefficients; do not attribute other keys to it.
- constant/thermophysicalProperties is an unused placeholder required by the \
OpenFOAM framework; its values are overwritten by the FGM model, and it has no \
fuel/oxidiser entries. Do not infer the mixture from it — all thermochemistry \
comes from the FGM table.
- The operating pressure p_0 must match the FGM table's pressure (101325 Pa for \
the shipped table); pressure is baked into the table, so a different value is not \
supported without regenerating it. Do NOT build a case at a different pressure — \
build_case rejects other values.
- The FGM table is generated at a fixed equivalence ratio of 0.7; mixture \
fraction is not a table dimension. Set up cases at phi = 0.7 (the default) and do \
NOT build a case at a different equivalence ratio — the current table does not \
support it (it would need regenerating), and build_case rejects other values.
- The two FGM control variables (the table's two dimensions, CV1 and CV2) are the \
progress variable Yc and the enthalpy ht. The progress variable is defined as \
Yc = Y_CO2 / M_CO2 — the CO2 mass fraction divided by the molar mass of CO2 \
(≈ 44 kg/kmol), so its units are kmol/kg. It is based on CO2 ONLY, is NOT \
normalised to 0–1, and is NOT a sum of CO/CO2/H2/H2O mass fractions. Never state \
any other progress-variable definition; this is the definition used to generate \
the table and to set the case values. If unsure, read \
agent/knowledge/solver_notes.md rather than guessing.
- FGMFoam writes results as OpenFOAM TIME-STEP directories inside each case \
directory (for example testcases/laminarSlotBurner_T300/0.1/), reconstructed from \
the per-rank processor*/ subdirectories by reconstructPar at the end of the run. \
How many snapshots are saved is set by writeInterval = sim_time / sim_write_steps \
in system/controlDict (sim_write_steps defaults to 100, giving evenly-spaced \
snapshots — do NOT claim only the final time step is written unless \
sim_write_steps is 1). There is NO postProcessing/ directory: the functions block \
in the shipped controlDict is entirely commented out, so NEVER tell the user that \
results appear under postProcessing/ (that path only exists if function objects \
are enabled, which they are not). The per-case solver log is log.FGMFoam in the \
case directory.

build_case parameters:
{_param_lines()}

When the user wants to set up a case, collect any missing REQUIRED parameters \
with concise follow-up questions, then call build_case. Do not ask for OPTIONAL \
parameters: leave them out and build_case fills them in and reports which \
values it chose. Never guess unstated values. When the user has given every \
REQUIRED parameter, call build_case in that same turn: do not list the defaults \
and do not ask for confirmation first. Only build_case, though: even when the \
same message also says "and run it", the run waits for the user's next message. \
For a campaign, do not ask about cores or the CPU plan yourself: call \
run_campaign as soon as the case type and the swept values are clear, leaving \
n_procs_per_case and n_concurrent_cases out. It then returns the machine's core \
count and the sensible options, and that is the question to ask. \
After a successful build, briefly summarise what was written and where, \
and offer to run it.

Building and running are SEPARATE steps. build_case only prepares the case; it \
does not start the solver. NEVER call run_case in the same turn as build_case or \
change_settings: after a build or a settings change you MUST stop and report \
(case ready / what changed) and merely OFFER to run — do not run it yourself. \
Only call run_case when the user's message ITSELF asks to run, execute, or launch \
the case. IMPORTANT: a duration such as "run 0.01 s", "run for 0.5 s", or "run \
0.1 s" specifies the SIMULATED TIME (the sim_time parameter) — it is part of a \
setup request, NOT a request to execute the solver; do not call run_case because \
the word "run" appears next to a time. When the user does explicitly ask to run, \
call run_case (it launches in the background), then tell them the log file to tail \
and that they can ask you to check progress with run_status. The same \
explicit-request rule applies to restart_case — only restart/continue/extend when \
asked. change_settings only edits an existing case's controlDict (endTime / deltaT \
/ write snapshots) and never starts the solver; after changing a setting, report \
what changed and that a run or restart is needed to apply it."""


def format_context(hits: list[Hit]) -> str:
    """Render retrieved chunks into a citable context block."""
    if not hits:
        return "(no relevant source files matched this query)"
    blocks = [
        f"Source {i} — {hit.chunk.citation} "
        f"(kind={hit.chunk.kind}, score={hit.score:.2f}):\n{hit.chunk.text}"
        for i, hit in enumerate(hits, 1)
    ]
    return "\n\n".join(blocks)


# Cues that the user is asking how something is implemented / how it works. When
# one of these appears AND a source file was retrieved, the agent is forced to
# read the source (rather than answer from partial excerpts).
#
# A retrieval-only gate ("force when the top retrieved chunk is C/C++") was tried
# and does not work for this corpus: the README and knowledge notes describe the
# internals in prose and outrank the actual .C/.H source for "how does X work"
# questions, so source is rarely the top hit; meanwhile a source chunk appears
# somewhere in nearly every query's top-k (including setup/run), so "any source
# in the hits" alone over-fires. The keyword intent signal, AND-ed with "a source
# file is somewhere in the hits", separates implementation questions from
# setup/run requests reliably here.
_IMPLEMENTATION_CUES: tuple[str, ...] = (
    "how does", "how is", "how are", "how the", "how it works",
    "implement", "implementation", "under the hood", "mechanism",
    "algorithm", "step by step", "step-by-step", "exactly how",
    "walk me through", "walk through", "go through the code", "source code",
    "explain", "describe", "detail", "tell me about", "procedure",
    "break down", "break it down", "lookup", "look up", "look-up",
    # Asking about an equation/term should read the source, not paraphrase it
    # (the model would otherwise risk fabricating the equation).
    "equation", "transport equation", "source term", "discretis", "discretiz",
)


def _wants_implementation(text: str) -> bool:
    lowered = text.lower()
    return any(cue in lowered for cue in _IMPLEMENTATION_CUES)


# "Separate steps" guard: preparing a case (build/settings change) and running it
# are deliberately distinct turns. Within ONE turn, once a prepare tool has run, a
# run/restart tool is NOT executed — it is deferred so the solver can never launch
# unprompted just because a (small) model chained the calls. Running then needs an
# explicit follow-up message. Enforced in code, so it holds for any model/provider.
_PREPARE_TOOLS: frozenset[str] = frozenset({"build_case", "change_settings"})
_RUN_TOOLS: frozenset[str] = frozenset(
    {"run_case", "run_cases", "restart_case", "start_campaign"}
)
_DEFERRED_RUN_RESULT: dict[str, Any] = {
    "status": "deferred",
    "message": (
        "Not run. Building/changing a case and running it are SEPARATE steps and "
        "cannot happen in the same turn. Tell the user the case is ready (or the "
        "settings are updated) and that they can run it by asking in a new message."
    ),
}

# "Feasibility-question" guard: a question about whether something is POSSIBLE or
# SUPPORTED ("is it possible to run at 10 bar?", "can the solver do X?") must be
# ANSWERED, never executed. A small model otherwise routes such a question to an
# action tool AND fabricates the missing parameters — verified: a 10-bar
# feasibility question made the 3B call build_case with invented velocity/
# temperatures/times, build a case at the DEFAULT pressure, then claim it was 10
# bar. When the current message matches a feasibility marker, every action tool is
# deferred (nothing is built/run/changed/deleted) and the model answers in words
# from the retrieved constraints. Enforced in code, so it holds for any model.
# Markers are info-seeking phrasings chosen NOT to fire on polite imperatives
# ("can you run it?", "can you set up a case ...?") — those carry real intent.
# switch_model is deliberately NOT here. The guard exists to stop a model
# fabricating physical parameters and then building/running/deleting a case;
# switching the answering model invents nothing, touches no case, and is undone
# by switching back. Guarding it only broke compound messages like "which models
# can I use? then switch to the local one", where the marker matched the
# question and blocked the explicit instruction beside it.
_ACTION_TOOLS: frozenset[str] = (
    _PREPARE_TOOLS | _RUN_TOOLS | frozenset({"kill_jobs", "delete_cases"})
)
_FEASIBILITY_MARKERS: tuple[str, ...] = (
    "is it possible", "is that possible", "would it be possible", "is there a way",
    "possible to", "am i able", "are we able", "is this supported",
    "is it supported", "does it support", "do you support", "can it ",
    "can the solver", "can this solver", "can i run", "can i use", "can i set",
    "can i build", "can i change", "can i do", "what happens if", "is it ok to",
    "is it okay to", "am i allowed", "are we allowed",
)
_DEFERRED_FEASIBILITY_RESULT: dict[str, Any] = {
    "status": "not_performed",
    "message": (
        "This message asks whether something is POSSIBLE or SUPPORTED — it is a "
        "question, not a request to perform the action. Do NOT build, run, change, "
        "or delete anything, and never invent parameter values the user did not "
        "give. Answer in words from the retrieved context (e.g. the pressure and "
        "equivalence-ratio constraints), and only act if the user explicitly asks "
        "in a later message."
    ),
}


def _is_feasibility_question(text: str) -> bool:
    lowered = text.lower()
    return any(marker in lowered for marker in _FEASIBILITY_MARKERS)


# "Do not start" guard: "prepare it, but do not start it yet" must hold even when
# a small model calls run_cases anyway or leaves launch unset on run_campaign
# (verified with the 3B: both happened, with invented case names). When the
# current message says not to start, every run tool is deferred and run_campaign
# is forced to launch=false. Enforced in code, so it holds for any model.
_DONT_START_MARKERS: tuple[str, ...] = (
    "do not start", "don't start", "dont start", "do not run", "don't run",
    "dont run", "do not launch", "don't launch", "dont launch", "not start it",
    "not run it", "without starting", "without running", "without launching",
    "only prepare", "just prepare", "prepare only", "only build", "just build",
    "build only", "only set it up", "just set it up", "wait for my", "until i say",
)
_DEFERRED_DONT_START_RESULT: dict[str, Any] = {
    "status": "not_started",
    "message": (
        "Not run: the user said NOT to start anything in this message. Report what "
        "was prepared and that it can be started with a later, explicit request."
    ),
}


def _forbids_starting(text: str) -> bool:
    lowered = text.lower().replace("’", "'")
    return any(marker in lowered for marker in _DONT_START_MARKERS)


# Case-type confirmation: a model may only build (or sweep) a case type that the
# conversation actually asked for — by template name, a plugin keyword, or a
# parameter unique to it. Otherwise the tool is not run and the model is told to
# ask. build_case already asks when case_type is MISSING; this covers a case type
# the model INVENTED (the 3B turned "a simulation at 350 K" into a slot burner).
_CASE_TYPE_TOOLS: frozenset[str] = frozenset({"build_case", "run_campaign"})


def _requested_template(name: str, args: dict[str, Any]) -> Any:
    """The template a build_case/run_campaign call would act on, if it names one."""
    case_type = args.get("case_type")
    if case_type:
        try:
            return case_templates.find_template(str(case_type))
        except ValueError:
            return None
    if name == "build_case" and args.get("case_name"):
        return case_templates.template_for_case_name(str(args["case_name"]))
    return None


def _case_type_from_user(user_input: str, said: str) -> Optional[str]:
    """The one case type the user named — in this message, else in the whole
    conversation — or None if they named none or several."""
    for text in (user_input, said):
        named = [t for t in case_templates.TEMPLATES if t.mentioned_in(text)]
        if len(named) == 1:
            return named[0].name
        if named:
            return None
    return None


def _variable_from_user(template_name: str, user_input: str) -> Optional[str]:
    """The sweepable parameter the user's message points at, or None.

    Either the message names exactly one ("plate thickness" -> plate_thickness),
    or its values carry a unit that exactly one sweepable parameter has
    ("300, 400 and 500 K" -> T_inlet for the flamelet). Never a guess from the
    case alone: "equivalence ratios 0.6, 0.7, 0.8" has neither, so it stays
    unset and the tool asks — it must not become a temperature sweep.
    """
    template = case_templates.TEMPLATE_BY_NAME.get(template_name)
    if template is None:
        return None
    # Strongest signal: the parameter whose name is followed by the value list
    # ("sweeps the plate thickness over 2, 3 and 4 mm"). Other parameters named
    # with single values in the same message are the fixed conditions.
    # ("plate" belongs to both T_plate and plate_thickness: the longest phrase
    # that precedes a list decides.)
    with_list: dict[str, int] = {}
    for variable in template.sweepable:
        spec = validator.SPEC_BY_NAME.get(variable)
        phrases = {variable.replace("_", " "), variable, *(spec.phrases if spec else ())}
        for phrase in phrases:
            if re.search(r"(?i:\b" + re.escape(phrase) + r"\b)[^\d\n]{0,30}?"
                         r"\d+(?:\.\d+)?\s*(?:mm|cm|m/s|m|K|s)?\s*(?:,|and|or|to)\s*\d",
                         user_input):
                with_list[variable] = max(with_list.get(variable, 0), len(phrase))
    if with_list:
        longest = max(with_list.values())
        best = [v for v, n in with_list.items() if n == longest]
        if len(best) == 1:
            return best[0]
    lowered = user_input.lower()
    named = [v for v in template.sweepable
             if v.replace("_", " ").lower() in lowered or v.lower() in lowered]
    if len(named) == 1:
        return named[0]
    if named:
        return None
    units = {unit for _, unit in _VALUE.findall(user_input)}
    if _UNIT_LIST.search(user_input):
        units |= {"mm", "cm"} & set(_UNIT_LIST.findall(user_input)[0])
    by_unit = [v for v in template.sweepable
               if (spec := validator.SPEC_BY_NAME.get(v)) is not None
               and (spec.unit in units or (spec.unit == "m" and units & {"mm", "cm"}))]
    return by_unit[0] if len(by_unit) == 1 else None


# Announced but not made: "I will call the change_settings function ..." with no
# call at all. Verified with the 3B. Treated like a false claim (one nudge).
_PROMISED_CALL = re.compile(
    r"\bI(?:'ll| will| am going to)\s+(?:now\s+)?(?:call|use|run|invoke|build|rebuild|set up|prepare|create|"
    r"change|update|start|launch|sweep)\b"
    r"|\b(?:call|calling|invoke|invoking)\s+the\s+`?\w+`?\s+(?:function|tool)\b",
    re.IGNORECASE,
)


# Claim check: a reply that says a case was built, or a run started, when no tool
# did so this turn is corrected in place. Verified with the 3B: it answered
# "rebuild it at 450 K" with only a thermochemistry call and replied that the case
# "has been set up and built successfully".
_BUILD_CLAIM = re.compile(
    r"(has|have) been (successfully )?(re)?(built|set up|created|prepared|updated)"
    r"|successfully (re)?(built|set up|created|prepared)"
    # "ready to be built" / "ready to go" promise a build, they do not claim one.
    r"|is (now )?(ready|built|set up)\b(?!\s+to\b)",
    re.IGNORECASE,
)
_RUN_CLAIM = re.compile(
    r"(is|are) (now |currently )?running|has been (started|launched)"
    r"|(simulation|campaign|run|solver) (has )?(started|launched)",
    re.IGNORECASE,
)
_BUILD_REQUEST = re.compile(
    r"\b(set ?up|build|rebuild|create|make|change|prepare|modify|update)\b", re.IGNORECASE
)
_RUN_SUCCESS = {"started", "running", "launched", "queued", "completed", "ok"}
_CLAIM_SUBJECT = re.compile(r"\b(case|cases|simulation|campaign|run|solver)\b", re.IGNORECASE)


# A negation or a condition just before the match means it is not a claim:
# "nothing has been built or started yet", "nothing is running", "before
# anything is built", "once it is set up".
_NEGATED = re.compile(r"\b(?:nothing|not|no|never|none|before|until|once|if|unless)\b|n't",
                      re.IGNORECASE)


def _claims(pattern: re.Pattern[str], reply: str) -> bool:
    """True when a sentence about a case/simulation/campaign matches ``pattern``
    ("the FGM table has been created with Chem1D" is not a claim about a case),
    unless the match is negated."""
    for sentence in re.split(r"(?<=[.!?])\s+|\n+", reply):
        if not _CLAIM_SUBJECT.search(sentence):
            continue
        for match in pattern.finditer(sentence):
            if not _NEGATED.search(sentence[max(0, match.start() - 25):match.start()]):
                return True
    return False


def _claim_correction(
    user_input: str, reply: str, outcomes: list[tuple[str, str, dict[str, Any]]]
) -> Optional[str]:
    """A correction to append when the reply claims what no tool did, else None.

    ``outcomes`` are (tool name, status, result) for every call this turn.
    """
    built = any(n in ("build_case", "change_settings") and st in ("built", "ok", "updated")
                or n == "run_campaign" and st == "prepared"
                for n, st, _ in outcomes)
    ran = any(n in _RUN_TOOLS | {"run_campaign", "run_status"} and st in _RUN_SUCCESS
              for n, st, _ in outcomes)
    problems = [
        (n, st, r) for n, st, r in outcomes
        if n in _ACTION_TOOLS | {"run_campaign"} and st not in ("built", "prepared", "ok")
    ]
    lines = []
    if (_claims(_BUILD_CLAIM, reply) and not built
            and (problems or _BUILD_REQUEST.search(user_input))):
        lines.append("no case was built or changed in this step")
    if _claims(_RUN_CLAIM, reply) and not ran:
        lines.append("nothing was started in this step")
    if not lines:
        return None
    detail = "; ".join(
        f"{n} returned '{st}'"
        + (f" (missing: {', '.join(r['missing'])})" if r.get("missing") else "")
        for n, st, r in problems
    )
    return (
        "\n\n> ⚠️ **Correction from the agent:** " + " and ".join(lines) + "."
        + (f" The tools reported: {detail}." if detail else "")
        + " The text above is wrong on that point."
    )


# Requested-value check: a value the user gave with a unit ("450 K", "5 mm") must
# reach the tool call that acts on it. Verified with the 3B: asked to rebuild at
# 450 K, it called change_settings with the OLD run controls, dropped the 450 K,
# and replied "updated as requested". When an action tool ran this turn but no
# argument carries a requested value, the model is told once, before it answers,
# so it can make the call it missed.
_VALUE = re.compile(
    r"(?<![\w.])(\d+(?:\.\d+)?(?:e-?\d+)?)\s*(K|mm|cm|m/s|m|s|bar|Pa)(?![\w/])"
)
_TO_SI: dict[str, float] = {"mm": 1e-3, "cm": 1e-2, "bar": 1e5}
_VALUE_TOOLS: frozenset[str] = frozenset({"build_case", "change_settings", "run_campaign"})


def _numbers(value: Any) -> list[float]:
    if isinstance(value, bool):
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, str):
        try:
            return [float(value)]
        except ValueError:
            return []
    if isinstance(value, list):
        return [n for item in value for n in _numbers(item)]
    if isinstance(value, dict):
        return [n for item in value.values() for n in _numbers(item)]
    return []


def _unapplied_values(user_input: str, calls: list[tuple[str, dict[str, Any]]]) -> list[str]:
    """Values with units in the message that no action-tool argument carries."""
    acted = [args for name, args in calls if name in _VALUE_TOOLS]
    if not acted:
        return []
    used = [n for args in acted for n in _numbers(args)]
    missing = []
    for number, unit in _VALUE.findall(user_input):
        target = float(number) * _TO_SI.get(unit, 1.0)
        if not any(abs(n - target) <= 0.01 * max(abs(target), 1e-12) for n in used):
            missing.append(f"{number} {unit}")
    return missing


# A knowledge note that ranks above every source chunk and names a source file
# decides which file the forced read opens. Code is sometimes hard for a small
# embedding model to match ("table lookup out of range" vs interpolation code
# full of eta1/w1p); the note is easy to match, and naming the file in it is how a
# note steers grounding — the same mechanism that corrects the agent's beliefs.
_SOURCE_PATH = re.compile(r"`(solver/[\w/.-]+\.(?:C|H|c|h))`")


_NOTE_HEADING = re.compile(r"^#+\s+(.+)$", re.MULTILINE)


def _source_named_by_notes(hits: list[Hit], question: str) -> Optional[str]:
    """The source file the top-ranked notes point to for this question.

    Only notes ranked above every source chunk count. Of those that name a file,
    the one whose section heading shares the most of the question's words decides,
    and a tie keeps the better-ranked one. Ranking alone is not enough: "how does
    the table lookup work, do the Yc bounds vary with enthalpy?" ranked the
    enthalpy boundary-condition note (which mentions the lookup) above the note
    headed "Table lookup", and the forced read opened the boundary condition.
    A note often names several files, so the one whose source contains the most
    of the question's words wins (a tie keeps the earlier one).
    """
    words = set(query_terms(question))
    best: tuple[int, list[str]] = (-1, [])
    for hit in hits:
        if hit.chunk.kind == "cpp":
            break
        named = [path for path in dict.fromkeys(_SOURCE_PATH.findall(hit.chunk.text))
                 if (config_module.REPO_ROOT / path).is_file()]
        if not named:
            continue
        heading = _NOTE_HEADING.search(hit.chunk.text)
        score = len(words & set(query_terms(heading.group(1)))) if heading else 0
        if score > best[0]:
            best = (score, named)
    named = best[1]
    if not named:
        return None

    def overlap(path: str) -> int:
        text = (config_module.REPO_ROOT / path).read_text(errors="replace").lower()
        return sum(word in text for word in words)

    return max(named, key=overlap)


# Invented values: a REQUIRED parameter (one the user must choose, such as T_plate
# or inlet_velocity) is accepted only if its value appears somewhere in what the
# user wrote — directly or after a unit conversion (5 mm -> 0.005). Verified with
# the 3B: build_case answered "missing T_inlet, T_plate", and in the same turn the
# model called it again with T_inlet=300 and T_plate=300 of its own making.
# Values inherited by a rebuild are not passed by the model, so are not affected.
_SCALES: tuple[float, ...] = (1.0, 1e-3, 1e-2, 1e5)
_ANY_NUMBER = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?(?:e-?\d+)?)")


def _invented_required(name: str, args: dict[str, Any], said: str) -> list[str]:
    """Required parameters whose value the user never wrote."""
    if name not in ("build_case", "run_campaign"):
        return []
    template = _requested_template(name, args) or case_templates.infer_template(
        validator.canonical_name(k) for k in args)
    if template is None:
        return []
    swept = validator.canonical_name(str(args.get("variable") or ""))
    given = [float(n) * scale for n in _ANY_NUMBER.findall(said) for scale in _SCALES]
    invented = []
    for key, value in args.items():
        key = validator.canonical_name(key)
        if key not in template.required or key == swept:
            continue
        numbers = _numbers(value)
        if numbers and not all(
            any(abs(n - g) <= 0.01 * max(abs(n), 1e-12) for g in given) for n in numbers
        ):
            invented.append(key)
    return invented


_UNIT_LIST = re.compile(
    r"((?:\d+(?:\.\d+)?\s*(?:,|and|or|to)\s*)+\d+(?:\.\d+)?)\s*(mm|cm)\b"
)


def _convert_stated_units(args: dict[str, Any], said: str) -> list[str]:
    """Convert length arguments the model passed in the unit the USER wrote.

    "a plate thickness of 5 mm" -> the model passes plate_thickness=5, which
    tools.py would refuse as out of range. When the user's own text has that
    number with a length unit, the conversion is unambiguous, so it is done here
    (in place) and reported. Returns the descriptions of what was converted.
    """
    stated: dict[float, tuple[str, float]] = {}
    for number, unit in _VALUE.findall(said):
        if unit in ("mm", "cm"):
            stated.setdefault(float(number), (unit, _TO_SI[unit]))
    # "2, 3 and 4 mm": the unit written once applies to the whole list.
    for numbers, unit in _UNIT_LIST.findall(said):
        for number in re.findall(r"\d+(?:\.\d+)?", numbers):
            stated.setdefault(float(number), (unit, _TO_SI[unit]))
    lengths = {spec.name for spec in validator.PARAM_SPECS if spec.unit == "m"}
    converted = []
    for key in list(args):
        if validator.canonical_name(key) not in lengths and not (
            key == "values"
            and validator.canonical_name(str(args.get("variable") or "")) in lengths
        ):
            continue
        value = args[key]
        items = value if isinstance(value, list) else [value]
        nums = _numbers(items)
        if len(nums) != len(items) or not nums or all(
            abs(n) <= tools_module.MAX_LENGTH_M for n in nums
        ):
            continue
        if not all(n in stated for n in nums):
            continue
        new = [n * stated[n][1] for n in nums]
        args[key] = new if isinstance(value, list) else new[0]
        converted.append(f"{key} {value} {stated[nums[0]][0]} -> {args[key]} m")
    return converted


# Stated values: "a plate thickness of 5 mm", "the slot (hole size) 6 mm wide",
# "an unburnt temperature of 450 K". A parameter phrase (ParamSpec.phrases or its
# name with spaces), at most a few filler words, then a number and a unit that
# fits the parameter. Verified with the 3B: told "rebuild with a plate thickness
# of 5 mm", it called change_settings with only the old run controls three times,
# nudges included. Python makes sure the call carries what the user stated.
_UNIT_SI: dict[str, tuple[str, float]] = {
    "mm": ("m", 1e-3), "cm": ("m", 1e-2), "m": ("m", 1.0), "K": ("K", 1.0),
    "m/s": ("m/s", 1.0), "s": ("s", 1.0), "ms": ("s", 1e-3),
}
_STATED_TOOLS: frozenset[str] = frozenset({"build_case", "change_settings"})


def _stated_parameters(text: str) -> dict[str, float]:
    """{parameter: SI value} for every 'phrase ... number unit' in ``text``."""
    found: dict[str, tuple[int, float]] = {}
    for spec in validator.PARAM_SPECS:
        phrases = {spec.name.replace("_", " "), spec.name, *spec.phrases}
        for phrase in sorted(phrases, key=len, reverse=True):
            # Phrase case-insensitive; unit case-sensitive (K is kelvin).
            pattern = re.compile(
                r"(?i:\b" + re.escape(phrase) + r"\b)[^\d\n.;]{0,25}?"
                r"(\d+(?:\.\d+)?(?:e-?\d+)?)\s*(mm|cm|m/s|ms|m|K|s)(?![\w/])"
            )
            for match in pattern.finditer(text):
                unit, scale = _UNIT_SI[match.group(2)]
                if unit != spec.unit:
                    continue
                # The longest (most specific) phrase decides a contested value.
                if spec.name not in found or len(phrase) > found[spec.name][0]:
                    found[spec.name] = (len(phrase), float(match.group(1)) * scale)
    return {name: value for name, (_, value) in found.items()}


# CPU-plan answers. run_campaign asks how to split the cores and stops; the user
# answers in the next message. The 3B then made no call at all — three times out
# of three, nudge included. The answer completes a call the tool itself asked
# for, so the chat loop makes it (like the forced source read) and the model
# only reports the result.
_PROCS_PER_CASE = re.compile(
    r"(\d+)\s*(?:mpi\s+)?(?:processors?|procs?|cores?|cpus?|ranks?|processes)\s*"
    r"(?:per|for each|each|a)\s*case", re.IGNORECASE)
_CONCURRENT = re.compile(
    r"(\d+)\s*(?:cases?\s*)?(?:at (?:the same|a|one) time|at once|concurrent(?:ly)?"
    r"|in parallel|simultaneous(?:ly)?)", re.IGNORECASE)


def _cpu_plan(text: str) -> Optional[tuple[int, int]]:
    procs, concurrent = _PROCS_PER_CASE.search(text), _CONCURRENT.search(text)
    if not procs or not concurrent:
        return None
    return int(procs.group(1)), int(concurrent.group(1))


# Context budget for local models (see providers.context_tokens). Characters per
# token for this mix of prose, C++ and JSON, chosen on the conservative side
# (measured ~4.3 for the system prompt + tool schemas).
_CHARS_PER_TOKEN: float = 3.5
# Tokens kept free for the model's answer.
_ANSWER_RESERVE: int = 1500


def _estimate_tokens(*parts: Any) -> int:
    return int(sum(len(p if isinstance(p, str) else json.dumps(p, ensure_ascii=False))
                   for p in parts) / _CHARS_PER_TOKEN)


def _tool_content(result: dict[str, Any]) -> str:
    """A tool result as the model reads it.

    Non-ASCII stays as it is: the default escaping turned every source file's
    OpenFOAM banner (box-drawing characters) into ``\\u2588`` sequences, several
    tokens per character, on every read.
    """
    return json.dumps(result, ensure_ascii=False)


def _densest_window(lines: list[str], question: str, budget_chars: int) -> list[tuple[int, int]]:
    """The run of lines, about ``budget_chars`` long, with the most question words.

    Used when a file is read without retrieved line ranges (a note named it), so
    the excerpt shows e.g. update() for a question about the update step, rather
    than the licence header and includes at the top of the file.
    """
    words = [w for w in query_terms(question) if len(w) > 3]
    if not words or not lines:
        return []
    hits = [sum(w in line.lower() for w in words) for line in lines]
    average = max(1, sum(len(line) + 1 for line in lines) // len(lines))
    span = max(20, min(len(lines), budget_chars // average))
    best, best_start = -1, 0
    total = sum(hits[:span])
    for start in range(0, len(lines) - span + 1):
        if start:
            total += hits[start + span - 1] - hits[start - 1]
        if total > best:
            best, best_start = total, start
    return [(best_start + 1, best_start + span)] if best > 0 else []


def _fit_file_read(result: dict[str, Any], budget_tokens: int,
                   focus: list[tuple[int, int]], question: str = "") -> dict[str, Any]:
    """Shrink a read_case_file result to ``budget_tokens``.

    Keeps the line ranges in ``focus`` (the retrieved chunks of this file) with
    some margin; without them, the part of the file densest in the question's
    words; else the start of the file. The note says what was shown.
    """
    content = result.get("content")
    if not isinstance(content, str) or _estimate_tokens(content) <= budget_tokens:
        return result
    lines = content.splitlines()
    budget_chars = int(max(budget_tokens, 200) * _CHARS_PER_TOKEN)
    if not focus and question:
        focus = _densest_window(lines, question, budget_chars)
    windows = sorted({(max(1, a - 15), min(len(lines), b + 15)) for a, b in focus}) or [
        (1, len(lines))
    ]
    kept: list[str] = []
    used = 0
    shown: list[str] = []
    for start, end in windows:
        block: list[str] = []
        for number in range(start, end + 1):
            line = lines[number - 1]
            if used + len(line) + 1 > budget_chars:
                break
            block.append(line)
            used += len(line) + 1
        if block:
            last = start + len(block) - 1
            kept.append(f"// --- lines {start}-{last} ---")
            kept.extend(block)
            shown.append(f"{start}-{last}")
        if used >= budget_chars:
            break
    fitted = dict(result)
    fitted["content"] = "\n".join(kept)
    fitted["note"] = (
        f"excerpt (lines {', '.join(shown)} of {len(lines)}) — the whole file does "
        "not fit this model's context window"
    )
    return fitted


def _message_text(message: Any) -> str:
    """Return a plain string whether content is a string or a list of chunks."""
    content = getattr(message, "content", None)
    if content is None:
        return ""
    if isinstance(content, str):
        return content
    parts: list[str] = []
    for chunk in content:
        text = getattr(chunk, "text", None)
        if text is None and isinstance(chunk, dict):
            text = chunk.get("text")
        if text:
            parts.append(text)
    return "".join(parts)


class ChatAgent:
    """Holds conversation state and runs one grounded, tool-enabled turn."""

    def __init__(
        self,
        config: Config,
        index: VectorIndex,
        on_tool_call: Optional[Callable[[str, dict[str, Any]], None]] = None,
        on_event: Optional[Callable[[str], None]] = None,
        on_guard: Optional[Callable[[str, str, str, bool], None]] = None,
    ) -> None:
        self.config = config
        self.index = index
        self.client = make_client(config)
        # Active chat model — mutable so it can be switched live (e.g. via the
        # REPL's /model command); use_config() below swaps the provider too.
        self.chat_model = config.chat_model
        tools_module.set_session_model(config.provider, config.chat_model)
        # A new conversation has no "last case" to default to or rebuild from.
        tools_module.reset_session_cases()
        self.history: list[dict[str, Any]] = []
        # Optional diagnostic hooks so the user can see what the agent is doing:
        # on_tool_call(name, args) fires before each tool runs; on_event(message)
        # reports other activity (retrieval, forced reads).
        self._on_tool_call = on_tool_call
        self._on_event = on_event
        # on_guard(kind, name, message, acted) fires whenever a guard detects
        # something, so an eval can count failures per guard.
        self._on_guard = on_guard

    def _emit(self, message: str) -> None:
        if self._on_event is not None:
            self._on_event(message)

    def _guard(self, kind: str, name: str, message: str) -> bool:
        """Report a guard that detected a problem; return True if it may act.

        ``kind`` is "repair" (the model's call runs with corrected arguments),
        "block" (the tool does not run; the model gets a substitute result) or
        "nudge" (the draft reply is retried once, or corrected). With
        FGM_GUARDS=shadow the detection is only logged and nothing is changed,
        so an eval measures what the model itself does.
        """
        acted = self.config.guards == "on"
        self._emit(message if acted else f"[shadow, not applied] {message}")
        if self._on_guard is not None:
            self._on_guard(kind, name, message, acted)
        return acted

    def switch_model(self, choice: str = "") -> dict[str, Any]:
        """Switch the model answering this conversation, or list the options.

        The single implementation behind both the conversational ``switch_model``
        tool and the REPL's ``/model`` command, so the two can never diverge.
        Resolution and every check stay in Python — the model only decides *when*
        to switch, exactly as with case building.

        Returns a result dict; failures are returned, never raised, so a bad
        model name ends the turn with an explanation rather than an exception.
        """
        available = [
            {"name": preset.name, "description": preset.summary}
            for preset in providers.PRESETS
        ]
        if not choice.strip():
            return {
                "status": "current",
                "current": config_module.describe(self.config),
                "available": available,
            }

        preset, bare_model = providers.split_choice(choice)
        try:
            if preset is None and bare_model:
                # A bare model id changes ONLY the model, keeping the current
                # provider, endpoint and key — re-resolving from scratch would
                # fall back to the default provider and ask for a key the user
                # does not need (e.g. "qwen2.5-3b-16k" while running locally).
                problem = providers.implausible_model(self.config.provider, bare_model)
                if problem:
                    raise ConfigError(problem)
                new_config = dataclasses.replace(self.config, chat_model=bare_model)
            else:
                new_config = config_module.load_config(model=choice)
        except ConfigError as exc:
            return {
                "status": "error",
                "requested": choice,
                "current": config_module.describe(self.config),
                "available": available,
                "error": str(exc),
            }

        # The stored vectors belong to one embedding model: switching to a chat
        # model that would need a different one is refused rather than silently
        # invalidating the index mid-conversation.
        if new_config.embed_model != self.config.embed_model:
            return {
                "status": "error",
                "requested": choice,
                "current": config_module.describe(self.config),
                "error": (
                    f"{choice!r} needs embedding model "
                    f"{new_config.embed_model!r}, but the knowledge index was "
                    f"built with {self.config.embed_model!r}. Rebuild it first: "
                    "agent/.venv/bin/python agent/index_build.py --force"
                ),
            }

        previous = config_module.describe(self.config)
        self.use_config(new_config)
        result = {
            "status": "switched",
            "previous": previous,
            "current": config_module.describe(new_config),
            "note": "Applies from the next reply, for this session only; the "
                    "conversation so far is kept.",
        }
        # A name that matches no preset is passed through as a raw model id, so
        # that models released after this code was written stay reachable. That
        # also means a typo cannot be caught here — say so rather than report a
        # clean switch and fail on the next request.
        if preset is None:
            result["warning"] = (
                f"{choice!r} is not a known preset, so it was taken as a raw "
                f"model id for the current provider ({new_config.provider}). "
                "Model ids are not validated here; if it is wrong, the next "
                "reply will fail with a provider error and you can switch back."
            )
        return result

    def use_config(self, config: Config) -> None:
        """Switch to another model/provider mid-session, keeping the conversation.

        Rebuilds the provider client (a preset switch can move from Mistral to
        Claude to a local server) but leaves ``self.history`` alone: history is
        stored in the provider-neutral OpenAI shape and translated per request,
        so a conversation survives switching provider partway through.

        The RAG index is untouched — it depends on the *embedding* model, which
        config deliberately holds steady across chat-model changes.
        """
        self.config = config
        self.client = make_client(config)
        self.chat_model = config.chat_model
        tools_module.set_session_model(config.provider, config.chat_model)

    def _complete(
        self,
        messages: list[dict[str, Any]],
        *,
        with_tools: bool,
        tool_choice: Any = "auto",
    ) -> Any:
        """Call the chat model, recovering from unparseable structured replies.

        mistral-large occasionally returns citations as structured reference
        chunks the SDK cannot parse; retrying with a plain-text nudge recovers
        (a no-op for OpenAI-compatible providers). Transient errors (rate limit,
        5xx, connection) are retried with backoff. ``tool_choice`` may be "auto",
        a provider "force some tool" token, or a named tool_choice dict.
        Provider-specific exception handling lives in the ``llm`` adapter.
        """
        attempt_messages = messages
        last_error: Exception | None = None
        for attempt in range(MAX_COMPLETION_ATTEMPTS):
            try:
                return self.client.complete(
                    messages=attempt_messages,
                    model=self.chat_model,
                    temperature=self.config.temperature,
                    tools=tools_module.TOOLS if with_tools else None,
                    tool_choice=tool_choice if with_tools else None,
                )
            except Exception as exc:  # noqa: BLE001 - classified via the adapter
                last_error = exc
                if self.client.is_response_validation_error(exc):
                    attempt_messages = [
                        *messages,
                        {
                            "role": "system",
                            "content": "Reply in plain text only. Do not use any "
                            "structured citation or reference format; write any "
                            "citations as ordinary parenthetical text.",
                        },
                    ]
                    continue
                if self.client.is_retryable(exc):
                    time.sleep(self.client.backoff_seconds(exc, attempt))
                    continue
                raise
        raise RuntimeError(
            f"Chat request failed after {MAX_COMPLETION_ATTEMPTS} attempts "
            f"(last error: {last_error})"
        )

    def _context_window(self) -> Optional[int]:
        return providers.context_tokens(self.chat_model)

    def _read_budget(self, system: str) -> Optional[int]:
        """Tokens a file read may still use this turn, or None when unlimited."""
        window = self._context_window()
        if window is None:
            return None
        used = _estimate_tokens(system, tools_module.TOOLS, self.history)
        return max(window - used - _ANSWER_RESERVE, 300)

    def _fitted_history(self, system: str) -> list[dict[str, Any]]:
        """The history, minus the oldest whole turns if it would not fit.

        Ollama drops the START of an over-long prompt — which is where the system
        prompt sits — so it is better to drop the oldest turns ourselves. A turn
        is cut whole (from one user message to the next), so a tool call never
        loses its result. The current turn is always kept.
        """
        window = self._context_window()
        if window is None:
            return self.history
        limit = window - _ANSWER_RESERVE
        fixed = _estimate_tokens(system, tools_module.TOOLS)
        history = self.history
        starts = [i for i, m in enumerate(history)
                  if m.get("role") == "user" and isinstance(m.get("content"), str)]
        dropped = 0
        while fixed + _estimate_tokens(history) > limit and len(starts) > 1:
            starts.pop(0)
            history = self.history[starts[0]:]
            dropped += 1
        if dropped:
            self._emit(f"dropped the {dropped} oldest turn(s) to fit the "
                       f"{window}-token context window")
        return history

    def _complete_pending_campaign(self, pending: dict[str, Any],
                                   plan: tuple[int, int], user_input: str) -> None:
        """Make the run_campaign call the user's CPU-plan answer completes."""
        args = {**pending, "n_procs_per_case": plan[0], "n_concurrent_cases": plan[1]}
        if _forbids_starting(user_input):
            args["launch"] = False
        if self._on_tool_call is not None:
            self._on_tool_call("run_campaign", args)
        result = tools_module.execute_tool("run_campaign", args)
        call_id = f"cpuplan-{len(self.history)}"
        self.history.append({
            "role": "assistant", "content": "",
            "tool_calls": [{"id": call_id, "type": "function", "function": {
                "name": "run_campaign", "arguments": json.dumps(args)}}],
        })
        self.history.append({"role": "tool", "tool_call_id": call_id,
                             "name": "run_campaign", "content": _tool_content(result)})

    def _inject_forced_read(self, hits: list[Hit], system: str = "",
                            question: str = "") -> None:
        """Deterministically read the most relevant source file into the history.

        Picks the highest-ranked C/C++ hit, reads it with read_case_file, and
        records a well-formed assistant tool-call + tool-result so the model
        answers from the actual source. Unlike a tool_choice-forced read this does
        NOT depend on the provider honouring tool_choice (local Ollama models treat
        it as advisory), so the source is guaranteed to be in context for any
        provider. Tools stay available, so the model can read further files itself.
        """
        cpp = next((hit for hit in hits if hit.chunk.kind == "cpp"), None)
        path = _source_named_by_notes(hits, question) or (cpp and cpp.chunk.rel_path)
        if not path:
            return
        result = None
        # A header (.H/.h) holds declarations; prefer its sibling implementation
        # (.C/.c), which has the actual code/equations, when that file exists.
        if path.endswith((".H", ".h")):
            impl = path[:-2] + (".C" if path.endswith(".H") else ".c")
            probe = tools_module.execute_tool("read_case_file", {"path": impl})
            if probe.get("status") == "ok":
                path, result = impl, probe
        if result is None:
            result = tools_module.execute_tool("read_case_file", {"path": path})
        budget = self._read_budget(system)
        if budget is not None:
            focus = [(h.chunk.start_line, h.chunk.end_line)
                     for h in hits if h.chunk.rel_path == path]
            result = _fit_file_read(result, budget, focus, question)
        args = {"path": path}
        self._emit("implementation question -> reading the top source file for grounding")
        if self._on_tool_call is not None:
            self._on_tool_call("read_case_file", args)
        call_id = f"forcedread-{len(self.history)}"
        self.history.append({
            "role": "assistant",
            "content": "",
            "tool_calls": [{
                "id": call_id,
                "type": "function",
                "function": {"name": "read_case_file", "arguments": json.dumps(args)},
            }],
        })
        self.history.append({
            "role": "tool",
            "tool_call_id": call_id,
            "name": "read_case_file",
            "content": _tool_content(result),
        })

    def ask(self, user_input: str) -> str:
        """Process one user message and return the assistant's reply text.

        On any failure (e.g. an exhausted rate-limit retry) the whole turn is
        rolled back, so the history never keeps a dangling user message or a
        tool-call without its results and the next turn stays valid.
        """
        snapshot = len(self.history)
        self.history.append({"role": "user", "content": user_input})
        try:
            return self._run_turn(user_input)
        except Exception:
            del self.history[snapshot:]
            raise

    def _retrieval_query(self, user_input: str) -> str:
        """Build the index query from the previous user turn + the current message.

        A context-dependent follow-up ("I don't understand", "why?") carries no
        topic of its own, so embedding it alone retrieves noise. Prepending the
        previous user turn lets such follow-ups inherit the topic and surface the
        right files. Only the immediately preceding user turn is used, to keep the
        query focused (more history dilutes the embedding). ``self.history`` already
        holds the current message as its last entry, so it is excluded here.
        """
        prior_user = [
            m["content"]
            for m in self.history[:-1]
            if m.get("role") == "user" and isinstance(m.get("content"), str)
        ]
        if not prior_user:
            return user_input
        return prior_user[-1] + "\n" + user_input

    def _run_turn(self, user_input: str) -> str:
        hits = retrieve(self.config, self._retrieval_query(user_input), self.index)
        if hits:
            top = ", ".join(hit.chunk.citation for hit in hits[:3])
            self._emit(f"retrieved {len(hits)} context chunks (top: {top})")
        else:
            self._emit("retrieved no matching context")
        system = SYSTEM_PROMPT + "\n\n# Retrieved context\n" + format_context(hits)

        # For "how does X work / how is it implemented" questions that surface a
        # source file, force the model to read the implementation on the first
        # step instead of answering from the partial retrieved excerpts. The
        # keyword intent is judged on the CURRENT message only, so a topic switch
        # ("now set up a case ...") after an implementation question is never
        # falsely forced — even though retrieval above uses the previous turn too.
        force_read = _wants_implementation(user_input) and (
            any(hit.chunk.kind == "cpp" for hit in hits)
            or _source_named_by_notes(hits, user_input) is not None
        )
        if force_read:
            # Deterministically read the most relevant source file into the
            # conversation BEFORE the model answers, so the answer is grounded in
            # the real source regardless of whether the provider enforces
            # tool_choice (local Ollama models do not).
            self._inject_forced_read(hits, system, user_input)

        # The user answered a pending campaign's CPU-plan question: complete it.
        pending = tools_module.pending_campaign()
        plan = _cpu_plan(user_input) if pending else None
        if plan and not _is_feasibility_question(user_input) and self._guard(
                "repair", "cpu_plan_completion",
                f"CPU plan answered ({plan[0]} per case, {plan[1]} at a time) "
                "-> completing the pending campaign"):
            self._complete_pending_campaign(pending, plan, user_input)

        # Tracks whether a case was prepared (built/changed) this turn, to enforce
        # the "separate steps" guard below (no auto-run in the same turn).
        prepared_this_turn = False
        # A feasibility/possibility QUESTION must be answered, not executed — block
        # every action tool this turn so the model cannot build/run/delete (and
        # fabricate parameters) in response to "is it possible to ...?".
        feasibility_question = _is_feasibility_question(user_input)
        # "Prepare it but do not start it": no run tool this turn, and a campaign
        # is built with launch=false whatever the model passes.
        dont_start = _forbids_starting(user_input)
        # Everything the user has said so far, to confirm a case type was asked for.
        said = "\n".join(m["content"] for m in self.history
                         if m.get("role") == "user" and isinstance(m.get("content"), str))
        # (tool, status, result) of every call this turn, for the claim check.
        outcomes: list[tuple[str, str, dict[str, Any]]] = []
        # (tool, args) of every executed call, for the requested-value check.
        executed: list[tuple[str, dict[str, Any]]] = []
        nudge: Optional[dict[str, Any]] = None
        # A nudge guard already reported this reply in shadow mode, so the
        # correction of the same claim is not counted a second time.
        flagged = False

        for step in range(MAX_TOOL_ROUNDS + 1):
            messages = [{"role": "system", "content": system},
                        *self._fitted_history(system)]
            response = self._complete(messages, with_tools=True)
            msg = response.choices[0].message

            if not msg.tool_calls:
                text = _message_text(msg)
                missing = _unapplied_values(user_input, executed)
                claimed = _claim_correction(user_input, text, outcomes) or (
                    not executed and _PROMISED_CALL.search(text) is not None)
                if claimed and not missing and nudge is None and step < MAX_TOOL_ROUNDS:
                    # The reply claims an action no tool performed: one chance to
                    # actually perform it before the correction is appended.
                    if self._guard("nudge", "false_claim",
                                   "reply claimed an action no tool performed — "
                                   "asking the model to do it or retract"):
                        nudge = {"role": "system", "content": (
                            "Check before answering: your draft says something was "
                            "(or will be) built, changed, prepared or started, but no "
                            "tool did that in this turn. If the user asked for it, call "
                            "the tool now; otherwise answer without claiming it."
                        )}
                        self.history.append(nudge)
                        continue
                    flagged = True
                if missing and nudge is None and step < MAX_TOOL_ROUNDS:
                    # Give the model one chance to make the call it missed; the
                    # draft reply is discarded and the note is not kept afterwards.
                    if self._guard("nudge", "unapplied_value",
                                   "requested value(s) not used by any tool: "
                                   + ", ".join(missing)
                                   + " — asking the model to apply them"):
                        nudge = {"role": "system", "content": (
                            "Check before answering: the user asked for "
                            + ", ".join(missing) + ", but no tool call in this turn "
                            "used that value. Call the right tool with it now (lengths "
                            "in metres: 5 mm = 0.005; conditions or geometry of an "
                            "existing case: build_case with its case_name), or say "
                            "plainly why it was not applied."
                        )}
                        self.history.append(nudge)
                        continue
                    flagged = True
                if nudge is not None:
                    self.history = [m for m in self.history if m is not nudge]
                correction = _claim_correction(user_input, text, outcomes)
                if missing and nudge is not None:
                    # Nudged and still not applied: say so rather than let it slip.
                    correction = (correction or "") + (
                        "\n\n> ⚠️ **Not applied:** " + ", ".join(missing)
                        + " was not passed to any tool in this step, so it has not "
                        "been applied. Please ask again, naming the parameter."
                    )
                if correction and not flagged and self._guard(
                        "nudge", "claim_correction",
                        "reply claimed an action no tool performed — corrected"):
                    text += correction
                self.history.append({"role": "assistant", "content": text})
                return text

            # Record the assistant's tool request, then run each tool.
            self.history.append({
                "role": "assistant",
                "content": _message_text(msg),
                "tool_calls": [
                    {
                        "id": tc.id,
                        "type": "function",
                        "function": {
                            "name": tc.function.name,
                            "arguments": tc.function.arguments,
                        },
                    }
                    for tc in msg.tool_calls
                ],
            })
            for tc in msg.tool_calls:
                name = tc.function.name
                try:
                    args = json.loads(tc.function.arguments or "{}")
                except json.JSONDecodeError:
                    args = {}
                if (name in _CASE_TYPE_TOOLS and not args.get("case_type")
                        and (named := _case_type_from_user(user_input, said))
                        and self._guard("repair", "case_type_fill",
                                        f"case_type={named} taken from your message")):
                    # The user named the case; the model left it out.
                    args["case_type"] = named
                if (name == "run_campaign" and not args.get("variable")
                        and args.get("case_type")
                        and (variable := _variable_from_user(str(args["case_type"]),
                                                             user_input))
                        and self._guard("repair", "variable_fill",
                                        f"variable={variable} taken from your message")):
                    args["variable"] = variable
                if name in _STATED_TOOLS | {"run_campaign"}:
                    swept = validator.canonical_name(str(args.get("variable") or ""))
                    for key, value in _stated_parameters(user_input).items():
                        if key == swept:
                            continue
                        current = _numbers(args.get(key))
                        if current and abs(current[0] - value) <= 0.01 * abs(value):
                            continue
                        if self._guard("repair", "stated_value",
                                       f"{key}={value:g} taken from your message"
                                       + (f" (the model passed {current[0]:g})"
                                          if current else "")):
                            args[key] = value
                if (name in _ACTION_TOOLS and feasibility_question
                        and self._guard("block", "feasibility_question",
                                        f"skipped {name} (feasibility question — "
                                        "answering, not acting)")):
                    # Feasibility-question guard: the user is asking IF something is
                    # possible/supported, so never act — defer and let the model
                    # answer from the retrieved constraints.
                    result = _DEFERRED_FEASIBILITY_RESULT
                elif (name in _RUN_TOOLS and dont_start
                      and self._guard("block", "dont_start",
                                      f"skipped {name} (the user said not to start it)")):
                    result = _DEFERRED_DONT_START_RESULT
                elif (name in _CASE_TYPE_TOOLS
                      and (template := _requested_template(name, args)) is not None
                      and not template.mentioned_in(said)
                      and case_templates.infer_template(
                          validator.canonical_name(k) for k in args) is not template
                      and self._guard("block", "case_type_not_named",
                                      f"skipped {name} (case type {template.name} was "
                                      "not named by the user — asking)")):
                    # The model chose a case type the user never asked for.
                    result = tools_module.missing_case_type_result()
                elif ((invented := _invented_required(name, args, said))
                      and self._guard("block", "invented_values",
                                      f"skipped {name} (values the user did not give: "
                                      + ", ".join(invented) + ")")):
                    # The model filled in values the user never gave.
                    result = {
                        "status": "missing_parameters",
                        "missing": invented,
                        "message": (
                            "Not built: the user did not give "
                            + ", ".join(invented) + ". Ask the user for them — "
                            "never invent required values."
                        ),
                    }
                elif (name in _RUN_TOOLS and prepared_this_turn
                      and self._guard("block", "separate_steps",
                                      f"skipped {name} (run is a separate step "
                                      "after build/change)")):
                    # Separate-steps guard: do NOT run in the same turn as a build/
                    # change. Defer instead of executing (the solver never launches
                    # unprompted), and let the model report "ready — run it next".
                    result = _DEFERRED_RUN_RESULT
                else:
                    if self._on_tool_call is not None:
                        self._on_tool_call(name, args)
                    if name in tools_module.SESSION_TOOLS:
                        # Mutates this session (the provider client), so it is
                        # handled here rather than by the stateless tool layer.
                        # Any remaining steps of this turn use the new model.
                        result = self.switch_model(str(args.get("model") or ""))
                    else:
                        # Converted on a copy, so shadow mode leaves args untouched.
                        converted = dict(args)
                        notes = _convert_stated_units(converted, said)
                        if notes and self._guard(
                                "repair", "unit_conversion",
                                "converted " + "; ".join(notes) + " (the unit you wrote)"):
                            args = converted
                        if (name == "run_campaign" and dont_start
                                and args.get("launch") is not False
                                and self._guard("repair", "dont_start_launch",
                                                "run_campaign forced to launch=false "
                                                "(the user said not to start it)")):
                            args["launch"] = False
                        result = tools_module.execute_tool(name, args)
                        executed.append((name, args))
                        if name == "read_case_file":
                            budget = self._read_budget(system)
                            if budget is not None:
                                result = _fit_file_read(result, budget, [], user_input)
                    if name in _PREPARE_TOOLS:
                        prepared_this_turn = True
                outcomes.append((name, str(result.get("status") or
                                           ("error" if "error" in result else "")), result))
                self.history.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "name": tc.function.name,
                    "content": _tool_content(result),
                })

        # Tool-call budget exhausted; ask the model for a plain summary.
        if nudge is not None:
            self.history = [m for m in self.history if m is not nudge]
        self.history.append({
            "role": "user",
            "content": "Please summarise the result so far in plain text.",
        })
        final = self._complete(
            [{"role": "system", "content": system}, *self._fitted_history(system)],
            with_tools=False,
        )
        text = _message_text(final.choices[0].message)
        self.history.append({"role": "assistant", "content": text})
        return text
