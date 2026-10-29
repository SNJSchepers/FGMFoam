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

"""Tests for the guards that keep small models on track — each one added after
the eval (agent/evals/) caught a local model doing exactly this. Tool-layer
guards are called directly; chat-loop guards run through ChatAgent with a
scripted fake model. Offline: no model, no GPU, no solver.

    agent/.venv/bin/python agent/tests/test_guards.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import campaign as campaign_module  # noqa: E402
import chat  # noqa: E402
import config as config_module  # noqa: E402
import of_writer  # noqa: E402
import providers  # noqa: E402
import tools  # noqa: E402

# Everything is written to a scratch directory, and nothing can launch a solver.
SCRATCH = Path(tempfile.mkdtemp(prefix="fgm-guards-"))
of_writer.TESTCASES_DIR = SCRATCH / "testcases"
of_writer.TESTCASES_DIR.mkdir()
campaign_module.CAMPAIGN_DIR = SCRATCH / "campaigns"
LAUNCHED: list[str] = []


def _launch_stub(name):
    def run(args):
        LAUNCHED.append(name)
        return {"status": "started"}
    return run


for _name in ("run_case", "run_cases", "restart_case", "start_campaign"):
    tools._DISPATCH[_name] = _launch_stub(_name)
    setattr(tools, f"_run_{_name}", tools._DISPATCH[_name])

SLOT = {"case_type": "laminarSlotBurner", "inlet_velocity": 0.4,
        "T_inlet": 300, "T_plate": 400}


# --------------------------------------------------------------------------- #
# Tool layer                                                                  #
# --------------------------------------------------------------------------- #

def test_lengths_given_in_mm_are_refused() -> None:
    tools.reset_session_cases()
    result = tools.execute_tool("build_case", {**SLOT, "plate_thickness": 5})
    assert result["status"] == "error" and "METRES" in result["message"]
    result = tools.execute_tool("run_campaign", {
        **SLOT, "variable": "plate_thickness", "values": [2, 3, 4],
        "n_procs_per_case": 1, "n_concurrent_cases": 1, "launch": False})
    assert result["status"] == "error" and "plate_thickness" in result["parameters"]


def test_a_rebuild_keeps_what_the_user_did_not_change() -> None:
    tools.reset_session_cases()
    first = tools.execute_tool("build_case", {**SLOT, "case_name": "geo"})
    assert first["status"] == "built"
    second = tools.execute_tool("build_case", {"case_type": "laminarSlotBurner",
                                               "case_name": "geo",
                                               "plate_thickness": 0.005})
    assert second["status"] == "built", second
    assert second["parameters"]["T_plate"] == 400
    assert second["parameters"]["plate_thickness"] == 0.005
    assert second["inherited_from"]["case_name"] == "laminarSlotBurner_geo"
    record = json.loads((Path(second["case_dir"]) / "agent-provenance.json").read_text())
    assert record["parameter_sources"]["T_plate"] == "inherited"
    assert record["parameter_sources"]["plate_thickness"] == "user"
    # The changed value now carries on to the next rebuild too.
    third = tools.execute_tool("build_case", {"case_type": "laminarSlotBurner",
                                              "hole_size": 0.006})
    assert third["parameters"]["plate_thickness"] == 0.005


def test_change_settings_rebuilds_for_geometry_instead_of_ignoring_it() -> None:
    tools.reset_session_cases()
    tools.execute_tool("build_case", {**SLOT, "case_name": "cs",
                                      "plate_thickness": 0.005})
    tools.reset_session_cases()   # the base must come from the case itself
    result = tools.execute_tool("change_settings", {"case_name": "laminarSlotBurner_cs",
                                                    "hole_size": 0.006})
    assert result["status"] == "built", result
    assert result["parameters"]["hole_size"] == 0.006
    assert result["parameters"]["plate_thickness"] == 0.005
    assert result["parameters"]["T_plate"] == 400
    assert "rebuilt" in result["note"]
    # Run controls alone are still edited in place.
    edited = tools.execute_tool("change_settings", {"case_name": "laminarSlotBurner_cs",
                                                    "sim_time": 0.02})
    assert edited["status"] == "updated"


def test_nothing_is_inherited_across_case_types_or_conversations() -> None:
    tools.reset_session_cases()
    tools.execute_tool("build_case", {**SLOT, "case_name": "a"})
    flamelet = tools.execute_tool("build_case", {"case_type": "laminarFlamelet",
                                                 "case_name": "b", "T_inlet": 350})
    assert "inherited_from" not in flamelet
    tools.reset_session_cases()
    fresh = tools.execute_tool("build_case", {"case_type": "laminarSlotBurner",
                                              "inlet_velocity": 0.5})
    assert fresh["status"] == "missing_parameters"


def test_case_tools_default_to_the_last_built_case() -> None:
    tools.reset_session_cases()
    built = tools.execute_tool("build_case", {**SLOT, "case_name": "last"})
    checked = tools.execute_tool("check_case", {})
    assert checked.get("case_name") == built["case_name"], checked


def test_no_case_named_and_none_built_asks_instead_of_guessing() -> None:
    tools.reset_session_cases()
    shutil.rmtree(of_writer.case_path(of_writer.DEFAULT_CASE_NAME), ignore_errors=True)
    result = tools.execute_tool("check_case", {})
    assert result["status"] == "error" and "which case" in result["error"]


def test_case_type_is_read_from_the_case_name_prefix() -> None:
    tools.reset_session_cases()
    result = tools.execute_tool("build_case", {"case_name": "laminarFlamelet_cond",
                                               "T_inlet": 300})
    assert result["status"] == "built" and result["case_type"] == "laminarFlamelet"


def test_a_campaign_never_defaults_its_case_type() -> None:
    result = tools.execute_tool("run_campaign", {"variable": "T_inlet",
                                                 "values": [300, 400]})
    assert result["status"] == "missing_case_type"


def test_the_tool_never_guesses_the_swept_variable() -> None:
    # Once guessed from "the flamelet has one sweepable parameter", an
    # equivalence-ratio sweep became T_inlet = 0.6-0.8 K.
    result = tools.execute_tool("run_campaign", {"case_type": "laminarFlamelet",
                                                 "values": [0.6, 0.7, 0.8]})
    assert result["status"] == "error"


def test_temperatures_outside_kelvin_range_are_refused() -> None:
    tools.reset_session_cases()
    result = tools.execute_tool("build_case", {"case_type": "laminarFlamelet",
                                               "T_inlet": 0.7})
    assert result["status"] == "error" and "KELVIN" in result["message"]
    result = tools.execute_tool("run_campaign", {
        "case_type": "laminarFlamelet", "variable": "T_inlet", "values": [0.6, 0.7],
        "n_procs_per_case": 1, "n_concurrent_cases": 1, "launch": False})
    assert result["status"] == "error" and "T_inlet" in result["parameters"]


def test_campaign_cases_do_not_inherit_from_the_last_build() -> None:
    tools.reset_session_cases()
    tools.execute_tool("build_case", {**SLOT, "case_name": "base",
                                      "plate_thickness": 0.005})
    result = tools.execute_tool("run_campaign", {
        **SLOT, "variable": "T_plate", "values": [350, 450], "campaign_name": "noinh",
        "n_procs_per_case": 1, "n_concurrent_cases": 1, "launch": False})
    assert result["status"] == "prepared", result
    record = json.loads((of_writer.case_path(result["cases"][0]["case"]) /
                         "agent-provenance.json").read_text())
    assert record["parameters"]["plate_thickness"] != 0.005


def test_context_window_is_read_from_the_local_tag() -> None:
    assert providers.context_tokens("qwen2.5-3b-16k") == 16384
    assert providers.context_tokens("claude-sonnet-5") is None


# --------------------------------------------------------------------------- #
# Chat loop, with a scripted model                                            #
# --------------------------------------------------------------------------- #

def _call(name, args, call_id="c1"):
    return SimpleNamespace(id=call_id, function=SimpleNamespace(
        name=name, arguments=json.dumps(args)))


class ScriptedClient:
    """Returns the next scripted reply: a list of tool calls, or a final text."""

    def __init__(self, script):
        self.script = list(script)

    def complete(self, **kwargs):
        step = self.script.pop(0)
        if isinstance(step, str):
            message = SimpleNamespace(content=step, tool_calls=None)
        else:
            message = SimpleNamespace(content="", tool_calls=step)
        return SimpleNamespace(choices=[SimpleNamespace(message=message)])


def _agent(script, model="qwen2.5-3b-16k", guards="on"):
    chat.retrieve = lambda *a, **k: []
    cfg = config_module.load_config(require_key=False, model="local-3b", guards=guards)
    agent = chat.ChatAgent(cfg, index=None)
    agent.client = ScriptedClient(script)
    agent.chat_model = model
    agent.events = []
    agent._on_event = agent.events.append
    # (kind, name, acted) of every guard that detected something.
    agent.guards = []
    agent._on_guard = lambda kind, name, message, acted: agent.guards.append(
        (kind, name, acted))
    return agent


def test_do_not_start_blocks_runs_and_forces_launch_false() -> None:
    LAUNCHED.clear()
    agent = _agent([
        [_call("run_campaign", {"case_type": "laminarFlamelet", "variable": "T_inlet",
                                "values": [300, 400], "campaign_name": "nostart",
                                "n_procs_per_case": 1, "n_concurrent_cases": 1})],
        [_call("run_cases", {"case_names": ["x"]}, "c2")],
        "Prepared.",
    ])
    agent.ask("Run a flamelet campaign at 300 and 400 K, but do not start it yet.")
    assert LAUNCHED == []
    results = [json.loads(m["content"]) for m in agent.history if m["role"] == "tool"]
    assert results[0]["status"] == "prepared"
    assert results[1]["status"] == "not_started"


def test_an_invented_case_type_is_sent_back_as_a_question() -> None:
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", "T_inlet": 350})],
        "Which case?",
    ])
    agent.ask("Set up a simulation with an unburnt temperature of 350 K.")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "missing_case_type"


def test_a_case_type_the_user_named_or_implied_is_accepted() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", "inlet_velocity": 0.4,
                              "T_inlet": 300, "T_plate": 400})],
        "Built.",
    ])
    # No case name in the message, but inlet velocity only exists for the burner.
    agent.ask("Set up a case at 0.4 m/s inlet velocity, 300 K inlet, 400 K plate.")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "built"


def test_a_claimed_build_that_did_not_happen_is_corrected() -> None:
    agent = _agent([
        [_call("compute_thermochemistry", {"equivalence_ratio": 0.7, "T_inlet": 450})],
        "The case has been rebuilt successfully at 450 K.",
        # Nudged once to actually do it; repeating the claim gets the correction.
        "The case has been rebuilt successfully at 450 K.",
    ])
    reply = agent.ask("Rebuild the flamelet case at 450 K.")
    assert "Correction from the agent" in reply


def _note_hit(heading: str) -> object:
    """The retrieved chunk of agent/knowledge/solver_notes.md under ``heading``."""
    from rag.chunker import Chunk
    from rag.index import Hit
    lines = (config_module.REPO_ROOT / "agent/knowledge/solver_notes.md").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("## " + heading))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    return Hit(Chunk("agent/knowledge/solver_notes.md", start + 1, end, "markdown",
                     "\n".join(lines[start:end])), 0.5)


def test_the_note_whose_heading_matches_steers_the_forced_read() -> None:
    # As retrieved for this question: the boundary-condition note ranked first.
    hits = [_note_hit("Enthalpy boundary condition"), _note_hit("Table lookup")]
    question = "Okay how does the table lookup work, we have Yc bounds which vary with enthalpy?"
    assert chat._source_named_by_notes(hits, question) == "solver/src/lookUp/FGM/FGMlib.c"
    # A question about the boundary condition still opens the boundary condition.
    bc = "How does the lookUpEnthalpy boundary condition work?"
    assert chat._source_named_by_notes(hits[::-1], bc).endswith(
        "lookUpEnthalpyFvPatchScalarField.C")


def test_a_word_in_most_paths_earns_no_path_bonus() -> None:
    # "solver" is in every solver/ path: it must not lift all source code above
    # the knowledge note that answers "what is the solver algorithm?".
    from rag.chunker import Chunk
    from rag.index import _lexical_scores
    chunks = [Chunk(f"solver/src/lib{i}/file{i}.C", 1, 2, "cpp", "code") for i in range(6)]
    chunks.append(Chunk("agent/knowledge/notes.md", 1, 2, "markdown", "algorithm"))
    scores = _lexical_scores("What is the solver algorithm", chunks)
    assert scores[-1] > scores[0]
    # A word that does name a file still earns it.
    assert _lexical_scores("file3", chunks)[3] > _lexical_scores("file3", chunks)[0]


def test_model_is_not_a_topic_word() -> None:
    # "model" is a substring of every combustionModel/ path; as a query word it
    # put that directory's licence headers above the diffusion-model note.
    from rag.index import _lexical_query_terms
    assert _lexical_query_terms("Which diffusion model is used?") == ["diffusion"]


def test_a_note_may_name_its_file_with_line_numbers() -> None:
    text = "loop (`solver/applications/FGMFoam/FGMFoam.C:145-209`) does, in order:"
    assert chat._SOURCE_PATH.findall(text) == ["solver/applications/FGMFoam/FGMFoam.C"]


def test_promises_and_denials_are_not_claims() -> None:
    # Mistral Large 4 phrasings that once tripped the claim check.
    request = "Set up a slot burner case at 0.4 m/s."
    for reply in ("The campaign is ready to be built, but I need your CPU plan.",
                  "The campaign is prepared and saved; nothing has been started.",
                  "Tell me the plate temperature; nothing is running yet.",
                  "I'd rather stop you before anything is built."):
        assert chat._claim_correction(request, reply, []) is None, reply
    assert chat._claim_correction(request, "The case has been built.", [])


def test_a_false_claim_first_gets_a_chance_to_do_the_work() -> None:
    tools.reset_session_cases()
    agent = _agent([
        "The flamelet case has been built at 300 K.",               # claims, did nothing
        [_call("build_case", {"case_type": "laminarFlamelet", "T_inlet": 300})],
        "The flamelet case has been built at 300 K.",               # now true
    ])
    reply = agent.ask("Build a flamelet case at 300 K.")
    assert "Correction" not in reply
    assert not any(m.get("role") == "system" for m in agent.history)


def test_a_requested_value_no_tool_used_is_sent_back_once() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_type": "laminarFlamelet", "case_name": "v",
                              "T_inlet": 300})],
        "Done.",
        [_call("change_settings", {"case_name": "laminarFlamelet_v", "sim_time": 0.1},
               "c2")],
        "Updated as requested.",                                     # 450 K dropped
        [_call("build_case", {"case_name": "laminarFlamelet_v", "T_inlet": 450}, "c3")],
        "Rebuilt at 450 K.",
    ])
    agent.ask("Set up a flamelet case named v at 300 K.")
    reply = agent.ask("Rebuild it at 450 K instead.")
    assert any("not used by any tool" in e for e in agent.events)
    built = [json.loads(m["content"]) for m in agent.history if m["role"] == "tool"][-1]
    assert built["status"] == "built" and built["parameters"]["T_inlet"] == 450
    assert reply == "Rebuilt at 450 K."


def test_required_values_the_user_never_gave_are_not_accepted() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", "inlet_velocity": 0.5})],
        [_call("build_case", {"case_type": "laminarSlotBurner", "inlet_velocity": 0.5,
                              "T_inlet": 300, "T_plate": 300}, "c2")],
        "Which temperatures?",
    ])
    agent.ask("Set up a slot burner case with an inlet velocity of 0.5 m/s.")
    results = [json.loads(m["content"]) for m in agent.history if m["role"] == "tool"]
    assert results[1]["status"] == "missing_parameters"
    assert set(results[1]["missing"]) == {"T_inlet", "T_plate"}


def test_values_the_user_gave_in_other_units_are_accepted() -> None:
    from chat import _invented_required
    said = "slot burner, 0.4 m/s, inlet 300 K, plate 400 K, plate thickness 5 mm"
    assert _invented_required("build_case", {
        "case_type": "laminarSlotBurner", "inlet_velocity": 0.4, "T_inlet": 300,
        "T_plate": 400, "plate_thickness": 0.005}, said) == []
    assert _invented_required("build_case", {
        "case_type": "laminarSlotBurner", "inlet_velocity": 0.4, "T_inlet": 300,
        "T_plate": 350}, said) == ["T_plate"]


def test_lengths_in_the_unit_the_user_wrote_are_converted() -> None:
    from chat import _convert_stated_units
    args = {"case_name": "x", "plate_thickness": 5, "hole_size": 0.006}
    notes = _convert_stated_units(args, "rebuild with a plate thickness of 5 mm")
    assert args["plate_thickness"] == 0.005 and args["hole_size"] == 0.006 and notes
    sweep = {"variable": "plate_thickness", "values": [2, 3, 4]}
    _convert_stated_units(sweep, "sweep the plate thickness over 2, 3 and 4 mm")
    assert sweep["values"] == [0.002, 0.003, 0.004]
    # No unit in the user's text: left alone (and refused by the tool).
    bare = {"plate_thickness": 5}
    assert not _convert_stated_units(bare, "plate thickness 5") and bare["plate_thickness"] == 5


def test_a_case_type_the_user_named_is_filled_in() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_name": "cond", "T_inlet": 300})],
        "Built.",
    ])
    agent.ask("Set up a laminar flamelet case named cond at 300 K.")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "built" and result["case_type"] == "laminarFlamelet"


def test_an_announced_call_that_was_not_made_is_sent_back() -> None:
    tools.reset_session_cases()
    tools.execute_tool("build_case", {"case_type": "laminarFlamelet", "case_name": "p",
                                      "T_inlet": 300})
    agent = _agent([
        "I will call the `change_settings` function with sim_time 0.02.",
        [_call("change_settings", {"case_name": "laminarFlamelet_p", "sim_time": 0.02})],
        "Done: end time is now 0.02 s.",
    ])
    reply = agent.ask("Change the end time of laminarFlamelet_p to 0.02 s.")
    assert reply.startswith("Done")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "updated"


def test_a_list_for_one_parameter_is_refused_not_dropped() -> None:
    tools.reset_session_cases()
    result = tools.execute_tool("build_case", {**SLOT, "plate_thickness": [0.002, 0.003]})
    assert result["status"] == "error" and "run_campaign" in result["message"]


def test_the_swept_variable_is_read_from_the_message() -> None:
    from chat import _variable_from_user
    assert _variable_from_user("laminarSlotBurner",
                               "sweep the plate thickness over 2, 3 and 4 mm") == "plate_thickness"
    assert _variable_from_user("laminarSlotBurner", "sweep it over 2, 3 and 4") is None
    # Other parameters named with single values are the fixed conditions.
    assert _variable_from_user(
        "laminarSlotBurner", "a slot burner campaign that sweeps the plate thickness "
        "over 2, 3 and 4 mm, at inlet velocity 0.4 m/s, inlet temperature 300 K and "
        "plate temperature 400 K") == "plate_thickness"
    assert _variable_from_user("laminarFlamelet",
                               "flamelets at unburnt temperatures of 300, 400 and 500 K") == "T_inlet"
    assert _variable_from_user("laminarFlamelet",
                               "a flamelet campaign over equivalence ratios 0.6, 0.7 and 0.8") is None


def test_stated_values_are_read_with_their_units() -> None:
    from chat import _stated_parameters as stated
    assert stated("Now rebuild that case with a Plate Thickness of 5 mm.") == {
        "plate_thickness": 0.005}
    assert stated("make the slot (hole size) 6 mm wide") == {"hole_size": 0.006}
    assert stated("an inlet velocity of 0.4 m/s, an inlet temperature of 300 K and "
                  "a plate temperature of 400 K") == {
        "inlet_velocity": 0.4, "T_inlet": 300.0, "T_plate": 400.0}
    # A unit that does not fit the parameter never matches it.
    assert "plate_thickness" not in stated("a plate temperature of 400 K")
    # Lists are campaigns, not one stated value.
    assert stated("unburnt temperatures of 300, 400 and 500 K") == {}


def test_a_stated_value_reaches_the_tool_the_model_chose() -> None:
    tools.reset_session_cases()
    tools.execute_tool("build_case", {**SLOT, "case_name": "st"})
    agent = _agent([
        # The 3B's pattern: change_settings with only the old run controls.
        [_call("change_settings", {"case_name": "laminarSlotBurner_st", "sim_time": 0.1})],
        "Updated.",
    ])
    agent.ask("Now rebuild the slot burner with a plate thickness of 5 mm.")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "built", result
    assert result["parameters"]["plate_thickness"] == 0.005
    assert result["parameters"]["T_plate"] == 400


def test_a_cpu_plan_answer_completes_the_pending_campaign() -> None:
    tools.reset_session_cases()
    LAUNCHED.clear()
    agent = _agent([
        [_call("run_campaign", {"case_type": "laminarFlamelet", "variable": "T_inlet",
                                "values": [300, 400], "campaign_name": "pend"})],
        "How should I split the cores?",
        "Prepared, not started.",          # the 3B's pattern: no call of its own
    ])
    agent.ask("Run a campaign of laminar flamelets at 300 and 400 K.")
    agent.ask("Use 1 processor per case and run 1 case at the same time. "
              "Only prepare it, do not start it yet.")
    results = [json.loads(m["content"]) for m in agent.history if m["role"] == "tool"]
    assert results[0]["status"] == "missing_cpu_plan"
    assert results[1]["status"] == "prepared", results[1]
    assert LAUNCHED == [] and tools.pending_campaign() == {}


def test_stated_conditions_reach_a_campaign_but_not_the_swept_variable() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("run_campaign", {"case_type": "laminarSlotBurner",
                                "variable": "plate_thickness", "values": [2, 3, 4],
                                "variable_values": [{"T_inlet": 300}],
                                "n_procs_per_case": 1, "n_concurrent_cases": 1})],
        "Prepared.",
    ])
    agent.ask("Prepare, but do not start, a slot burner campaign that sweeps the plate "
              "thickness over 2, 3 and 4 mm, at inlet velocity 0.4 m/s, inlet "
              "temperature 300 K and plate temperature 400 K.")
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "prepared", result
    assert sorted(result["values"]) == [0.002, 0.003, 0.004]


def test_an_excerpt_without_ranges_shows_the_part_the_question_is_about() -> None:
    read = tools.execute_tool("read_case_file", {
        "path": "solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C"})
    fitted = chat._fit_file_read(read, 2000, [], "what does the update step write "
                                 "into the flow fields, and negative source terms")
    assert "thermo:psi" in fitted["content"]


def test_old_turns_are_dropped_before_the_prompt_overflows() -> None:
    agent = _agent(["ok"], model="tiny-2k")
    agent.history = [{"role": "user", "content": "x" * 20000},
                     {"role": "assistant", "content": "y"}]
    agent.ask("hello")
    assert any("dropped" in e for e in agent.events)


# --------------------------------------------------------------------------- #
# Shadow mode (FGM_GUARDS=shadow): detect and report, never act               #
# --------------------------------------------------------------------------- #

def test_guards_report_that_they_acted() -> None:
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", "T_inlet": 350})],
        "Which case?",
    ])
    agent.ask("Set up a simulation with an unburnt temperature of 350 K.")
    assert agent.guards == [("block", "case_type_not_named", True)]


def test_shadow_mode_lets_a_blocked_call_through() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", "T_inlet": 350})],
        "Built.",
    ], guards="shadow")
    agent.ask("Set up a simulation with an unburnt temperature of 350 K.")
    assert agent.guards == [("block", "case_type_not_named", False)]
    assert any(e.startswith("[shadow, not applied] skipped build_case") for e in agent.events)
    # The tool itself ran, and its own validation still answers.
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] != "missing_case_type"


def test_shadow_mode_does_not_enforce_do_not_start() -> None:
    LAUNCHED.clear()
    agent = _agent([
        [_call("run_campaign", {"case_type": "laminarFlamelet", "variable": "T_inlet",
                                "values": [300, 400], "campaign_name": "shadowstart",
                                "n_procs_per_case": 1, "n_concurrent_cases": 1})],
        [_call("run_cases", {"case_names": ["x"]}, "c2")],
        "Prepared.",
    ], guards="shadow")
    agent.ask("Run a flamelet campaign at 300 and 400 K, but do not start it yet.")
    assert ("repair", "dont_start_launch", False) in agent.guards
    assert ("block", "dont_start", False) in agent.guards
    assert "run_cases" in LAUNCHED            # the (stubbed) launch really happened


def test_shadow_mode_leaves_the_arguments_alone() -> None:
    tools.reset_session_cases()
    agent = _agent([
        [_call("build_case", {"case_type": "laminarSlotBurner", **SLOT,
                              "plate_thickness": 5})],
        "Built.",
    ], guards="shadow")
    agent.ask("Build a slot burner at 0.4 m/s, 300 K inlet, 400 K plate, "
              "with a plate thickness of 5 mm.")
    kinds = {name for _, name, _ in agent.guards}
    assert {"stated_value"} <= kinds and all(not acted for *_, acted in agent.guards)
    # Neither the stated value nor the unit conversion was applied: the tool
    # refuses 5 (metres) exactly as it would without the agent around it.
    result = json.loads(next(m for m in agent.history if m["role"] == "tool")["content"])
    assert result["status"] == "error" and "METRES" in result["message"]


def test_shadow_mode_reports_a_false_claim_once_and_does_not_retry() -> None:
    agent = _agent(["The flamelet case has been built at 300 K."], guards="shadow")
    reply = agent.ask("Build a flamelet case at 300 K.")
    assert agent.guards == [("nudge", "false_claim", False)]
    assert "Correction" not in reply
    assert not any(m.get("role") == "system" for m in agent.history)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    try:
        for test in tests:
            test()
            print(f"  ok  {test.__name__}")
    finally:
        shutil.rmtree(SCRATCH, ignore_errors=True)
    print(f"\n{len(tests)} passed")
