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

"""Tests for the case plugin system: discovery, validation, the measurement
vocabulary, and — the point of it all — a case added from outside the agent
working end to end with no change to the agent's code.

    agent/.venv/bin/python agent/tests/test_plugins.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

AGENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(AGENT))

import case_templates  # noqa: E402
import cases  # noqa: E402
import measures  # noqa: E402
from postrun import Sample  # noqa: E402

FLAMELET = case_templates.TEMPLATE_BY_NAME["laminarFlamelet"]


def _problems(case) -> str:
    try:
        cases.validate(case, "test")
    except cases.CasePluginError as exc:
        return str(exc)
    return ""


def test_both_shipped_cases_are_plugins() -> None:
    names = [t.name for t in case_templates.TEMPLATES]
    assert names == ["laminarSlotBurner", "laminarFlamelet"], "ordered by `order`"
    assert case_templates.DEFAULT_TEMPLATE.name == "laminarSlotBurner"
    for name in names:
        assert (AGENT / "cases" / name / "case.py").is_file()


def test_the_shipped_plugins_are_valid() -> None:
    for template in case_templates.TEMPLATES:
        assert _problems(template) == "", template.name


def test_validation_names_each_inconsistency() -> None:
    broken = replace(
        FLAMELET, name="bad_name",
        extra_placeholders={"MYTOKEN": "not_a_parameter"},
        sweepable=("hole_size",),
        metric=None,
    )
    problems = _problems(broken)
    assert "no '_'" in problems
    assert "MYTOKEN" in problems
    assert "'hole_size' is sweepable" in problems
    assert "metric" in problems


def test_a_parameter_without_a_spec_is_rejected() -> None:
    problems = _problems(replace(FLAMELET, optional=FLAMELET.optional + ("swirl",)))
    assert "'swirl' has no ParamSpec" in problems


def test_a_missing_template_directory_is_rejected() -> None:
    problems = _problems(replace(FLAMELET, template_path=Path("/no/such/case")))
    assert "template directory not found" in problems


def test_crossing_position_measure() -> None:
    samples = [Sample(0.1, [0.0, 1.0, 2.0], {"T": [300.0, 1500.0, 1800.0]})]
    result = measures.crossing_position("T")(samples, 900.0)
    assert result["status"] == "measured"
    assert abs(result["crossing_position"] - 0.5) < 1e-12


def test_field_extremum_measure() -> None:
    samples = [Sample(0.1, [0.0, 1.0, 2.0], {"T": [300.0, 1900.0, 1800.0]})]
    result = measures.field_extremum("T", "max")(samples, 0.0)
    assert result["max_T"] == 1900.0 and result["max_T_position"] == 1.0


PLUGIN = '''from case_api import CaseTemplate, ParamSpec
from measures import field_extremum
CASE = CaseTemplate(
    name="myFlame",
    summary="a user-supplied case, reusing the flamelet's OpenFOAM files",
    required=("T_inlet",),
    optional=("sim_time", "sim_timestep", "sim_write_steps", "n_procs",
              "equivalence_ratio", "p_0", "ignition_label"),
    extra_placeholders={"INITIALTEMPERATURE": "T_inlet"},
    template_path=TEMPLATE,
    param_specs=(ParamSpec("ignition_label", float, False, 1.0, "-", "unused"),),
    sweepable=("T_inlet",),
    analyse=field_extremum("T", "max"),
    metric={"key": "max_T", "label": "Peak T [K]", "title": "Peak temperature"},
)
'''


def test_a_case_added_from_outside_the_agent_works_end_to_end() -> None:
    """Drop a folder on FGM_CASE_PATH; nothing in the agent is edited."""
    with tempfile.TemporaryDirectory() as tmp:
        plugin_dir = Path(tmp) / "myFlame"
        plugin_dir.mkdir()
        (plugin_dir / "case.py").write_text(
            PLUGIN.replace("TEMPLATE", repr(str(FLAMELET.template_dir)))
        )
        script = """
import json, tools, validator, campaign, case_templates, of_writer
out = {}
out["listed"] = [t.name for t in case_templates.TEMPLATES]
out["param"] = "ignition_label" in validator.SPEC_BY_NAME
r = tools.execute_tool("build_case", {"case_type": "myFlame", "T_inlet": 300.0,
                       "sim_timestep": 2e-6, "sim_time": 1e-4, "case_name": "plugtest"})
out["build"] = r["status"]; out["name"] = r.get("case_name")
out["sweep"] = campaign.check_variable("myFlame", "T_inlet")
out["metric"] = campaign.METRIC["myFlame"]["key"]
of_writer.delete_case("myFlame_plugtest")
print(json.dumps(out))
"""
        env = dict(os.environ, FGM_CASE_PATH=tmp, PYTHONPATH=str(AGENT))
        proc = subprocess.run([sys.executable, "-c", script], env=env,
                              capture_output=True, text=True, timeout=300)
        assert proc.returncode == 0, proc.stderr[-800:]
        out = json.loads(proc.stdout.strip().splitlines()[-1])
    assert "myFlame" in out["listed"]
    assert out["param"], "the plugin's own parameter must be known to the validator"
    assert out["build"] == "built" and out["name"] == "myFlame_plugtest"
    assert out["sweep"] is None, "sweepable variables come from the plugin"
    assert out["metric"] == "max_T"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
