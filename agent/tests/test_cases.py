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

"""Tests for multi-template case setup.

Covers the rules that keep two (or more) template cases from interfering: each
declares its own required parameters and placeholder tokens, generated cases are
namespaced by template, and no template directory can ever be overwritten or
deleted.

The build tests write a real case under ``testcases/`` and remove it again —
they exercise the same code path the agent uses, and need no API key or network.

Runs with pytest, or standalone:

    agent/.venv/bin/python agent/tests/test_cases.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import case_templates  # noqa: E402
import of_writer  # noqa: E402
import tools  # noqa: E402
import validator  # noqa: E402
from case_templates import find_template  # noqa: E402

FLAMELET = case_templates.TEMPLATE_BY_NAME["laminarFlamelet"]
SLOT_BURNER = case_templates.TEMPLATE_BY_NAME["laminarSlotBurner"]

# Placeholder tokens that must never survive into a built case.
_TOKEN_PREFIXES = ("INLET", "PLATE", "INITIAL", "FINAL", "SIMULATION", "NPROCS")


def _unsubstituted(case_dir: Path) -> list[str]:
    """Working files still containing a placeholder token (0.orig/ excluded)."""
    offenders = []
    for sub in ("0", "system"):
        for path in sorted((case_dir / sub).rglob("*")):
            if not path.is_file():
                continue
            text = path.read_text()
            if any(token in text for token in _TOKEN_PREFIXES):
                offenders.append(str(path.relative_to(case_dir)))
    return offenders


def test_templates_declare_distinct_parameters() -> None:
    """The 1D flamelet has no inlet and no plate, so it must not ask for them."""
    assert "inlet_velocity" in SLOT_BURNER.required
    assert "T_plate" in SLOT_BURNER.required
    assert "inlet_velocity" not in FLAMELET.parameters
    assert "T_plate" not in FLAMELET.parameters
    # Both still need the unburnt temperature: it is the main physical knob and
    # has no defensible default. The run controls do have one (per case), so
    # they are optional — see test_run_controls_default_to_each_case_shipped_values.
    assert "T_inlet" in FLAMELET.required and "T_inlet" in SLOT_BURNER.required
    for name in ("sim_time", "sim_timestep"):
        for template in (FLAMELET, SLOT_BURNER):
            assert name in template.optional and name not in template.required


def test_templates_spell_the_temperature_differently() -> None:
    assert SLOT_BURNER.placeholders["INLETTEMPERATURE"] == "T_inlet"
    assert FLAMELET.placeholders["INITIALTEMPERATURE"] == "T_inlet"
    assert "INLETTEMPERATURE" not in FLAMELET.placeholders
    assert "INITIALTEMPERATURE" not in SLOT_BURNER.placeholders


def test_placeholders_match_the_real_template_files() -> None:
    """Guards against a template file gaining a token nothing substitutes."""
    for template in case_templates.TEMPLATES:
        found: set[str] = set()
        for sub in ("0.orig", "system.orig"):
            for path in (template.template_dir / sub).rglob("*"):
                if not path.is_file() or path.name.startswith("."):
                    continue
                for word in path.read_text().split():
                    token = word.strip(";(){}$\"',")
                    if token.isupper() and token.startswith(_TOKEN_PREFIXES):
                        found.add(token)
        declared = set(template.placeholders)
        assert found <= declared, (
            f"{template.name}: tokens in the case files with no replacement: "
            f"{sorted(found - declared)}"
        )


def test_validator_requirements_are_per_template() -> None:
    # The flamelet needs nothing beyond the temperature: the run controls fall
    # back to the values it ships with.
    assert validator.find_missing({"T_inlet": 300.0}, FLAMELET) == []
    # The slot burner still needs its inlet and plate, which have no default.
    assert sorted(validator.find_missing({"T_inlet": 300.0}, SLOT_BURNER)) == [
        "T_plate", "inlet_velocity",
    ]


def test_defaults_exclude_parameters_the_template_does_not_use() -> None:
    filled = validator.apply_defaults(
        {"T_inlet": 300.0, "sim_time": 0.1, "sim_timestep": 2e-6}, FLAMELET
    )
    assert "inlet_velocity" not in filled and "T_plate" not in filled
    assert filled["n_procs"] == 4 and filled["p_0"] == 101325.0


def test_temperature_aliases_resolve() -> None:
    """"initial temperature" is the natural phrasing for the flamelet."""
    for alias in ("T_initial", "T_unburnt", "T_0"):
        assert validator.coerce_types({alias: 300.0}) == {"T_inlet": 300.0}


def test_case_names_are_namespaced_by_template() -> None:
    assert of_writer.case_path("hot", FLAMELET).name == "laminarFlamelet_hot"
    assert of_writer.case_path("hot", SLOT_BURNER).name == "laminarSlotBurner_hot"
    # An already-prefixed name is respected whichever template is passed, so a
    # case can be addressed by name alone once it exists.
    assert of_writer.case_path("laminarFlamelet_hot", SLOT_BURNER).name == (
        "laminarFlamelet_hot"
    )


def test_template_directories_can_never_be_deleted() -> None:
    for name in ("laminarFlamelet", "laminarSlotBurner"):
        resolved = of_writer.case_path(name)
        assert resolved.name not in case_templates.PROTECTED_NAMES
        assert not resolved.exists(), "a template name must resolve to a prefixed case"


def test_unknown_case_type_is_rejected() -> None:
    try:
        find_template("laminarBunsen")
    except ValueError as exc:
        assert "laminarFlamelet" in str(exc)
    else:
        raise AssertionError("expected a ValueError for an unknown case type")


def test_ambiguous_request_asks_instead_of_defaulting() -> None:
    """"Set up a case at 300 K" must not silently become a slot burner.

    Quietly building the wrong geometry is worse than one extra question, and
    enforcing it in Python means it holds for every model.
    """
    result = tools.execute_tool("build_case", {
        "T_inlet": 300.0, "sim_time": 0.1, "sim_timestep": 2e-6,
    })
    assert result["status"] == "missing_case_type"
    assert result["missing"] == ["case_type"]
    assert {c["case_type"] for c in result["choices"]} == {
        "laminarSlotBurner", "laminarFlamelet",
    }
    # Nothing may be written while the question is outstanding.
    assert not (of_writer.TESTCASES_DIR / "laminarSlotBurner_agent").exists()


def test_case_type_is_inferred_when_parameters_can_only_mean_one() -> None:
    """An inlet velocity implies the slot burner — asking would be pedantic."""
    assert case_templates.infer_template(
        ["inlet_velocity", "T_inlet", "sim_time"]
    ) is SLOT_BURNER
    assert case_templates.infer_template(["T_plate", "sim_time"]) is SLOT_BURNER
    # Parameters common to both templates decide nothing.
    assert case_templates.infer_template(
        ["T_inlet", "sim_time", "sim_timestep", "n_procs", "p_0"]
    ) is None
    assert case_templates.infer_template([]) is None


def test_inference_survives_parameter_aliases() -> None:
    """The flamelet's natural phrasing must not be mistaken for a slot burner."""
    result = tools.execute_tool("build_case", {
        "T_initial": 300.0, "sim_time": 0.1, "sim_timestep": 2e-6,
    })
    assert result["status"] == "missing_case_type"


def test_explicit_case_type_always_wins_over_inference() -> None:
    result = tools.execute_tool("build_case", {
        "case_type": "laminarFlamelet",
        "T_inlet": 300.0, "sim_time": 0.1, "sim_timestep": 2e-6,
        "case_name": "pytest_explicit",
    })
    try:
        assert result["status"] == "built"
        assert result["case_type"] == "laminarFlamelet"
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_explicit")


def test_discriminating_parameters_are_actually_distinct() -> None:
    """Guards the inference: a template with no unique parameter can never win."""
    for template in case_templates.TEMPLATES:
        others = {
            name
            for other in case_templates.TEMPLATES if other.name != template.name
            for name in other.parameters
        }
        unique = case_templates.discriminating_parameters(template)
        assert not (unique & others)
    # The slot burner is distinguishable; the flamelet is the parameter subset,
    # so it is reached by naming it, not by inference.
    assert case_templates.discriminating_parameters(SLOT_BURNER)


def test_run_controls_default_to_each_case_shipped_values() -> None:
    """The recommended step differs per case — 1e-5 s vs 2e-6 s."""
    assert SLOT_BURNER.defaults == {"sim_time": 0.1, "sim_timestep": 1.0e-5}
    assert FLAMELET.defaults == {"sim_time": 0.1, "sim_timestep": 2.0e-6}
    for template in case_templates.TEMPLATES:
        filled = validator.apply_defaults({"T_inlet": 300.0}, template)
        assert filled["sim_timestep"] == template.defaults["sim_timestep"]
        assert filled["sim_time"] == template.defaults["sim_time"]
    # And they are no longer demanded up front.
    assert "sim_timestep" not in validator.find_missing({"T_inlet": 300.0}, FLAMELET)


def test_shipped_defaults_match_the_run_scripts() -> None:
    """The recommendation must stay the value the case actually ships with."""
    import re

    for template in case_templates.TEMPLATES:
        script = (template.template_dir / "run").read_text()
        for key, param in (("simtime", "sim_time"), ("simtimestep", "sim_timestep")):
            match = re.search(rf"^{key}=([0-9.eE+-]+)", script, re.M)
            assert match, f"{template.name}: no {key} in run"
            assert float(match.group(1)) == template.defaults[param], (
                f"{template.name}: {param} drifted from testcases/"
                f"{template.name}/run"
            )


def test_the_step_is_estimated_rather_than_taken_from_the_shipped_value() -> None:
    """The slot burner's velocity scale is known without solving a flame, so
    the step is computed from the mesh and the flow."""
    result = tools.execute_tool("build_case", {
        "inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0,
        "case_name": "pytest_dt",
    })
    try:
        estimate = result["sim_timestep_estimated"]
        assert result["parameters"]["sim_timestep"] == estimate["value"]
        assert estimate["target_courant"] == 0.1
        # The shipped rule of thumb is not repeated on top of a real estimate.
        assert "timestep_warning" not in result
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_dt")


def test_a_faster_inlet_gets_a_smaller_step() -> None:
    steps = {}
    for velocity in (0.2, 2.0):
        result = tools.execute_tool("build_case", {
            "inlet_velocity": velocity, "T_inlet": 300.0, "T_plate": 300.0,
            "case_name": "pytest_dt",
        })
        steps[velocity] = result["parameters"]["sim_timestep"]
        of_writer.delete_case("laminarSlotBurner_pytest_dt")
    assert steps[2.0] < steps[0.2]


def test_the_hot_mixture_rule_of_thumb_still_exists_as_a_fallback() -> None:
    """Used only where no estimate can be made; the advice itself is unchanged."""
    advice = case_templates.timestep_advice(FLAMELET, {"T_inlet": 600.0})
    assert advice is not None
    assert advice["suggested"] < FLAMELET.defaults["sim_timestep"]
    assert "adjustTimeStep no" in advice["message"]
    assert "Courant" in advice["message"]


def test_no_warning_at_the_reference_temperature_or_when_chosen() -> None:
    for args, why in (
        ({"T_inlet": 300.0}, "at the reference temperature"),
        ({"T_inlet": 600.0, "sim_timestep": 5e-7}, "step chosen explicitly"),
    ):
        result = tools.execute_tool("build_case", {
            "case_type": "laminarFlamelet", "case_name": "pytest_quiet", **args,
        })
        try:
            assert "timestep_warning" not in result, why
        finally:
            of_writer.delete_case("laminarFlamelet_pytest_quiet")


def test_recommended_values_are_reported_as_recommendations() -> None:
    """The model must be able to say which numbers it chose, not the user.

    The time step is estimated per case rather than taken from the shipped
    value, so here it is `sim_time` that falls back to the recommendation.
    """
    result = tools.execute_tool("build_case", {
        "case_type": "laminarFlamelet", "T_inlet": 300.0, "sim_timestep": 2e-6,
        "case_name": "pytest_reported",
    })
    try:
        used = result["used_recommended_values"]
        assert used == {"sim_time": 0.1}, "sim_timestep was the user's choice"
        assert "sim_timestep_estimated" not in result
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_reported")


def test_generated_cases_are_never_indexed() -> None:
    """A built case must not enter the knowledge corpus.

    It is rebuilt on every build, so indexing it fills the corpus with near
    duplicates of the template and marks the index stale each time. This once
    held only for the slot burner, because the prefix was hardcoded.
    """
    from rag.corpus import collect_corpus

    result = tools.execute_tool("build_case", {
        "case_type": "laminarFlamelet", "T_inlet": 300.0,
        "case_name": "pytest_corpus",
    })
    try:
        assert result["status"] == "built"
        indexed = {f.rel_path for f in collect_corpus()}
        assert not any("pytest_corpus" in path for path in indexed)
        # Every template's generated cases, not just one prefix.
        for template in case_templates.TEMPLATES:
            assert not any(
                path.startswith(f"testcases/{template.prefix}") for path in indexed
            ), template.name
        # The templates themselves are still indexed.
        assert any(
            path.startswith("testcases/laminarFlamelet/") for path in indexed
        )
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_corpus")


def test_geometry_defaults_reproduce_the_shipped_mesh() -> None:
    """Parameterising blockMeshDict must not change the case it ships as.

    These defaults are exactly the original hardcoded geometry
    (l_plate = 5 lF, h = 10 lF, h_plate = 4 lF, l_out = 30 lF, dx = lF/10), so a
    default build
    yields the same 43000-cell mesh as before.
    """
    filled = validator.apply_defaults(
        {"inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0}, SLOT_BURNER
    )
    lF = filled["flame_thickness"]
    assert lF == 6.2278e-04
    assert abs(filled["plate_thickness"] - 5 * lF) < 1e-12
    assert abs(filled["hole_pitch"] - 2 * 10 * lF) < 1e-12
    # hole_size = 2*(h - h_plate) = 2*(10 lF - 4 lF)
    assert abs(filled["hole_size"] - 2 * 6 * lF) < 1e-12
    assert abs(filled["outlet_length"] - 30 * lF) < 1e-12


def test_geometry_is_substituted_into_the_mesh_dictionary() -> None:
    result = tools.execute_tool("build_case", {
        "inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0,
        "plate_thickness": 6.2278e-3, "hole_size": 3.7e-3, "outlet_length": 0.035,
        "case_name": "pytest_geom",
    })
    try:
        assert result["status"] == "built"
        text = (Path(result["case_dir"]) / "system" / "blockMeshDict").read_text()
        assert "PLATETHICKNESS" not in text and "HOLESIZE" not in text
        assert "FLAMETHICKNESS" not in text and "HOLEPITCH" not in text
        assert "OUTLETLENGTH" not in text
        assert "l_out   0.035;" in text
        assert "0.0062278" in text          # the doubled plate thickness
        assert "0.0037" in text             # the reduced hole size
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_geom")


def test_the_flamelet_has_no_geometry_parameters() -> None:
    """Its blockMeshDict is not parameterised, so those must not leak in."""
    for name in ("plate_thickness", "hole_pitch", "hole_size", "outlet_length",
                 "flame_thickness"):
        assert name not in FLAMELET.parameters
    filled = validator.apply_defaults({"T_inlet": 300.0}, FLAMELET)
    assert not any(k in filled for k in ("plate_thickness", "hole_size"))


def test_build_flamelet_case_substitutes_everything() -> None:
    result = tools.execute_tool("build_case", {
        "case_type": "laminarFlamelet",
        "T_initial": 300.0, "sim_time": 0.1, "sim_timestep": 2e-6,
        "case_name": "pytest_tmp",
    })
    try:
        assert result["status"] == "built", result
        assert result["case_name"] == "laminarFlamelet_pytest_tmp"
        case_dir = Path(result["case_dir"])
        assert _unsubstituted(case_dir) == []
        # Conditions applied; mesh and BC types come from the template untouched.
        assert "uniform 300.0" in (case_dir / "0" / "T").read_text()
        assert "endTime           0.1;" in (case_dir / "system" / "controlDict").read_text()
        assert (case_dir / "constant" / "lookUp" / "database.fgm").is_file()
        # The flamelet has no inlet velocity to write.
        assert "uniform (0.0 0.0 0.0)" in (case_dir / "0" / "U").read_text()
    finally:
        of_writer.delete_case("laminarFlamelet_pytest_tmp")


def test_build_slot_burner_case_is_unaffected() -> None:
    result = tools.execute_tool("build_case", {
        "inlet_velocity": 0.2, "T_inlet": 300.0, "T_plate": 300.0,
        "sim_time": 0.1, "sim_timestep": 1e-5, "case_name": "pytest_tmp",
    })
    try:
        assert result["status"] == "built", result
        assert result["case_name"] == "laminarSlotBurner_pytest_tmp"
        case_dir = Path(result["case_dir"])
        assert _unsubstituted(case_dir) == []
        assert "uniform (0.2 0.0 0.0)" in (case_dir / "0" / "U").read_text()
    finally:
        of_writer.delete_case("laminarSlotBurner_pytest_tmp")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test in tests:
        test()
        print(f"  ok  {test.__name__}")
    print(f"\n{len(tests)} passed")
