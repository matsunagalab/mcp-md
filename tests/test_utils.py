"""Unit tests for mdzen.utils module.

Tests cover:
- validate_step_prerequisites: Workflow step prerequisite validation
- format_duration: Human-readable duration formatting
- safe_dict/safe_list: ADK state deserialization helpers
"""

import pytest

from mdzen.utils import (
    format_duration,
    safe_dict,
    safe_list,
)
from mdzen.workflow import (
    validate_step_prerequisites,
    get_current_step_info,
    SETUP_STEPS,
    STEP_CONFIG,
)


class TestValidateStepPrerequisites:
    """Tests for validate_step_prerequisites function."""

    def test_prepare_complex_no_prerequisites(self):
        """prepare_complex has no prerequisites."""
        is_valid, errors = validate_step_prerequisites("prepare_complex", {})
        assert is_valid is True
        assert errors == []

    def test_solvate_requires_merged_pdb(self):
        """solvate requires merged_pdb."""
        is_valid, errors = validate_step_prerequisites("solvate", {})
        assert is_valid is False
        assert any("merged_pdb" in e for e in errors)

    def test_solvate_with_merged_pdb(self):
        """solvate passes when merged_pdb is present."""
        is_valid, errors = validate_step_prerequisites(
            "solvate", {"merged_pdb": "/path/to/merged.pdb"}
        )
        assert is_valid is True
        assert errors == []

    def test_build_topology_requires_solvated_pdb_and_box(self):
        """build_topology requires solvated_pdb and box_dimensions."""
        is_valid, errors = validate_step_prerequisites("build_topology", {})
        assert is_valid is False
        assert any("solvated_pdb" in e for e in errors)
        assert any("box_dimensions" in e for e in errors)

    def test_build_topology_with_all_prerequisites(self):
        """build_topology passes with all prerequisites."""
        outputs = {
            "solvated_pdb": "/session/solvated.pdb",
            "box_dimensions": {"x": 80.0, "y": 80.0, "z": 80.0},
        }
        is_valid, errors = validate_step_prerequisites("build_topology", outputs)
        assert is_valid is True
        assert errors == []

    def test_run_simulation_requires_topology_files(self):
        """run_simulation requires prmtop and rst7."""
        is_valid, errors = validate_step_prerequisites("run_simulation", {})
        assert is_valid is False
        assert any("prmtop" in e for e in errors)
        assert any("rst7" in e for e in errors)

    def test_run_simulation_with_topology_files(self):
        """run_simulation passes with topology files."""
        outputs = {
            "prmtop": "/session/system.parm7",
            "rst7": "/session/system.rst7",
        }
        is_valid, errors = validate_step_prerequisites("run_simulation", outputs)
        assert is_valid is True
        assert errors == []

    def test_unknown_step_has_no_prerequisites(self):
        """Unknown steps should have no prerequisites."""
        is_valid, errors = validate_step_prerequisites("unknown_step", {})
        assert is_valid is True
        assert errors == []


class TestFormatDuration:
    """Tests for format_duration function."""

    def test_format_seconds_only(self):
        """Duration under 60 seconds should show seconds only."""
        assert format_duration(45.5) == "45.5s"
        assert format_duration(0.1) == "0.1s"

    def test_format_minutes_and_seconds(self):
        """Duration over 60 seconds should show minutes and seconds."""
        assert format_duration(90) == "1m 30s"
        assert format_duration(125.5) == "2m 6s"

    def test_format_exact_minute(self):
        """Exact minute boundaries should format correctly."""
        assert format_duration(60) == "1m 0s"
        assert format_duration(120) == "2m 0s"

    def test_format_hours(self):
        """Duration over 60 minutes should show hours."""
        assert format_duration(3600) == "1h 0m"
        assert format_duration(3720) == "1h 2m"


class TestSafeDict:
    """Tests for safe_dict function."""

    def test_dict_input(self):
        """Dict input should be returned as-is."""
        d = {"a": 1, "b": 2}
        assert safe_dict(d) == d

    def test_json_string_input(self):
        """JSON string input should be parsed."""
        s = '{"a": 1, "b": 2}'
        assert safe_dict(s) == {"a": 1, "b": 2}

    def test_none_input(self):
        """None should return empty dict."""
        assert safe_dict(None) == {}

    def test_none_with_default(self):
        """None with default should return default."""
        assert safe_dict(None, {"default": True}) == {"default": True}

    def test_invalid_json(self):
        """Invalid JSON should return default."""
        assert safe_dict("not json") == {}

    def test_empty_string(self):
        """Empty string should return default."""
        assert safe_dict("") == {}


class TestSafeList:
    """Tests for safe_list function."""

    def test_list_input(self):
        """List input should be returned as-is."""
        lst = [1, 2, 3]
        assert safe_list(lst) == lst

    def test_json_string_input(self):
        """JSON string input should be parsed."""
        s = '[1, 2, 3]'
        assert safe_list(s) == [1, 2, 3]

    def test_none_input(self):
        """None should return empty list."""
        assert safe_list(None) == []

    def test_none_with_default(self):
        """None with default should return default."""
        assert safe_list(None, ["default"]) == ["default"]

    def test_invalid_json(self):
        """Invalid JSON should return default."""
        assert safe_list("not json") == []

    def test_empty_string(self):
        """Empty string should return default."""
        assert safe_list("") == []


class TestGetCurrentStepInfo:
    """Tests for get_current_step_info function."""

    def test_empty_completed_steps(self):
        """Empty completed steps should return first step."""
        info = get_current_step_info([])
        assert info["current_step"] == SETUP_STEPS[0]
        assert info["step_index"] == 1

    def test_some_steps_completed(self):
        """Should return next incomplete step."""
        info = get_current_step_info(["prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"
        assert info["step_index"] == 3

    def test_all_steps_completed(self):
        """All steps completed should return is_complete=True."""
        info = get_current_step_info(SETUP_STEPS.copy())
        assert info["current_step"] is None
        assert info["next_tool"] is None
        assert info["is_complete"] is True

    def test_handles_duplicates(self):
        """Should handle duplicate step entries."""
        info = get_current_step_info(["prepare_complex", "prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"


class TestWorkflowMappings:
    """Tests for workflow step mappings."""

    def test_step_to_tool_mapping(self):
        """All steps should have tool mappings in STEP_CONFIG."""
        for step in SETUP_STEPS:
            assert step in STEP_CONFIG, f"Missing config for step: {step}"
            assert "tool" in STEP_CONFIG[step], f"Missing tool in config for step: {step}"

    def test_step_config_completeness(self):
        """All steps should have complete configuration."""
        required_keys = ["tool", "inputs", "servers", "allowed_tools", "estimate"]
        for step in SETUP_STEPS:
            for key in required_keys:
                assert key in STEP_CONFIG[step], f"Missing '{key}' in config for step: {step}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
