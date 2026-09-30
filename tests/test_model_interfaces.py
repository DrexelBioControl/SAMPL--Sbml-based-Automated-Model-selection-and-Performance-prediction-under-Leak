"""Shared model defaults, run-specific choices, and explicit precedence."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from code2_configuration import load_model_interface


class ModelInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.file = self.root / "interfaces.json"
        self.shared = {
            "model_a": {
                "input_control": {"type": "species_initial_concentration", "id": "input"},
                "default_initial_species": {"gate": 25, "reporter": 500, "fuel": 0},
                "default_observable_units": "nM",
                "default_fixed_parameters": {"decay": 0.1, "other": 2.0},
            },
            "model_b": {
                "input_control": None,
                "default_initial_species": {"other_species": 3.0},
            },
        }
        self.write_shared()
        self.config = {
            "model_interfaces_file": "interfaces.json",
            "observable_id": "reporter",
            "analysis_control": {"type": "parameter", "id": "theta"},
            "initial_species_overrides": {},
            "fixed_parameters": {},
        }

    def write_shared(self):
        self.file.write_text(json.dumps(self.shared), encoding="utf-8")

    def load(self, model="model_a"):
        return load_model_interface(self.config, model, self.root)

    def test_defaults_and_analysis_choices_are_combined(self):
        result = self.load()
        self.assertEqual(result["fixed_initial_species"], {"gate": 25, "reporter": 500, "fuel": 0})
        self.assertEqual(result["observable_id"], "reporter")
        self.assertEqual(result["observable_units"], "nM")
        self.assertEqual(result["analysis_control"]["id"], "theta")

    def test_species_overrides_replace_only_named_defaults(self):
        self.config["initial_species_overrides"] = {"gate": 40, "new_species": 8}
        result = self.load()
        self.assertEqual(result["fixed_initial_species"], {"gate": 40, "reporter": 500, "fuel": 0, "new_species": 8})

    def test_parameter_overrides_retain_other_fixed_defaults(self):
        self.config["fixed_parameters"] = {"decay": 0.25}
        self.assertEqual(self.load()["fixed_parameters"], {"decay": 0.25, "other": 2.0})

    def test_units_can_be_overridden_for_one_run(self):
        self.config["observable_units"] = "uM"
        self.assertEqual(self.load()["observable_units"], "uM")

    def test_two_runs_can_choose_different_observables_and_controls(self):
        first = self.load()
        self.config["observable_id"] = "other_output"
        self.config["analysis_control"] = {"type": "species_initial_concentration", "id": "fuel"}
        second = self.load()
        self.assertEqual(first["fixed_initial_species"], second["fixed_initial_species"])
        self.assertNotEqual(first["observable_id"], second["observable_id"])
        self.assertNotEqual(first["analysis_control"], second["analysis_control"])

    def test_selected_model_uses_only_its_own_defaults(self):
        result = self.load("model_b")
        self.assertEqual(result["fixed_initial_species"], {"other_species": 3.0})
        self.assertIsNone(result["input_control"])
        self.assertEqual(result["observable_units"], "model concentration units")

    def test_loading_and_editing_result_do_not_change_config_or_shared_file(self):
        config_before = copy.deepcopy(self.config)
        shared_before = self.file.read_bytes()
        result = self.load()
        result["analysis_control"]["id"] = "edited"
        result["fixed_initial_species"]["gate"] = 999
        self.assertEqual(self.config, config_before)
        self.assertEqual(self.file.read_bytes(), shared_before)
        self.assertEqual(self.load()["fixed_initial_species"]["gate"], 25)

    def test_absolute_and_project_relative_paths_give_same_result(self):
        relative = self.load()
        self.config["model_interfaces_file"] = str(self.file)
        absolute = load_model_interface(self.config, "model_a", self.root / "unrelated")
        self.assertEqual(relative, absolute)

    def test_missing_model_and_missing_file_are_clear_errors(self):
        with self.assertRaisesRegex(ValueError, "No shared interface"):
            self.load("missing")
        self.config["model_interfaces_file"] = "missing.json"
        with self.assertRaises(FileNotFoundError):
            self.load()

    def test_old_inline_configuration_is_rejected_even_with_file_pointer(self):
        self.config["model_interfaces"] = self.shared
        with self.assertRaisesRegex(ValueError, "Move the inline"):
            self.load()

    def test_analysis_choices_cannot_silently_hide_in_shared_defaults(self):
        self.shared["model_a"]["analysis_control"] = {"type": "parameter", "id": "hidden"}
        self.write_shared()
        with self.assertRaisesRegex(ValueError, "Unsupported fields"):
            self.load()

    def test_invalid_run_choices_and_override_objects_are_rejected(self):
        original = copy.deepcopy(self.config)
        for key, value in [("observable_id", ""), ("analysis_control", None),
                           ("initial_species_overrides", []), ("fixed_parameters", []),
                           ("model_interfaces_file", "")]:
            with self.subTest(key=key):
                self.config = dict(original, **{key: value})
                with self.assertRaises(ValueError):
                    self.load()
        self.config = original
        self.shared["model_a"]["default_initial_species"] = []
        self.write_shared()
        with self.assertRaisesRegex(ValueError, "default_initial_species"):
            self.load()


if __name__ == "__main__":
    unittest.main()
