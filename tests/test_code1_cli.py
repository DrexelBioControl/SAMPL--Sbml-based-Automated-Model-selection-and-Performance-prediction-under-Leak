"""Checks Code 1's run choices without starting a scientific fit."""

import argparse
import ast
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py"

# Runs the actual condition-setting lines without importing the full fitting script.
condition_setup = ast.parse(SCRIPT.read_text(encoding="utf-8"))
condition_setup.body = [
    node for node in condition_setup.body
    if isinstance(node, ast.Assign)
    and any(isinstance(target, ast.Name) and target.id in ("ALL_CONDITIONS", "FIT_CONDITIONS")
            for target in node.targets)
]

# Reads only the setup helpers; importing Code 1 would start fitting.
tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
tree.body = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name in ("read_command_line", "find_candidate_models", "read_output_prefix")
]
helpers = {"argparse": argparse}
exec(compile(tree, str(SCRIPT), "exec"), helpers)
read_command_line = helpers["read_command_line"]
find_candidate_models = helpers["find_candidate_models"]
read_output_prefix = helpers["read_output_prefix"]


class Code1CommandLineTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.base = Path(self.folder.name).resolve()
        self.config = self.base / "configs" / "fit_all_models.json"

    def test_omitted_fit_conditions_uses_every_condition(self):
        config = {"all_conditions": ["A", "B", "C"]}
        settings = {"config": config}
        exec(compile(condition_setup, str(SCRIPT), "exec"), settings)
        self.assertEqual(settings["FIT_CONDITIONS"], ["A", "B", "C"])
        self.assertIsNot(settings["FIT_CONDITIONS"], settings["ALL_CONDITIONS"])
        self.assertNotIn("fit_conditions", config)

    def test_explicit_fitting_subset_keeps_other_conditions_held_out(self):
        config = {"all_conditions": ["A", "B", "C"], "fit_conditions": ["A", "C"]}
        settings = {"config": config}
        exec(compile(condition_setup, str(SCRIPT), "exec"), settings)
        self.assertEqual(settings["FIT_CONDITIONS"], ["A", "C"])
        self.assertEqual(settings["ALL_CONDITIONS"], ["A", "B", "C"])
        self.assertIsNot(settings["FIT_CONDITIONS"], config["fit_conditions"])

    def test_omitted_options_preserve_manual_defaults(self):
        args = read_command_line(self.base, self.config, "load", [])
        self.assertEqual(args.config, self.config)
        self.assertEqual(args.mode, "load")
        self.assertIsNone(args.model)
        self.assertIsNone(args.model_name)

    def test_launcher_can_set_fold_config_and_one_model(self):
        args = read_command_line(self.base, self.config, "load", [
            "--config", "folds/without_A.json", "--mode", "fit",
            "--model", "models/chosen.xml",
        ])
        self.assertEqual(args.config, self.base / "folds" / "without_A.json")
        self.assertEqual(args.mode, "fit")
        self.assertEqual(args.model, [self.base / "models" / "chosen.xml"])

    def test_absolute_paths_and_repeated_models_are_preserved(self):
        config = self.base / "settings.json"
        model = self.base / "chosen.xml"
        args = read_command_line(self.base, self.config, "fit", [
            "--config", str(config), "--model", str(model),
            "--model", "other.xml",
        ])
        self.assertEqual(args.config, config)
        self.assertEqual(args.model, [model, self.base / "other.xml"])

    def test_invalid_mode_is_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            read_command_line(self.base, self.config, "fit", ["--mode", "unknown"])

    def test_archived_model_can_keep_original_name(self):
        args = read_command_line(self.base, self.config, "fit", [
            "--model", "results/selected_model.xml", "--model-name", "model1_activegate_fuel",
        ])
        self.assertEqual(args.model_name, "model1_activegate_fuel")
        self.assertEqual(args.model, [self.base / "results" / "selected_model.xml"])

    def test_model_name_override_requires_one_model(self):
        for options in ([], ["--model", "A.xml", "--model", "B.xml"]):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    read_command_line(self.base, self.config, "fit", options + ["--model-name", "A"])

    def test_model_name_override_cannot_create_a_folder(self):
        for name in ("", " ", "..", "../A", "A\\B", "C:A", "A\nB", "A?B"):
            with self.subTest(name=name), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    read_command_line(self.base, self.config, "fit", [
                        "--model", "selected_model.xml", "--model-name", name,
                    ])

    def test_output_prefix_is_optional_and_keeps_user_spelling(self):
        self.assertEqual(read_output_prefix({}), "")
        self.assertEqual(read_output_prefix({"output_prefix": "experiment_A_"}), "experiment_A_")
        self.assertEqual(read_output_prefix({"output_prefix": "test"}), "test")

    def test_output_prefix_cannot_create_folders_or_invalid_filenames(self):
        for prefix in (None, 4, "../", "C:\\results\\", "a/b", "bad?", "bad*", "a\nb"):
            with self.subTest(prefix=prefix), self.assertRaisesRegex(ValueError, "output_prefix"):
                read_output_prefix({"output_prefix": prefix})

    def test_default_discovery_and_selected_model_subset(self):
        for name in ("z.xml", "A.xml", "ignore.txt"):
            (self.base / name).touch()
        self.assertEqual(
            find_candidate_models(self.base, "*.xml"),
            [self.base / "A.xml", self.base / "z.xml"],
        )
        self.assertEqual(
            find_candidate_models(self.base, "*.xml", [self.base / "z.xml"]),
            [self.base / "z.xml"],
        )

    def test_selected_files_must_exist_and_be_files(self):
        for path in (self.base / "absent.xml", self.base):
            with self.subTest(path=path), self.assertRaises(FileNotFoundError):
                find_candidate_models(self.base, "*.xml", [path])

    def test_model_names_cannot_collide_in_output_files(self):
        first = self.base / "same.xml"
        second = self.base / "SAME.sbml"
        first.touch()
        second.touch()
        with self.assertRaisesRegex(ValueError, "distinct filename stems"):
            find_candidate_models(self.base, "*.xml", [first, second])

    def test_empty_model_folder_is_reported(self):
        with self.assertRaises(FileNotFoundError):
            find_candidate_models(self.base, "*.xml")

    def test_help_works_without_scientific_packages(self):
        result = subprocess.run(
            [sys.executable, "-S", "-B", str(SCRIPT), "--help"],
            cwd=str(self.base), capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for option in ("--config", "--mode", "--model", "--model-name"):
            self.assertIn(option, result.stdout)


if __name__ == "__main__":
    unittest.main()
