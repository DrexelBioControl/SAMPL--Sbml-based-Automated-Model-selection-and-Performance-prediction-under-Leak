import copy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

import workflow_steps


def write_csv(filename, rows):
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class TestRunScript(unittest.TestCase):
    def test_uses_same_python_separate_arguments_and_captures_both_streams(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            log_file = project / "logs" / "step.log"

            def run(command, **kwargs):
                self.assertEqual(command, [
                    sys.executable, "-u", str(project / "Code 1.py"),
                    "--config", str(project / "some settings.json"),
                    "--mode", "fit", "--model", "model with spaces.xml",
                ])
                self.assertEqual(kwargs["cwd"], str(project))
                self.assertEqual(kwargs["env"]["MPLBACKEND"], "Agg")
                self.assertTrue(kwargs["check"])
                self.assertEqual(kwargs["stderr"], subprocess.STDOUT)
                self.assertNotIn("shell", kwargs)
                kwargs["stdout"].write("step output\n")

            with patch.object(workflow_steps.subprocess, "run", side_effect=run) as process:
                workflow_steps.run_script(
                    project, "Code 1.py", "some settings.json", log_file,
                    ["--mode", "fit", "--model", "model with spaces.xml"],
                )
            process.assert_called_once()
            self.assertEqual(log_file.read_text(encoding="utf-8"), "step output\n")

    def test_failed_child_reports_the_log_and_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            log_file = Path(temporary) / "step.log"

            def fail(command, **kwargs):
                kwargs["stdout"].write("the underlying error\n")
                raise subprocess.CalledProcessError(2, command)

            with patch.object(workflow_steps.subprocess, "run", side_effect=fail):
                with self.assertRaisesRegex(RuntimeError, "did not complete") as raised:
                    workflow_steps.run_script(temporary, "main.py", "settings.json", log_file)
            self.assertIn(str(log_file), str(raised.exception))
            self.assertIn("underlying error", log_file.read_text(encoding="utf-8"))

    def test_failed_process_start_is_readable(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(workflow_steps.subprocess, "run", side_effect=OSError("cannot start")):
                with self.assertRaisesRegex(RuntimeError, "cannot start"):
                    workflow_steps.run_script(temporary, "main.py", "settings.json", "step.log")


class TestLeaveOneOut(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name)
        self.model = self.project / "selected model.xml"
        self.model.write_text("unchanged model", encoding="utf-8")
        self.excel = self.project / "measurements.xlsx"
        self.excel.write_bytes(b"unchanged measurements")
        self.config = {
            "run_name": "identification",
            "excel_file": str(self.excel),
            "models_dir": str(self.project),
            "fit_conditions": ["A", "B", "C"],
            "all_conditions": ["A", "B", "C", "external_test"],
            "conditions": {
                name: {"initial_values": {"X": 1}}
                for name in ("A", "B", "C", "external_test")
            },
            "results_dir": str(self.project / "identification_results"),
            "figures_dir": str(self.project / "identification_figures"),
            "generated_models_dir": str(self.project / "generated_models"),
            "selection_metric": "held_out_normalized_rmse",
        }
        self.original_config = copy.deepcopy(self.config)
        self.original_result = self.project / "identification_results" / "selected_model.json"
        self.original_result.parent.mkdir()
        self.original_result.write_text("unchanged selection", encoding="utf-8")
        self.output = Path(self.config["results_dir"])
        self.summary = self.output / f"{self.model.stem}_loo_summary.csv"
        self.seen_configs = []

    def fake_fit(self, project_dir, script_name, config_file, log_file, extra_args=None):
        self.assertEqual(Path(project_dir), self.project)
        self.assertEqual(script_name, "Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py")
        self.assertEqual(extra_args, ["--mode", "fit", "--model", str(self.model), "--model-name", self.model.stem])
        config = json.loads(Path(config_file).read_text(encoding="utf-8"))
        self.seen_configs.append(config)
        held_out = [name for name in config["all_conditions"] if name not in config["fit_conditions"]]
        self.assertEqual(len(held_out), 1)
        self.assertNotIn("external_test", config["all_conditions"])
        self.assertEqual(config["selection_metric"], self.config["selection_metric"])
        self.assertFalse(config["save_model_selection"])
        self.assertEqual(config["generated_models_dir"], self.config["generated_models_dir"])
        self.assertEqual(config["excel_file"], str(self.excel))
        self.assertTrue(Path(config["results_dir"]).is_absolute())
        self.assertTrue(Path(config["figures_dir"]).is_absolute())
        self.assertEqual(Path(config["results_dir"]), Path(config_file).parent)
        self.assertEqual(config["figures_dir"], self.config["figures_dir"])
        self.assertEqual(Path(log_file).parent, Path(config_file).parent)
        error = {"A": 0.1, "B": 0.2, "C": 0.3, "D": 0.4, "E": 0.5, "F": 0.6}[held_out[0]]
        prefix = config["output_prefix"]
        (self.output / f"{prefix}{self.model.stem}_fit.npz").write_bytes(b"fold fit")
        write_csv(self.output / f"{prefix}{self.model.stem}_held_out_rmse.csv", [{
            "condition": held_out[0], "raw_rmse": 10 * error, "normalized_rmse": error,
        }])
        write_csv(self.output / f"{prefix}fit_scores.csv", [{
            "model": self.model.stem, "global_normalized_rmse": 0.05,
            "held_out_normalized_rmse": error, "optimizer_success": held_out[0] != "B",
        }])

    def run_loo(self, condition_names=None, fake_fit=None):
        with patch.object(workflow_steps, "run_script", side_effect=fake_fit or self.fake_fit):
            with redirect_stdout(io.StringIO()):
                return workflow_steps.run_leave_one_out(
                    self.project, self.config,
                    condition_names if condition_names is not None else ["A", "B", "C"],
                    self.model,
                )

    def test_each_fold_refits_only_remaining_conditions_without_changing_originals(self):
        summary_file = self.run_loo()
        self.assertEqual(summary_file, self.summary)
        self.assertEqual([row["fit_conditions"] for row in self.seen_configs], [
            ["B", "C"], ["A", "C"], ["A", "B"],
        ])
        self.assertTrue(all(row["all_conditions"] == ["A", "B", "C"] for row in self.seen_configs))
        self.assertEqual(len({row["results_dir"] for row in self.seen_configs}), 1)
        self.assertEqual([row["output_prefix"] for row in self.seen_configs], ["loo_01_", "loo_02_", "loo_03_"])
        self.assertFalse(any(path.is_dir() for path in self.output.iterdir()))
        self.assertEqual(self.config, self.original_config)
        self.assertEqual(self.excel.read_bytes(), b"unchanged measurements")
        self.assertEqual(self.original_result.read_text(encoding="utf-8"), "unchanged selection")
        self.assertEqual(self.model.read_text(encoding="utf-8"), "unchanged model")

    def test_summaries_keep_each_held_out_condition_error_without_averaging(self):
        summary_file = self.run_loo()
        summary = json.loads(summary_file.with_suffix(".json").read_text(encoding="utf-8"))
        self.assertEqual(summary["condition_pool"], ["A", "B", "C"])
        self.assertEqual(summary["number_of_folds"], 3)
        self.assertNotIn("mean_condition_normalized_rmse", summary)
        self.assertNotIn("score_description", summary)
        self.assertEqual([row["normalized_rmse"] for row in summary["folds"]], [0.1, 0.2, 0.3])
        self.assertEqual([row["optimizer_success"] for row in summary["folds"]], [True, False, True])
        with summary_file.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["held_out_condition"] for row in rows], ["A", "B", "C"])
        self.assertEqual(json.loads(rows[0]["fit_conditions"]), ["B", "C"])

    def test_custom_pool_excludes_every_condition_outside_it(self):
        self.run_loo(["C", "A"])
        self.assertEqual([row["fit_conditions"] for row in self.seen_configs], [["A"], ["C"]])
        self.assertEqual(self.seen_configs[0]["all_conditions"], ["C", "A"])

    def test_six_condition_pool_refits_five_and_withholds_one_in_every_fold(self):
        names = ["A", "B", "C", "D", "E", "F"]
        self.config["all_conditions"] = names
        self.config["fit_conditions"] = names.copy()
        self.config["conditions"].update({name: {} for name in ("D", "E", "F")})
        self.run_loo(names)
        self.assertEqual(len(self.seen_configs), 6)
        for held_out, config in zip(names, self.seen_configs):
            self.assertEqual(config["all_conditions"], names)
            self.assertEqual(len(config["fit_conditions"]), 5)
            self.assertEqual(set(config["fit_conditions"]), set(names) - {held_out})

    def test_invalid_pool_is_rejected_before_running_or_writing(self):
        for names in ([], ["A"], ["A", "A"], ["A", "missing"], ["A", 2]):
            with self.subTest(names=names):
                with self.assertRaises(ValueError):
                    self.run_loo(names)
                self.assertEqual(list(self.output.iterdir()), [self.original_result])
        self.assertEqual(self.seen_configs, [])

    def test_condition_must_have_settings_and_be_in_all_conditions(self):
        del self.config["conditions"]["B"]
        with self.assertRaisesRegex(ValueError, "missing"):
            self.run_loo()
        self.assertEqual(list(self.output.iterdir()), [self.original_result])

    def test_repeat_replaces_fold_outputs_but_preserves_the_identification(self):
        self.run_loo()
        self.summary.write_text("old summary", encoding="utf-8")
        self.run_loo()
        self.assertEqual(len(self.seen_configs), 6)
        self.assertIn("held_out_condition", self.summary.read_text(encoding="utf-8"))
        self.assertEqual(self.original_result.read_text(encoding="utf-8"), "unchanged selection")
        self.assertFalse(any(path.is_dir() for path in self.output.iterdir()))

    def test_previous_success_is_not_reported_after_a_failed_rerun(self):
        self.run_loo()
        with self.assertRaisesRegex(RuntimeError, "Cannot read workflow results"):
            self.run_loo(fake_fit=lambda *args, **kwargs: None)
        self.assertFalse(self.summary.exists())
        self.assertFalse(self.summary.with_suffix(".json").exists())

    def test_optional_prefix_applies_to_fold_files_and_summary(self):
        self.config["output_prefix"] = "experiment_A_"
        summary = self.run_loo()
        self.assertEqual(summary.name, f"experiment_A_{self.model.stem}_loo_summary.csv")
        self.assertEqual(self.seen_configs[0]["output_prefix"], "experiment_A_loo_01_")
        self.assertTrue((self.output / "experiment_A_loo_01_fitting_settings.json").is_file())

    def test_saved_model_archive_keeps_original_model_identity(self):
        archive = self.project / "selected_model.xml"
        archive.write_text("<model/>", encoding="utf-8")
        def fit(project_dir, script, config_file, log_file, extra_args):
            self.assertEqual(extra_args, ["--mode", "fit", "--model", str(archive), "--model-name", self.model.stem])
            self.fake_fit(project_dir, script, config_file, log_file,
                          ["--mode", "fit", "--model", str(self.model), "--model-name", self.model.stem])
        with patch.object(workflow_steps, "run_script", side_effect=fit):
            summary = workflow_steps.run_leave_one_out(self.project, self.config, ["A", "B"], archive, self.model.stem)
        self.assertEqual(summary, self.summary)

    def test_failed_fold_stops_without_a_completion_summary(self):
        def fail_second(*args, **kwargs):
            if len(self.seen_configs) == 1:
                raise RuntimeError("second fit failed")
            self.fake_fit(*args, **kwargs)

        with self.assertRaisesRegex(RuntimeError, "second fit failed"):
            self.run_loo(fake_fit=fail_second)
        self.assertEqual(len(self.seen_configs), 1)
        self.assertTrue((self.output / f"loo_01_{self.model.stem}_fit.npz").is_file())
        self.assertFalse((self.output / "loo_03_fitting_settings.json").exists())
        self.assertFalse((self.summary).exists())
        self.assertFalse((self.summary.with_suffix(".json")).exists())

    def test_wrong_or_extra_held_out_conditions_are_rejected(self):
        def wrong_condition(*args, **kwargs):
            self.fake_fit(*args, **kwargs)
            filename = self.output / f"{self.seen_configs[-1]['output_prefix']}{self.model.stem}_held_out_rmse.csv"
            write_csv(filename, [
                {"condition": "A", "raw_rmse": 1, "normalized_rmse": 0.1},
                {"condition": "external_test", "raw_rmse": 2, "normalized_rmse": 0.2},
            ])

        with self.assertRaisesRegex(RuntimeError, "Expected only held-out condition"):
            self.run_loo(fake_fit=wrong_condition)
        self.assertFalse((self.summary).exists())

    def test_extra_candidate_models_are_rejected(self):
        def extra_model(*args, **kwargs):
            self.fake_fit(*args, **kwargs)
            write_csv(self.output / f"{self.seen_configs[-1]['output_prefix']}fit_scores.csv", [
                {"model": self.model.stem, "global_normalized_rmse": 0.05},
                {"model": "other_model", "global_normalized_rmse": 0.05},
            ])

        with self.assertRaisesRegex(RuntimeError, "Expected only selected model"):
            self.run_loo(fake_fit=extra_model)

    def test_nonfinite_prediction_error_stops_the_run(self):
        def invalid_error(*args, **kwargs):
            self.fake_fit(*args, **kwargs)
            filename = self.output / f"{self.seen_configs[-1]['output_prefix']}{self.model.stem}_held_out_rmse.csv"
            write_csv(filename, [{"condition": "A", "raw_rmse": 1, "normalized_rmse": "nan"}])

        with self.assertRaisesRegex(RuntimeError, "Non-finite"):
            self.run_loo(fake_fit=invalid_error)
        self.assertFalse((self.summary).exists())

    def test_missing_error_table_stops_the_run(self):
        with self.assertRaisesRegex(RuntimeError, "Cannot read workflow results"):
            self.run_loo(fake_fit=lambda *args, **kwargs: None)
        self.assertFalse((self.summary).exists())


if __name__ == "__main__":
    unittest.main()
