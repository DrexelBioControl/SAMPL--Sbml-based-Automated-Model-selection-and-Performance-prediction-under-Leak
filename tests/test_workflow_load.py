import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from workflow_steps import read_previous_identification


class PreviousIdentificationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name)
        self.old_run = self.project / "old_run"
        self.results = self.old_run / "identification" / "results"
        self.model_name = "model1_activegate_fuel"
        self.fit_folder = self.results / self.model_name
        self.fit_folder.mkdir(parents=True)
        self.archive = self.results / "selected_model.xml"
        self.archive.write_text("<saved_model/>", encoding="utf-8")
        self.fit_file = self.fit_folder / f"{self.model_name}_results.npz"
        self.fit_file.write_bytes(b"saved fit bytes")
        self.parameters_file = self.fit_folder / "fitted_parameters.csv"
        self.parameters_file.write_text("parameter,value\nk,2\n", encoding="utf-8")
        self.selection_file = self.results / "selected_model.json"
        self.selected = {
            "model_name": self.model_name,
            "model_file": str(self.project / "missing_models" / f"{self.model_name}.xml"),
            "fit_results_file": str(self.project / "missing_results" / "old_fit.npz"),
            "fitted_parameters_file": str(self.project / "missing_results" / "parameters.csv"),
            "observable_id": "ROL",
            "species_aliases": {"input": "IN_temp"},
            "fit_conditions": ["condition_1", "condition_2"],
        }
        self.write_selection()

    def write_selection(self):
        self.selection_file.write_text(json.dumps(self.selected), encoding="utf-8")

    def snapshot(self, folder):
        return {
            str(path.relative_to(folder)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in folder.rglob("*") if path.is_file()
        }

    def test_accepts_results_run_identification_and_selection_paths(self):
        for path in (self.results, self.old_run, self.results.parent, self.selection_file):
            with self.subTest(path=path):
                previous = read_previous_identification(self.project, path)
                self.assertEqual(previous["selection_file"], self.selection_file)
                self.assertEqual(previous["fit_file"], self.fit_file)
        previous = read_previous_identification(self.project, "old_run")
        self.assertEqual(previous["model_file"], self.archive)

    def test_local_saved_files_take_priority_over_existing_recorded_paths(self):
        current_model = self.project / "current.xml"
        current_fit = self.project / "current.npz"
        current_parameters = self.project / "current.csv"
        for path in (current_model, current_fit, current_parameters):
            path.write_bytes(b"different file")
        self.selected.update(model_file=str(current_model), fit_results_file=str(current_fit),
                             fitted_parameters_file=str(current_parameters))
        self.write_selection()
        previous = read_previous_identification(self.project, self.results)
        self.assertEqual(previous["model_file"], self.archive)
        self.assertEqual(previous["fit_file"], self.fit_file)
        self.assertEqual(previous["parameters_file"], self.parameters_file)

    def test_resolves_recorded_relative_and_absolute_files_when_local_copy_missing(self):
        self.archive.unlink()
        self.fit_file.unlink()
        project_model = self.project / "recorded.xml"
        project_model.write_text("<recorded_model/>", encoding="utf-8")
        results_fit = self.results / "recorded.npz"
        results_fit.write_bytes(b"recorded fit")
        self.selected.update(model_file="recorded.xml", fit_results_file="recorded.npz")
        self.write_selection()
        previous = read_previous_identification(self.project, self.results)
        self.assertEqual(previous["model_file"], project_model)
        self.assertEqual(previous["fit_file"], results_fit)
        self.selected.update(model_file=str(project_model), fit_results_file=str(results_fit))
        self.write_selection()
        self.assertEqual(read_previous_identification(self.project, self.results)["fit_file"], results_fit)

    def test_missing_selection_explains_required_file(self):
        with self.assertRaisesRegex(ValueError, "selected_model.json"):
            read_previous_identification(self.project, "absent")

    def test_missing_saved_model_does_not_use_current_model_by_name(self):
        self.archive.unlink()
        current_models = self.project / "models"
        current_models.mkdir()
        (current_models / f"{self.model_name}.xml").write_text("<different/>", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "saved SBML model"):
            read_previous_identification(self.project, self.results)

    def test_missing_saved_fit_explains_failure(self):
        self.fit_file.unlink()
        with self.assertRaisesRegex(ValueError, "saved fitted results"):
            read_previous_identification(self.project, self.results)

    def test_missing_readable_parameters_is_allowed(self):
        self.parameters_file.unlink()
        previous = read_previous_identification(self.project, self.results)
        self.assertIsNone(previous["parameters_file"])

    def test_checks_required_record_fields(self):
        for key in ("model_name", "model_file", "fit_results_file"):
            original = self.selected.pop(key)
            self.write_selection()
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, key):
                read_previous_identification(self.project, self.results)
            self.selected[key] = original

    def test_rejects_unsafe_model_names(self):
        for name in ("", "../escape", "nested/model", "nested\\model", ".", "..", "a:b", "trailing."):
            self.selected["model_name"] = name
            self.write_selection()
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "model_name"):
                read_previous_identification(self.project, self.results)

    def test_rejects_malformed_selection_json(self):
        for contents in ("not json", "[]"):
            self.selection_file.write_text(contents, encoding="utf-8")
            with self.subTest(contents=contents), self.assertRaisesRegex(ValueError, "selection"):
                read_previous_identification(self.project, self.results)

    def test_reading_keeps_sources_and_metadata_unchanged(self):
        before = self.snapshot(self.project)
        previous = read_previous_identification(self.project, self.old_run)
        self.assertEqual(previous["selected"], self.selected)
        self.assertEqual(before, self.snapshot(self.project))

    def test_shared_flat_results_are_found(self):
        results = self.project / "results"
        results.mkdir()
        archive = results / "selected_model.xml"
        archive.write_bytes(self.archive.read_bytes())
        fit = results / f"{self.model_name}_fit.npz"
        fit.write_bytes(self.fit_file.read_bytes())
        parameters = results / f"{self.model_name}_parameters.csv"
        parameters.write_bytes(self.parameters_file.read_bytes())
        selection = results / "selected_model.json"
        selection.write_text(json.dumps(self.selected), encoding="utf-8")
        previous = read_previous_identification(self.project, "results")
        self.assertEqual(previous["model_file"], archive)
        self.assertEqual(previous["fit_file"], fit)
        self.assertEqual(previous["parameters_file"], parameters)

    def test_prefixed_selection_uses_matching_archive_and_recorded_fit_basename(self):
        archive = self.results / "experiment_A_selected_model.xml"
        archive.write_text("<experiment_A/>", encoding="utf-8")
        fit = self.results / f"experiment_A_{self.model_name}_fit.npz"
        fit.write_bytes(b"experiment A fit")
        parameters = self.results / f"experiment_A_{self.model_name}_parameters.csv"
        parameters.write_text("parameter,value\nk,3\n", encoding="utf-8")
        selected = dict(self.selected, fit_results_file=str(self.project / "moved" / fit.name),
                        fitted_parameters_file=str(self.project / "moved" / parameters.name))
        selection = self.results / "experiment_A_selected_model.json"
        selection.write_text(json.dumps(selected), encoding="utf-8")
        previous = read_previous_identification(self.project, selection)
        self.assertEqual(previous["model_file"], archive)
        self.assertEqual(previous["fit_file"], fit)
        self.assertEqual(previous["parameters_file"], parameters)
        self.assertEqual(previous["selected"]["model_name"], self.model_name)

    def test_load_pointer_does_not_pick_unrelated_selected_archive(self):
        source_model = self.project / "previous_model.xml"
        source_model.write_text("<previous/>", encoding="utf-8")
        pointer = self.results / "loaded_model.json"
        selected = dict(self.selected, model_file=str(source_model))
        pointer.write_text(json.dumps(selected), encoding="utf-8")
        previous = read_previous_identification(self.project, pointer)
        self.assertEqual(previous["model_file"], source_model)

    def test_missing_prefixed_fit_does_not_load_an_unprefixed_fit(self):
        selection = self.results / "experiment_A_selected_model.json"
        selection.with_suffix(".xml").write_bytes(self.archive.read_bytes())
        selected = dict(self.selected, fit_results_file=str(self.results / "experiment_A_missing_fit.npz"))
        selection.write_text(json.dumps(selected), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "saved fitted results"):
            read_previous_identification(self.project, selection)

    def test_recorded_archive_is_used_before_current_model(self):
        source_archive = self.project / "archived.xml"
        source_archive.write_text("<previous/>", encoding="utf-8")
        pointer = self.results / "loaded_model.json"
        selected = dict(self.selected, archived_model_file=str(source_archive))
        pointer.write_text(json.dumps(selected), encoding="utf-8")
        previous = read_previous_identification(self.project, pointer)
        self.assertEqual(previous["model_file"], source_archive)

    def test_load_pointer_uses_exact_source_fit_despite_same_name_in_destination(self):
        source = self.project / "source"
        source.mkdir()
        source_fit = source / self.fit_file.name
        source_fit.write_bytes(b"chosen fit")
        selected = dict(self.selected, model_file=str(self.archive),
                        fit_results_file=str(source_fit), loaded_from="previous_selection.json")
        pointer = self.results / "loaded_model.json"
        pointer.write_text(json.dumps(selected), encoding="utf-8")
        previous = read_previous_identification(self.project, pointer)
        self.assertEqual(previous["fit_file"], source_fit)
        source_fit.unlink()
        with self.assertRaisesRegex(ValueError, "saved fitted results"):
            read_previous_identification(self.project, pointer)

    def test_does_not_choose_a_prefixed_selection_automatically(self):
        self.selection_file.rename(self.results / "other_selected_model.json")
        with self.assertRaisesRegex(ValueError, "selected_model.json"):
            read_previous_identification(self.project, self.results)

    def test_moving_saved_folder_preserves_the_original_fit_and_model(self):
        moved = self.project / "moved_run"
        shutil.copytree(str(self.old_run), str(moved))
        previous = read_previous_identification(self.project, moved)
        self.assertEqual(previous["model_file"], moved / "identification" / "results" / "selected_model.xml")
        self.assertTrue(previous["fit_file"].is_file())


if __name__ == "__main__":
    unittest.main()
