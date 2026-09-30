"""Checks prompt choices and the handoff between fitting, LOO and prediction."""
import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import Run_SAMPL as launcher


# These tiny child scripts test the workflow wiring, not SBML fitting or prediction.
FITTING_STUB = '''
import argparse, csv, json, shutil
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--config'); p.add_argument('--mode'); p.add_argument('--model', action='append')
p.add_argument('--model-name')
a = p.parse_args()
c = json.loads(Path(a.config).read_text())
assert a.mode == 'fit'
models = a.model or sorted(str(x) for x in Path(c['models_dir']).glob('*.xml'))
model = Path(models[0])
name = a.model_name or model.stem
prefix = c.get('output_prefix', '')
out = Path(c['results_dir']); out.mkdir(parents=True, exist_ok=True)
fit = out / (prefix + name + '_fit.npz'); fit.write_text(json.dumps(c['fit_conditions']))
save_selection = c.get('save_model_selection', True)
if save_selection:
    archive = out / (prefix + 'selected_model.xml')
    shutil.copyfile(str(model), str(archive))
    (out / (prefix + 'selected_model.json')).write_text(json.dumps({
        'model_name': name, 'model_file': str(model), 'archived_model_file': str(archive),
        'fit_results_file': str(fit)}))
(out / (prefix + 'observed_run.json')).write_text(json.dumps({
    'mode': a.mode, 'models': models, 'model_name': name}))
table = 'model_comparison.csv' if save_selection else 'fit_scores.csv'
with (out / (prefix + table)).open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['model', 'global_normalized_rmse', 'optimizer_success'])
    w.writeheader(); w.writerow({'model': name, 'global_normalized_rmse': 0.01, 'optimizer_success': True})
with (out / (prefix + name + '_held_out_rmse.csv')).open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['condition', 'raw_rmse', 'normalized_rmse', 'n_prediction_points'])
    w.writeheader()
    for name in c['all_conditions']:
        if name not in c['fit_conditions']:
            w.writerow({'condition': name, 'raw_rmse': 1, 'normalized_rmse': 0.1, 'n_prediction_points': 10})
'''

PREDICTION_STUB = '''
import argparse, json
from pathlib import Path
p = argparse.ArgumentParser(); p.add_argument('--config'); a = p.parse_args()
c = json.loads(Path(a.config).read_text())
selected = json.loads(Path(c['selected_model_file']).read_text())
fitted = json.loads(Path(selected['fit_results_file']).read_text())
stem = c.get('output_prefix', '') + selected['model_name'] + '_' + Path(c['metric']['definition_file']).stem
image = Path(c['output_file']) if c.get('output_file') else Path(c['figures_dir']) / (stem + '.png')
image.parent.mkdir(parents=True, exist_ok=True); image.write_bytes(b'prediction image')
out = Path(c['results_dir']); out.mkdir(parents=True, exist_ok=True)
(out / (image.stem + '.json')).write_text(json.dumps({
    'fitted_conditions': fitted, 'selected_model_file': c['selected_model_file']}))
(out / (image.stem + '.npz')).write_bytes(b'prediction values')
(out / (image.stem + '.settings.json')).write_text(json.dumps(c))
'''


class WorkflowLauncherTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.project = Path(folder.name).resolve()
        for name in ('configs', 'models', 'data', 'metrics', 'results'):
            (self.project / name).mkdir()
        for name in ('a', 'b', 'c'):
            (self.project / 'models' / (name + '.xml')).write_text('<sbml/>')
        (self.project / 'data/data.xlsx').write_bytes(b'unchanged workbook')
        self.conditions = ['C' + str(i) for i in range(6)]
        self.config = {
            'run_name': 'example', 'excel_file': 'data/data.xlsx', 'models_dir': 'models',
            'model_pattern': '*.xml', 'fit_conditions': self.conditions[:5],
            'all_conditions': self.conditions.copy(), 'conditions': {name: {} for name in self.conditions},
            'results_dir': 'results', 'figures_dir': 'figures', 'selection_metric': 'held_out_normalized_rmse',
        }
        self.fit_file = self.project / 'configs/fitting.json'
        launcher.save_json(self.fit_file, self.config)
        self.metric = self.project / 'metrics/example.py'
        self.metric.write_text('# metric placeholder')
        self.prediction = {
            'metric': {'name': 'example', 'definition_file': 'metrics/example.py', 'settings': {}},
            'model_interfaces_file': 'configs/interfaces.json', 'selected_model_file': 'results/old.json',
            'figures_dir': 'figures', 'results_dir': 'results', 'output_prefix': '',
        }
        self.prediction_file = self.project / 'configs/example.json'
        launcher.save_json(self.prediction_file, self.prediction)
        (self.project / 'results/selected_model.json').write_text('old result to preserve')
        (self.project / 'Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py').write_text(FITTING_STUB)
        (self.project / 'Code2_MasterSBMLPrediction_GeneralInterfaces.py').write_text(PREDICTION_STUB)

    def test_detects_models_from_config_and_validates_conditions(self):
        with contextlib.redirect_stdout(io.StringIO()) as output:
            models = launcher.detect_inputs(self.project, self.config)
        self.assertEqual(len(models), 3)
        self.assertIn('6 conditions specified in the fitting JSON', output.getvalue())
        self.config['fit_conditions'].append('missing')
        with self.assertRaises(ValueError):
            launcher.detect_inputs(self.project, self.config)

    def test_launcher_accepts_omitted_fitting_subset(self):
        del self.config['fit_conditions']
        with contextlib.redirect_stdout(io.StringIO()):
            models = launcher.detect_inputs(self.project, self.config)
        prepared = launcher.prepare_fit_config(self.project, self.config)
        self.assertEqual(len(models), 3)
        self.assertEqual(prepared['fit_conditions'], self.conditions)
        self.assertNotIn('fit_conditions', self.config)

    def test_yes_no_retries_instead_of_treating_typo_as_no(self):
        with patch('builtins.input', side_effect=['perhaps', 'YES']), contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(launcher.ask_yes_no('Run LOO?'))

    def test_only_matching_metric_settings_are_offered(self):
        (self.project / 'metrics/MetricTemplate.py').write_text('# example only')
        wrong = copy.deepcopy(self.prediction)
        wrong['metric']['definition_file'] = 'metrics/another.py'
        launcher.save_json(self.project / 'configs/wrong.json', wrong)
        with patch('builtins.input', side_effect=['example', '']), contextlib.redirect_stdout(io.StringIO()) as output:
            selected = launcher.choose_prediction(self.project)
        self.assertEqual(selected, self.prediction_file)
        self.assertNotIn('wrong.json', output.getvalue())
        self.assertNotIn('fitting.json', output.getvalue())
        self.assertNotIn('MetricTemplate.py', output.getvalue())

    def test_run_snapshots_fit_all_conditions_and_keep_originals_unchanged(self):
        original = copy.deepcopy(self.config)
        run = self.project / 'results'
        prepared = launcher.prepare_fit_config(self.project, self.config)
        self.assertEqual(prepared['fit_conditions'], self.conditions)
        self.assertEqual(prepared['selection_metric'], 'global_normalized_rmse')
        self.assertEqual(self.config, original)
        prediction = launcher.prepare_prediction_config(self.project, self.prediction, run / 'selected.json', prepared)
        self.assertEqual(prediction['selected_model_file'], str(run / 'selected.json'))
        self.assertEqual(self.prediction['selected_model_file'], 'results/old.json')
        for key in ('results_dir', 'figures_dir', 'generated_models_dir'):
            self.assertEqual(prediction[key], str(self.project / key[:-4]))

    def test_main_asks_mode_then_four_choices_and_uses_all_conditions_for_loo(self):
        with patch.object(launcher, '__file__', str(self.project / 'Run_SAMPL.py')), \
                patch('sys.argv', ['Run_SAMPL.py', '--fit-config', str(self.fit_file)]), \
                patch('builtins.input', side_effect=['', 'yes', 'yes', 'example', '']) as answers, \
                patch.object(launcher, 'execute_workflow') as execute, \
                contextlib.redirect_stdout(io.StringIO()):
            launcher.main()
        self.assertEqual(answers.call_count, 5)
        execute.assert_called_once_with(
            self.project, self.fit_file, self.config, self.conditions, self.prediction_file,
            mode='fit', previous_identification=None)

    def test_mode_accepts_fit_default_and_retries_invalid_answers(self):
        with patch('builtins.input', return_value=''):
            self.assertEqual(launcher.choose_mode(), 'fit')
        with patch('builtins.input', side_effect=['maybe', ' LOAD ']), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(launcher.choose_mode(), 'load')

    def test_check_does_not_prompt_or_run_calculations(self):
        with patch.object(launcher, '__file__', str(self.project / 'Run_SAMPL.py')), \
                patch('sys.argv', ['Run_SAMPL.py', '--fit-config', str(self.fit_file), '--check']), \
                patch('builtins.input') as answers, patch.object(launcher, 'execute_workflow') as execute, \
                contextlib.redirect_stdout(io.StringIO()) as output:
            launcher.main()
        answers.assert_not_called()
        execute.assert_not_called()
        self.assertIn('3 candidate models detected', output.getvalue())

    def test_complete_workflow_refits_six_folds_and_predicts_from_all_six(self):
        # The user's main JSON no longer needs to contain a manual fitting subset.
        del self.config['fit_conditions']
        launcher.save_json(self.fit_file, self.config)
        original_files = {path: path.read_bytes() for path in (
            self.fit_file, self.prediction_file, self.project / 'data/data.xlsx')}
        with contextlib.redirect_stdout(io.StringIO()) as output:
            run = launcher.execute_workflow(
                self.project, self.fit_file, self.config, self.conditions, self.prediction_file)
        summary = launcher.read_json(run / 'workflow_summary.json')
        self.assertEqual(run, self.project / 'results')
        self.assertEqual(summary['status'], 'completed')
        self.assertEqual(summary['completed_steps'], ['model_identification', 'loo', 'performance_prediction'])
        fit_settings = launcher.read_json(run / 'fitting_settings.json')
        self.assertEqual(fit_settings['fit_conditions'], self.conditions)
        identification = launcher.read_json(run / 'observed_run.json')
        self.assertEqual(len(identification['models']), 3)
        for index, held_out in enumerate(self.conditions, 1):
            prefix = 'loo_%02d_' % index
            config = launcher.read_json(run / (prefix + 'fitting_settings.json'))
            self.assertEqual(config['fit_conditions'], [name for name in self.conditions if name != held_out])
            self.assertEqual(config['all_conditions'], self.conditions)
            self.assertFalse(config['save_model_selection'])
            observed = launcher.read_json(run / (prefix + 'observed_run.json'))
            self.assertEqual(observed['models'], [str(run / 'selected_model.xml')])
            self.assertEqual(observed['model_name'], 'a')
            self.assertEqual(json.loads((run / (prefix + 'a_fit.npz')).read_text()), config['fit_conditions'])
            self.assertFalse((run / (prefix + 'selected_model.json')).exists())
        self.assertEqual(json.loads((run / 'a_fit.npz').read_text()), self.conditions)
        self.assertEqual(launcher.read_json(run / 'selected_model.json')['model_name'], 'a')
        self.assertTrue((self.project / 'figures/a_example.png').exists())
        self.assertTrue((run / 'a_example.npz').exists())
        self.assertTrue((run / 'a_loo_summary.csv').exists())
        self.assertFalse(any(path.is_dir() for path in run.iterdir()))
        prediction = launcher.read_json(run / 'a_example.json')
        self.assertEqual(prediction['fitted_conditions'], self.conditions)
        self.assertEqual(prediction['selected_model_file'], str(run / 'selected_model.json'))
        for path, contents in original_files.items():
            self.assertEqual(path.read_bytes(), contents)
        self.assertIn('Completed LOO.', output.getvalue())
        self.assertIn('Done.', output.getvalue())

    def test_identification_only_creates_no_loo_or_prediction(self):
        with contextlib.redirect_stdout(io.StringIO()):
            run = launcher.execute_workflow(self.project, self.fit_file, self.config, [], None)
        self.assertFalse((run / 'loo').exists())
        self.assertFalse((run / 'prediction').exists())
        self.assertFalse(list(run.glob('loo_*')))
        self.assertFalse((run / 'a_example.npz').exists())
        self.assertEqual(launcher.read_json(run / 'workflow_summary.json')['completed_steps'], ['model_identification'])

    def make_previous_identification(self):
        with contextlib.redirect_stdout(io.StringIO()):
            run = launcher.execute_workflow(self.project, self.fit_file, self.config, [], None)
        return run

    def test_load_mode_skips_identification_but_runs_six_fresh_loo_fits(self):
        previous_run = self.make_previous_identification()
        previous = launcher.read_previous_identification(self.project, previous_run)
        original = {previous_run / name: (previous_run / name).read_bytes() for name in (
            'a_fit.npz', 'selected_model.xml', 'selected_model.json', 'identification.log', 'observed_run.json')}
        with contextlib.redirect_stdout(io.StringIO()) as output:
            run = launcher.execute_workflow(
                self.project, self.fit_file, self.config, self.conditions, self.prediction_file,
                mode='load', previous_identification=previous)
        summary = launcher.read_json(run / 'workflow_summary.json')
        self.assertEqual(summary['identification_mode'], 'load')
        self.assertEqual(summary['completed_steps'], ['load_identification', 'loo', 'performance_prediction'])
        self.assertEqual(run, previous_run)
        loaded = launcher.read_json(run / 'loaded_model.json')
        self.assertEqual(loaded['loaded_from'], str(previous_run / 'selected_model.json'))
        self.assertEqual(loaded['fit_results_file'], str(previous_run / 'a_fit.npz'))
        self.assertEqual(loaded['model_file'], str(previous_run / 'selected_model.xml'))
        for index, held_out in enumerate(self.conditions, 1):
            prefix = 'loo_%02d_' % index
            settings = launcher.read_json(run / (prefix + 'fitting_settings.json'))
            observed = launcher.read_json(run / (prefix + 'observed_run.json'))
            self.assertEqual(observed['mode'], 'fit')
            self.assertEqual(len(observed['models']), 1)
            self.assertEqual(observed['models'], [str(run / 'selected_model.xml')])
            self.assertEqual(observed['model_name'], 'a')
            self.assertEqual(settings['fit_conditions'], [name for name in self.conditions if name != held_out])
            self.assertEqual(Path(settings['excel_file']), self.project / 'data/data.xlsx')
        prediction = launcher.read_json(run / 'a_example.json')
        self.assertEqual(prediction['fitted_conditions'], self.conditions)
        self.assertEqual(prediction['selected_model_file'], str(run / 'loaded_model.json'))
        for path, data in original.items():
            self.assertEqual(path.read_bytes(), data)
        self.assertIn('Loaded previous model identification', output.getvalue())
        self.assertNotIn('Running model identification.', output.getvalue())

    def test_load_prediction_only_does_not_require_original_workbook_or_candidate_folder(self):
        previous_run = self.make_previous_identification()
        (self.project / 'data/data.xlsx').unlink()
        for model in (self.project / 'models').glob('*.xml'):
            model.unlink()
        # Any accidental identification or LOO call must fail this test.
        (self.project / 'Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py').write_text('raise RuntimeError("Unexpected fitting")')
        with patch.object(launcher, '__file__', str(self.project / 'Run_SAMPL.py')), \
                patch('sys.argv', ['Run_SAMPL.py', '--fit-config', str(self.fit_file)]), \
                patch('builtins.input', side_effect=['load', str(previous_run), 'no', 'yes', 'example', '']), \
                contextlib.redirect_stdout(io.StringIO()):
            launcher.main()
        run = self.project / 'results'
        self.assertEqual(launcher.read_json(run / 'a_example.json')['fitted_conditions'], self.conditions)
        self.assertEqual(launcher.read_json(run / 'workflow_summary.json')['status'], 'completed')

    def test_missing_previous_identification_does_not_create_a_run(self):
        with self.assertRaisesRegex(ValueError, 'previous model identification'):
            launcher.execute_workflow(self.project, self.fit_file, self.config, [], None, mode='load')
        self.assertFalse((self.project / 'results/workflow_runs').exists())
        self.assertFalse((self.project / 'results/workflow_summary.json').exists())

    def test_previous_folder_prompt_retries_invalid_path(self):
        previous_run = self.make_previous_identification()
        with patch('builtins.input', side_effect=['missing_results', str(previous_run)]), \
                contextlib.redirect_stdout(io.StringIO()) as output:
            previous = launcher.choose_previous_identification(self.project)
        self.assertEqual(previous['selected']['model_name'], 'a')
        self.assertIn('Cannot use that folder', output.getvalue())

    def test_failed_step_stops_without_announcing_completion(self):
        with patch.object(launcher, 'run_script', side_effect=RuntimeError('Fitting failed')), \
                contextlib.redirect_stdout(io.StringIO()) as output, self.assertRaises(RuntimeError):
            launcher.execute_workflow(self.project, self.fit_file, self.config, self.conditions, self.prediction_file)
        run = self.project / 'results'
        summary = launcher.read_json(run / 'workflow_summary.json')
        self.assertEqual(summary['status'], 'stopped')
        self.assertEqual(summary['completed_steps'], [])
        self.assertNotIn('Completed', output.getvalue())
        self.assertNotIn('Done.', output.getvalue())
        self.assertFalse((run / 'selected_model.json').exists())

    def test_prefix_keeps_an_existing_identification_and_names_all_new_outputs(self):
        original_results = self.make_previous_identification()
        preserved = {original_results / name: (original_results / name).read_bytes()
                     for name in ('selected_model.json', 'selected_model.xml', 'a_fit.npz')}
        self.config['output_prefix'] = 'experiment_A_'
        with contextlib.redirect_stdout(io.StringIO()):
            results = launcher.execute_workflow(
                self.project, self.fit_file, self.config, self.conditions, self.prediction_file)
        self.assertEqual(results, original_results)
        for path, contents in preserved.items():
            self.assertEqual(path.read_bytes(), contents)
        self.assertTrue((results / 'experiment_A_selected_model.json').is_file())
        self.assertTrue((results / 'experiment_A_a_fit.npz').is_file())
        self.assertTrue((results / 'experiment_A_a_example.npz').is_file())
        self.assertTrue((self.project / 'figures/experiment_A_a_example.png').is_file())
        for index in range(1, 7):
            self.assertTrue((results / ('experiment_A_loo_%02d_a_fit.npz' % index)).is_file())
        selected = launcher.read_previous_identification(
            self.project, results / 'experiment_A_selected_model.json')
        self.assertEqual(selected['fit_file'], results / 'experiment_A_a_fit.npz')
        self.assertEqual(selected['model_file'], results / 'experiment_A_selected_model.xml')
        self.assertFalse(any(path.is_dir() for path in results.iterdir()))

    def test_repeating_a_calculation_replaces_matching_files_without_new_folders(self):
        with contextlib.redirect_stdout(io.StringIO()):
            first = launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file)
        filenames = sorted(path.name for path in first.iterdir())
        (first / 'a_fit.npz').write_text('stale fitted values')
        (first / 'a_example.npz').write_bytes(b'stale prediction')
        (self.project / 'figures/a_example.png').write_bytes(b'stale image')
        with contextlib.redirect_stdout(io.StringIO()):
            second = launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file)
        self.assertEqual(second, first)
        self.assertEqual(sorted(path.name for path in second.iterdir()), filenames)
        self.assertEqual(json.loads((second / 'a_fit.npz').read_text()), self.conditions)
        self.assertEqual((second / 'a_example.npz').read_bytes(), b'prediction values')
        self.assertEqual((self.project / 'figures/a_example.png').read_bytes(), b'prediction image')
        self.assertFalse(any(path.is_dir() for path in second.iterdir()))

    def test_loading_external_results_saves_only_a_pointer_and_prediction(self):
        previous_results = self.make_previous_identification()
        previous = launcher.read_previous_identification(self.project, previous_results)
        preserved = {path: path.read_bytes() for path in previous_results.iterdir() if path.is_file()}
        self.config['results_dir'] = 'another_results'
        self.config['figures_dir'] = 'another_figures'
        with contextlib.redirect_stdout(io.StringIO()):
            results = launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file,
                mode='load', previous_identification=previous)
        self.assertEqual(results, self.project / 'another_results')
        pointer = launcher.read_json(results / 'loaded_model.json')
        self.assertEqual(pointer['fit_results_file'], str(previous_results / 'a_fit.npz'))
        self.assertEqual(pointer['archived_model_file'], str(previous_results / 'selected_model.xml'))
        self.assertFalse((results / 'a_fit.npz').exists())
        self.assertFalse((results / 'selected_model.xml').exists())
        self.assertFalse((results / 'identification.log').exists())
        self.assertEqual(launcher.read_json(results / 'a_example.json')['fitted_conditions'], self.conditions)
        for path, contents in preserved.items():
            self.assertEqual(path.read_bytes(), contents)

    def test_explicit_image_filename_is_retained_with_numbers_in_results(self):
        self.prediction['output_file'] = 'figures/my_amplification.png'
        launcher.save_json(self.prediction_file, self.prediction)
        with contextlib.redirect_stdout(io.StringIO()):
            results = launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file)
        self.assertTrue((self.project / 'figures/my_amplification.png').is_file())
        self.assertTrue((results / 'my_amplification.npz').is_file())
        self.assertTrue((results / 'my_amplification.settings.json').is_file())
        self.assertFalse((self.project / 'figures/my_amplification.npz').exists())
        self.assertEqual(launcher.read_json(self.prediction_file), self.prediction)

    def test_successful_exit_without_new_selection_cannot_reuse_stale_result(self):
        self.make_previous_identification()
        with patch.object(launcher, 'run_script') as run_script, \
                patch.object(launcher, 'run_leave_one_out') as loo, \
                contextlib.redirect_stdout(io.StringIO()) as output, self.assertRaises(ValueError):
            launcher.execute_workflow(
                self.project, self.fit_file, self.config, self.conditions, self.prediction_file)
        self.assertEqual(run_script.call_count, 1)
        loo.assert_not_called()
        summary = launcher.read_json(self.project / 'results/workflow_summary.json')
        self.assertEqual(summary['status'], 'stopped')
        self.assertEqual(summary['completed_steps'], [])
        self.assertNotIn('Completed', output.getvalue())
        self.assertNotIn('Done.', output.getvalue())

    def test_loo_cannot_replace_the_loaded_fit_with_a_colliding_fold_filename(self):
        results = self.make_previous_identification()
        previous = launcher.read_previous_identification(self.project, results)
        colliding_fit = results / 'loo_01_a_fit.npz'
        colliding_fit.write_bytes(previous['fit_file'].read_bytes())
        previous['fit_file'] = colliding_fit
        previous['selected']['fit_results_file'] = str(colliding_fit)
        preserved = {path: path.read_bytes() for path in results.iterdir() if path.is_file()}
        with patch.object(launcher, 'run_script') as run_script, \
                patch.object(launcher, 'run_leave_one_out') as loo, \
                self.assertRaisesRegex(ValueError, 'replace the loaded fit'):
            launcher.execute_workflow(
                self.project, self.fit_file, self.config, self.conditions, self.prediction_file,
                mode='load', previous_identification=previous)
        run_script.assert_not_called()
        loo.assert_not_called()
        self.assertEqual({path: path.read_bytes() for path in results.iterdir() if path.is_file()}, preserved)

    def test_failed_prediction_preserves_successful_settings_and_keeps_pending_settings(self):
        with contextlib.redirect_stdout(io.StringIO()):
            results = launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file)
        final_settings = results / 'a_example.settings.json'
        pending_settings = results / 'a_example.pending_settings.json'
        original_settings = final_settings.read_bytes()
        original_values = (results / 'a_example.npz').read_bytes()
        self.assertFalse(pending_settings.exists())
        self.prediction['metric']['settings']['changed_setting'] = 123
        launcher.save_json(self.prediction_file, self.prediction)
        (self.project / 'Code2_MasterSBMLPrediction_GeneralInterfaces.py').write_text(
            'raise RuntimeError("Deliberate prediction failure")')
        with contextlib.redirect_stdout(io.StringIO()) as output, self.assertRaises(RuntimeError):
            launcher.execute_workflow(
                self.project, self.fit_file, self.config, [], self.prediction_file)
        self.assertEqual(final_settings.read_bytes(), original_settings)
        self.assertEqual((results / 'a_example.npz').read_bytes(), original_values)
        self.assertEqual(launcher.read_json(pending_settings)['metric']['settings']['changed_setting'], 123)
        summary = launcher.read_json(results / 'workflow_summary.json')
        self.assertEqual(summary['status'], 'stopped')
        self.assertEqual(summary['completed_steps'], ['model_identification'])
        self.assertNotIn('Completed performance', output.getvalue())
        self.assertNotIn('Done.', output.getvalue())


if __name__ == '__main__':
    unittest.main()
