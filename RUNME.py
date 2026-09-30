#---------------------------
# SAMPL: FITTING, OPTIONAL LOO, AND PERFORMANCE PREDICTION
#---------------------------
# Input: fitting JSON, data/models or a saved identification, and metric settings.
# Output: figures/ for images; results/ for numbers, settings and logs.
#
# 1. Asks which fitting JSON to use, then whether to fit or load.
# 2. Asks whether to run LOO and performance prediction.
# 3. Fits or loads the selected model, then runs the requested extra steps.
# 4. Saves using fixed filenames; output_prefix can name a separate set of results.
#
# Run this file. The numerical calculations stay in Code 1 and Code 2.
#---------------------------


#------------------------------------------------------------------------------------------
import argparse                                                                           #|
import copy                                                                               #|
import json                                                                               #|
from pathlib import Path                                                                  #|
                                                                                          #|
from workflow_steps import (run_script, run_leave_one_out, read_previous_identification)  #|
#------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# 1) Reads files and choices using the project folder as the starting point.         #|
                                                                                     #|
def project_path(project_dir, filename):                                             #|
    """Return an absolute path, including when the filename is project-relative."""  #|
    return (Path(project_dir) / filename).resolve()                                  #|
#-------------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def read_json(filename):                                                      #|
    with Path(filename).open('r', encoding='utf-8') as handle:                #|
        return json.load(handle)                                              #|
#------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
def save_json(filename, values):                                               #|
    Path(filename).parent.mkdir(parents=True, exist_ok=True)                   #|
    Path(filename).write_text(json.dumps(values, indent=2), encoding='utf-8')  #|
#-------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------
def ask_yes_no(question, default=False):                                                            #|
    """Keeps asking until the answer is yes or no."""                                               #|
    default_text = 'yes' if default else 'no'                                                       #|
    while True:                                                                                     #|
        answer = input('{} [yes/no; default {}]: '.format(question, default_text)).strip().lower()  #|
        if not answer:                                                                              #|
            return default                                                                          #|
        if answer in ('yes', 'y'):                                                                  #|
            return True                                                                             #|
        if answer in ('no', 'n'):                                                                   #|
            return False                                                                            #|
        print('Please type yes or no.')                                                             #|
#----------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
def choose_file(question, files, default=None):                                        #|
    """Accepts a listed filename, its name without the extension, or its number."""    #|
    if not files:                                                                      #|
        raise ValueError('No suitable files were found for: ' + question)              #|
    default = default if default in files else files[0]                                #|
    for index, path in enumerate(files, 1):                                            #|
        print('  {}. {}'.format(index, path.name))                                     #|
    while True:                                                                        #|
        answer = input('{} [{}]: '.format(question, default.name)).strip().strip('"')  #|
        if not answer:                                                                 #|
            return default                                                             #|
        for index, path in enumerate(files, 1):                                        #|
            if answer.lower() in (str(index), path.name.lower(), path.stem.lower()):   #|
                return path                                                            #|
        print('Please enter one of the listed filenames or numbers.')                  #|
#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
def detect_inputs(project_dir, config, model_file=None):                                                 #|
    """Checks the workbook and lists the models Code 1 will fit."""                                      #|
    excel_file = project_path(project_dir, config['excel_file'])                                         #|
    if not excel_file.is_file():                                                                         #|
        raise FileNotFoundError('The configured Excel data file was not found: ' + str(excel_file))      #|
    if model_file is None:                                                                               #|
        models_dir = project_path(project_dir, config['models_dir'])                                     #|
        model_files = sorted(                                                                            #|
            (path for path in models_dir.glob(config.get('model_pattern', '*.xml')) if path.is_file()),  #|
            key=lambda path: path.name.lower(),                                                          #|
        )                                                                                                #|
        if not model_files:                                                                              #|
            raise ValueError('No candidate models matched the fitting settings in: ' + str(models_dir))  #|
    else:                                                                                                #|
        # LOO in load mode needs only the selected model and the configured data.                        #|
        model_files = [Path(model_file)]                                                                 #|
        if not model_files[0].is_file():                                                                 #|
            raise FileNotFoundError('The saved model file was not found: ' + str(model_file))            #|
    all_conditions = config['all_conditions']                                                            #|
    # Accepts settings with no manual fitting subset.                                                    #|
    fit_conditions = config.get('fit_conditions', all_conditions)                                        #|
    for label, names in (('fit_conditions', fit_conditions), ('all_conditions', all_conditions)):        #|
        if not names or len(names) != len(set(names)):                                                   #|
            raise ValueError(label + ' must contain a nonempty list of distinct conditions.')            #|
        if any(name not in config['conditions'] for name in names):                                      #|
            raise ValueError(label + ' contains a condition missing from the conditions settings.')      #|
    if not set(fit_conditions).issubset(all_conditions):                                                 #|
        raise ValueError('Every fit condition must also appear in all_conditions.')                      #|
    print('\nData detected: ' + str(excel_file))                                                         #|
    if model_file is None:                                                                               #|
        print('{} candidate models detected.'.format(len(model_files)))                                  #|
    else:                                                                                                #|
        print('Selected model available for fresh LOO fits:')                                            #|
    for path in model_files:                                                                             #|
        print('  ' + path.name)                                                                          #|
    print('{} conditions specified in the fitting JSON.'.format(len(all_conditions)))                    #|
    return model_files                                                                                   #|
#---------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------
# 2) Chooses the fitting JSON, fit/load mode and optional steps.                                           #|
                                                                                                           #|
def choose_fit_config(project_dir):                                                                        #|
    """Asks for a fitting JSON filename or path before starting the workflow."""                           #|
    print('\nFitting settings folder: ' + str(Path(project_dir) / 'configs'))                              #|
    while True:                                                                                            #|
        answer = input('Fitting settings filename [fit_all_models.json]: ').strip().strip('"').strip("'")  #|
        filename = Path(answer or 'fit_all_models.json').expanduser()                                      #|
        # Accepts a filename with or without the .json extension.                                          #|
        if not filename.suffix:                                                                            #|
            filename = filename.with_suffix('.json')                                                       #|
        # Looks for a bare filename in configs/; also accepts a full or project-relative path.             #|
        if filename.parent == Path('.'):                                                                   #|
            filename = Path('configs') / filename                                                          #|
        fit_file = project_path(project_dir, filename)                                                     #|
        try:                                                                                               #|
            settings = read_json(fit_file)                                                                 #|
            # Distinguishes fitting settings from metric settings and model defaults.                      #|
            required = ('excel_file', 'models_dir', 'all_conditions', 'conditions')                        #|
            if not isinstance(settings, dict) or any(key not in settings for key in required):             #|
                raise ValueError('Choose a fitting JSON containing: ' + ', '.join(required))               #|
        except (OSError, ValueError) as error:                                                             #|
            print('Cannot use that settings file: ' + str(error))                                          #|
            continue                                                                                       #|
        return fit_file                                                                                    #|
#-----------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
def choose_mode():                                                                             #|
    """Defaults to a new fit; load reuses the selected model and saved parameters."""          #|
    while True:                                                                                #|
        mode = input('\nModel identification mode [fit/load; default fit]: ').strip().lower()  #|
        if mode in ('', 'fit'):                                                                #|
            return 'fit'                                                                       #|
        if mode == 'load':                                                                     #|
            return 'load'                                                                      #|
        print('Please type fit or load.')                                                      #|
#-----------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------
def choose_previous_identification(project_dir):                                                        #|
    """Checks the previous results folder before accepting it."""                                       #|
    print('Enter a saved selection JSON or its results folder. Older workflow run folders also work.')  #|
    while True:                                                                                         #|
        folder = input('Previous results folder: ').strip().strip('"').strip("'")                       #|
        if not folder:                                                                                  #|
            print('Please enter the previous results folder.')                                          #|
            continue                                                                                    #|
        try:                                                                                            #|
            previous = read_previous_identification(project_dir, Path(folder).expanduser())             #|
        except (OSError, ValueError, KeyError) as error:                                                #|
            print('Cannot use that folder: ' + str(error))                                              #|
            continue                                                                                    #|
        print('Previous model: ' + previous['selected']['model_name'])                                  #|
        return previous                                                                                 #|
#--------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
def choose_prediction(project_dir):                                                            #|
    """Pairs a metric file with a run JSON that already points to that metric."""              #|
    metrics_dir = Path(project_dir) / 'metrics'                                                #|
    metrics = sorted(                                                                          #|
        (path for path in metrics_dir.glob('*.py')                                             #|
         if not path.name.startswith('_') and path.name != 'MetricTemplate.py'),               #|
        key=lambda path: path.name.lower(),                                                    #|
    )                                                                                          #|
    print('\nAvailable performance metrics:')                                                  #|
    metric_file = choose_file('Metric filename', metrics, metrics_dir / 'amplification.py')    #|
    settings_files = []                                                                        #|
    for path in sorted((Path(project_dir) / 'configs').glob('*.json')):                        #|
        try:                                                                                   #|
            settings = read_json(path)                                                         #|
        except (ValueError, OSError):                                                          #|
            continue                                                                           #|
        metric = settings.get('metric', {}) if isinstance(settings, dict) else {}              #|
        if isinstance(metric, dict) and metric.get('definition_file'):                         #|
            if project_path(project_dir, metric['definition_file']) == metric_file.resolve():  #|
                settings_files.append(path)                                                    #|
    if not settings_files:                                                                     #|
        raise ValueError(                                                                      #|
            'No JSON in configs/ points to {}. Set metric.definition_file in your '            #|
            'prediction settings to metrics/{}.'.format(metric_file.name, metric_file.name))   #|
    print('\nMatching prediction settings:')                                                   #|
    preferred = Path(project_dir) / 'configs' / (metric_file.stem + '_prediction.json')        #|
    settings_file = choose_file('Prediction settings filename', settings_files, preferred)     #|
    return settings_file                                                                       #|
#-----------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
# 3) Prepares the settings and shared output folders.                                            #|
                                                                                                 #|
def prepare_fit_config(project_dir, config):                                                     #|
    settings = copy.deepcopy(config)                                                             #|
    # Uses the complete dataset; the old manual LOO split is not reused here.                    #|
    settings['fit_conditions'] = list(settings['all_conditions'])                                #|
    if settings.get('selection_metric') == 'held_out_normalized_rmse':                           #|
        settings['selection_metric'] = 'global_normalized_rmse'                                  #|
    for key in ('excel_file', 'models_dir'):                                                     #|
        settings[key] = str(project_path(project_dir, settings[key]))                            #|
    for key, default in (('results_dir', 'results'), ('figures_dir', 'figures'),                 #|
                         ('generated_models_dir', 'generated_models')):                          #|
        settings[key] = str(project_path(project_dir, settings.get(key, default)))               #|
    # An optional name such as experiment_A_ keeps a separate set of files.                      #|
    prefix = settings.get('output_prefix', '')                                                   #|
    if not isinstance(prefix, str) or any(                                                       #|
        ord(character) < 32 or character in '<>:"/\\|?*' for character in prefix                 #|
    ):                                                                                           #|
        raise ValueError('output_prefix must be a filename prefix, for example experiment_A_.')  #|
    settings['output_prefix'] = prefix                                                           #|
    settings['save_model_selection'] = True                                                      #|
    return settings                                                                              #|
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------
def prepare_prediction_config(project_dir, config, selected_file, fit_settings):                                   #|
    settings = copy.deepcopy(config)                                                                               #|
    settings['selected_model_file'] = str(selected_file)                                                           #|
    settings['metric']['definition_file'] = str(project_path(project_dir, settings['metric']['definition_file']))  #|
    settings['model_interfaces_file'] = str(project_path(project_dir, settings['model_interfaces_file']))          #|
    for definition in settings.get('property_definitions', {}).values():                                           #|
        definition['definition_file'] = str(project_path(project_dir, definition['definition_file']))              #|
    # Uses the same folders and optional prefix as fitting and LOO.                                                #|
    for key in ('results_dir', 'figures_dir', 'generated_models_dir', 'output_prefix'):                            #|
        settings[key] = fit_settings[key]                                                                          #|
    if settings.get('output_file'):                                                                                #|
        filename = Path(settings['output_file']).name                                                              #|
        settings['output_file'] = str(Path(settings['figures_dir']) / filename)                                    #|
    return settings                                                                                                #|
#-------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------
# 4) Runs the steps in order and records which ones completed.                                                       #|
                                                                                                                     #|
def execute_workflow(project_dir, fit_file, config, loo_conditions, prediction_file,                                 #|
                     mode='fit', previous_identification=None):                                                      #|
    """Fits or loads an identification, then runs fresh LOO and/or prediction."""                                    #|
    if mode not in ('fit', 'load'):                                                                                  #|
        raise ValueError('Workflow mode must be fit or load.')                                                       #|
    if mode == 'load' and previous_identification is None:                                                           #|
        raise ValueError('Load mode needs a previous model identification.')                                         #|
    project_dir = Path(project_dir).resolve()                                                                        #|
    fit_settings = prepare_fit_config(project_dir, config)                                                           #|
    results_dir = Path(fit_settings['results_dir'])                                                                  #|
    figures_dir = Path(fit_settings['figures_dir'])                                                                  #|
    prefix = fit_settings['output_prefix']                                                                           #|
    # Protects a loaded fit if it happens to use one of the planned LOO filenames.                                   #|
    if mode == 'load':                                                                                               #|
        model_name = previous_identification['selected']['model_name']                                               #|
        for index in range(1, len(loo_conditions) + 1):                                                              #|
            fold_file = results_dir / '{}loo_{:02d}_{}_fit.npz'.format(prefix, index, model_name)                    #|
            if fold_file.resolve() == Path(previous_identification['fit_file']).resolve():                           #|
                raise ValueError('This LOO would replace the loaded fit. Choose a different output_prefix.')         #|
    results_dir.mkdir(parents=True, exist_ok=True)                                                                   #|
    figures_dir.mkdir(parents=True, exist_ok=True)                                                                   #|
    summary_file = results_dir / (prefix + 'workflow_summary.json')                                                  #|
    summary = {                                                                                                      #|
        'fitting_settings_source': str(fit_file),                                                                    #|
        'identification_mode': mode,                                                                                 #|
        'previous_selection_file': (                                                                                 #|
            str(previous_identification['selection_file']) if mode == 'load' else None),                             #|
        'prediction_settings_source': str(prediction_file) if prediction_file else None,                             #|
        'loo_conditions': list(loo_conditions),                                                                      #|
        'completed_steps': [], 'status': 'running',                                                                  #|
    }                                                                                                                #|
    save_json(summary_file, summary)                                                                                 #|
    print('\nResults: ' + str(results_dir), flush=True)                                                              #|
    print('Figures: ' + str(figures_dir), flush=True)                                                                #|
    print('Matching output filenames are replaced when a calculation is repeated.', flush=True)                      #|
    try:                                                                                                             #|
        # Captures the chosen metric settings now, before a potentially long fit.                                    #|
        prediction_config = read_json(prediction_file) if prediction_file else None                                  #|
                                                                                                                     #|
        # Fits all candidates, or reads an existing fit without copying its data.                                    #|
        if mode == 'fit':                                                                                            #|
            fit_snapshot = results_dir / (prefix + 'fitting_settings.json')                                          #|
            save_json(fit_snapshot, fit_settings)                                                                    #|
            selected_file = results_dir / (prefix + 'selected_model.json')                                           #|
            # Removes the old selection so a stopped fit cannot appear successful.                                   #|
            if selected_file.is_file():                                                                              #|
                selected_file.unlink()                                                                               #|
            print('\nRunning model identification. Fitting may take several minutes.', flush=True)                   #|
            run_script(project_dir, 'Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py',                       #|
                       fit_snapshot, results_dir / (prefix + 'identification.log'), extra_args=['--mode', 'fit'])    #|
            identification = read_previous_identification(project_dir, selected_file)                                #|
        else:                                                                                                        #|
            identification = previous_identification                                                                 #|
            selected_file = results_dir / (prefix + 'loaded_model.json')                                             #|
            # Saves only the file locations; the previous model and parameters stay in place.                        #|
            selected = dict(identification['selected'])                                                              #|
            selected['model_file'] = str(identification['model_file'])                                               #|
            selected['archived_model_file'] = str(identification['model_file'])                                      #|
            selected['fit_results_file'] = str(identification['fit_file'])                                           #|
            if identification.get('parameters_file'):                                                                #|
                selected['fitted_parameters_file'] = str(identification['parameters_file'])                          #|
            selected['loaded_from'] = str(identification['selection_file'])                                          #|
            save_json(selected_file, selected)                                                                       #|
        selected = identification['selected']                                                                        #|
        summary['selected_model_file'] = str(selected_file)                                                          #|
        summary['completed_steps'].append('model_identification' if mode == 'fit' else 'load_identification')        #|
        save_json(summary_file, summary)                                                                             #|
        message = 'Completed model identification: ' if mode == 'fit' else 'Loaded previous model identification: '  #|
        print(message + selected['model_name'], flush=True)                                                          #|
                                                                                                                     #|
        # Names the prediction files after the model and metric.                                                     #|
        prediction_snapshot = None                                                                                   #|
        if prediction_config is not None:                                                                            #|
            prediction_settings = prepare_prediction_config(                                                         #|
                project_dir, prediction_config, selected_file, fit_settings)                                         #|
            if prediction_settings.get('output_file'):                                                               #|
                prediction_stem = Path(prediction_settings['output_file']).stem                                      #|
            else:                                                                                                    #|
                metric_name = Path(prediction_settings['metric']['definition_file']).stem                            #|
                prediction_stem = prefix + selected['model_name'] + '_' + metric_name                                #|
            prediction_snapshot = results_dir / (prediction_stem + '.pending_settings.json')                         #|
                                                                                                                     #|
        # Refits only the selected model, leaving out one condition at a time.                                       #|
        if loo_conditions:                                                                                           #|
            loo_table = run_leave_one_out(                                                                           #|
                project_dir, fit_settings, loo_conditions,                                                           #|
                identification['model_file'], model_name=selected['model_name'])                                     #|
            summary['loo_summary_file'] = str(loo_table)                                                             #|
            summary['completed_steps'].append('loo')                                                                 #|
            save_json(summary_file, summary)                                                                         #|
            print('Completed LOO.', flush=True)                                                                      #|
                                                                                                                     #|
        if prediction_snapshot:                                                                                      #|
            # Keeps previous successful settings until the new prediction finishes.                                  #|
            save_json(prediction_snapshot, prediction_settings)                                                      #|
            print('\nRunning performance metric prediction.', flush=True)                                            #|
            run_script(project_dir, 'Code2_MasterSBMLPrediction_GeneralInterfaces.py',                               #|
                       prediction_snapshot, results_dir / (prediction_stem + '.log'))                                #|
            saved_settings = results_dir / (prediction_stem + '.settings.json')                                      #|
            prediction_snapshot.unlink()                                                                             #|
            summary['prediction_settings_file'] = str(saved_settings)                                                #|
            summary['completed_steps'].append('performance_prediction')                                              #|
            save_json(summary_file, summary)                                                                         #|
            print('Completed performance metric prediction.', flush=True)                                            #|
        summary['status'] = 'completed'                                                                              #|
        save_json(summary_file, summary)                                                                             #|
    except (Exception, KeyboardInterrupt) as error:                                                                  #|
        summary['status'] = 'stopped'                                                                                #|
        summary['error'] = str(error) or 'Interrupted by the user.'                                                  #|
        save_json(summary_file, summary)                                                                             #|
        raise                                                                                                        #|
    print('\nDone. Results: ' + str(results_dir), flush=True)                                                        #|
    return results_dir                                                                                               #|
#---------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------
# 5) Starts the interactive launcher when this file is run.                                                                    #|
                                                                                                                               #|
def main():                                                                                                                    #|
    project_dir = Path(__file__).resolve().parent                                                                              #|
    parser = argparse.ArgumentParser(description='Run SAMPL fitting, optional LOO, and performance prediction.')               #|
    parser.add_argument('--fit-config',                                                                                        #|
                        help='Fitting JSON path; skips the filename prompt when supplied.')                                    #|
    parser.add_argument('--check', action='store_true', help='Check data/model paths and exit without fitting.')               #|
    args = parser.parse_args()                                                                                                 #|
    # Uses the supplied path, or asks for the fitting settings filename.                                                       #|
    if args.fit_config:                                                                                                        #|
        fit_file = project_path(project_dir, args.fit_config)                                                                  #|
    else:                                                                                                                      #|
        fit_file = choose_fit_config(project_dir)                                                                              #|
    config = read_json(fit_file)                                                                                               #|
    print('Fitting settings: ' + str(fit_file))                                                                                #|
    if args.check:                                                                                                             #|
        detect_inputs(project_dir, config)                                                                                     #|
        print('Input paths and condition settings checked. No calculations were run.')                                         #|
        return                                                                                                                 #|
    mode = choose_mode()                                                                                                       #|
    previous = None                                                                                                            #|
    if mode == 'fit':                                                                                                          #|
        detect_inputs(project_dir, config)                                                                                     #|
    else:                                                                                                                      #|
        previous = choose_previous_identification(project_dir)                                                                 #|
    loo_conditions = []                                                                                                        #|
    if ask_yes_no('\nPerform LOO after model identification?'):                                                                #|
        loo_conditions = list(config['all_conditions'])                                                                        #|
        if len(loo_conditions) < 2:                                                                                            #|
            raise ValueError('LOO needs at least two conditions.')                                                             #|
        if mode == 'load':                                                                                                     #|
            print('Fresh LOO fits use the data and settings in: ' + str(fit_file))                                             #|
            detect_inputs(project_dir, config, model_file=previous['model_file'])                                              #|
        print('LOO will refit the selected model {} times, leaving out one condition each time.'.format(len(loo_conditions)))  #|
    prediction_file = None                                                                                                     #|
    if ask_yes_no('\nPerform performance metric prediction?', default=True):                                                   #|
        prediction_file = choose_prediction(project_dir)                                                                       #|
    execute_workflow(project_dir, fit_file, config, loo_conditions, prediction_file,                                           #|
                     mode=mode, previous_identification=previous)                                                              #|
#-------------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------
if __name__ == '__main__':                                                    #|
    try:                                                                      #|
        main()                                                                #|
    except (EOFError, KeyboardInterrupt):                                     #|
        print('\nStopped by the user.')                                       #|
        raise SystemExit(1)                                                   #|
    except (OSError, ValueError, KeyError, RuntimeError) as error:            #|
        print('\nStopped: ' + str(error))                                     #|
        raise SystemExit(1)                                                   #|
#------------------------------------------------------------------------------
