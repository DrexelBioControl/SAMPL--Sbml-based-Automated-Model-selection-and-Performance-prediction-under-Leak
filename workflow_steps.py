#---------------------------
# WORKFLOW STEPS
#---------------------------
# Input:  Code 1/2 scripts, settings and a new or previously selected model.
# Output: Logs and LOO results in the shared results folder.
#
# Runs each script using the same Python environment as the launcher.
# In load mode, reads the existing model and fit without copying them.
# For LOO, refits the selected model after withholding one condition.
# Each fold has its own filename prefix; the identification fit stays intact.
#---------------------------


#------------------------------------------------------------------------------
import copy                                                                   #|
import csv                                                                    #|
import json                                                                   #|
import math                                                                   #|
import os                                                                     #|
from pathlib import Path                                                      #|
import subprocess                                                             #|
import sys                                                                    #|
#------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# 1) Finds the files from a previous model identification.                                   #|
                                                                                             #|
# Finds a saved file relative to its results folder or the project folder.                   #|
def find_saved_file(filename, results_dir, project_dir):                                     #|
    if not isinstance(filename, str) or not filename.strip():                                #|
        return None                                                                          #|
    path = Path(filename).expanduser()                                                       #|
    candidates = [path] if path.is_absolute() else [results_dir / path, project_dir / path]  #|
    for candidate in candidates:                                                             #|
        if candidate.is_file():                                                              #|
            return candidate.resolve()                                                       #|
    return None                                                                              #|
#---------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------
# Reads a previous selection and locates its saved model and fitted parameters.                                                #|
def read_previous_identification(project_dir, results_folder):                                                                 #|
    project_dir = Path(project_dir).resolve()                                                                                  #|
    entered_path = Path(results_folder).expanduser()                                                                           #|
    if not entered_path.is_absolute():                                                                                         #|
        entered_path = project_dir / entered_path                                                                              #|
    entered_path = entered_path.resolve()                                                                                      #|
                                                                                                                               #|
    # Accepts a results folder, a workflow run folder, or the selection JSON itself.                                           #|
    if entered_path.is_file():                                                                                                 #|
        candidates = [entered_path]                                                                                            #|
    else:                                                                                                                      #|
        candidates = [                                                                                                         #|
            entered_path / "selected_model.json",                                                                              #|
            entered_path / "identification" / "results" / "selected_model.json",                                               #|
            entered_path / "results" / "selected_model.json",                                                                  #|
        ]                                                                                                                      #|
    selection_file = next((path for path in candidates if path.is_file()), None)                                               #|
    if selection_file is None:                                                                                                 #|
        raise ValueError(f"Cannot find selected_model.json in the previous results: {entered_path}")                           #|
    try:                                                                                                                       #|
        selected = json.loads(selection_file.read_text(encoding="utf-8-sig"))                                                  #|
    except (OSError, ValueError) as error:                                                                                     #|
        raise ValueError(f"Cannot read the previous model selection: {selection_file}") from error                             #|
    if not isinstance(selected, dict):                                                                                         #|
        raise ValueError(f"The previous model selection must contain a JSON object: {selection_file}")                         #|
    for key in ("model_name", "model_file", "fit_results_file"):                                                               #|
        if not isinstance(selected.get(key), str) or not selected[key].strip():                                                #|
            raise ValueError(f"The previous model selection is missing {key}: {selection_file}")                               #|
                                                                                                                               #|
    # Keeps the original model name for its aliases and parameter settings.                                                    #|
    model_name = selected["model_name"]                                                                                        #|
    if (                                                                                                                       #|
        model_name != model_name.strip() or model_name.endswith(".")                                                           #|
        or ".." in model_name or any(character in model_name for character in '\\/:*?"<>|\x00')                                #|
    ):                                                                                                                         #|
        raise ValueError(f"The saved model_name is not a usable filename: {model_name!r}")                                     #|
    results_dir = selection_file.parent                                                                                        #|
                                                                                                                               #|
    # Matches the archive to this selection, including an optional filename prefix.                                            #|
    model_file = None                                                                                                          #|
    if selection_file.name.endswith("selected_model.json"):                                                                    #|
        archive = selection_file.with_suffix(".xml")                                                                           #|
        if archive.is_file():                                                                                                  #|
            model_file = archive                                                                                               #|
    if model_file is None:                                                                                                     #|
        model_file = find_saved_file(selected.get("archived_model_file"), results_dir, project_dir)                            #|
    if model_file is None:                                                                                                     #|
        model_file = find_saved_file(selected["model_file"], results_dir, project_dir)                                         #|
    if model_file is None:                                                                                                     #|
        raise ValueError(f"Cannot find the saved SBML model for {model_name!r}: {results_dir}")                                #|
    # Tries the recorded filename here first, then the original folder layout.                                                 #|
    prefixed_selection = selection_file.name.endswith("selected_model.json") and selection_file.name != "selected_model.json"  #|
    fit_candidates = [results_dir / Path(selected["fit_results_file"]).name]                                                   #|
    # A named selection must not fall back to an unrelated unprefixed fit.                                                     #|
    if not prefixed_selection:                                                                                                 #|
        fit_candidates.extend([                                                                                                #|
            results_dir / f"{model_name}_fit.npz",                                                                             #|
            results_dir / model_name / f"{model_name}_results.npz",                                                            #|
        ])                                                                                                                     #|
    fit_file = next((path for path in fit_candidates if path.is_file()), None)                                                 #|
    if fit_file is None:                                                                                                       #|
        fit_file = find_saved_file(selected["fit_results_file"], results_dir, project_dir)                                     #|
    # A load pointer uses the exact source fit, even if another fit has the same name here.                                    #|
    if selected.get("loaded_from"):                                                                                            #|
        fit_file = find_saved_file(selected["fit_results_file"], results_dir, project_dir)                                     #|
    if fit_file is None:                                                                                                       #|
        raise ValueError(f"Cannot find the saved fitted results for {model_name!r}: {results_dir}")                            #|
    parameter_candidates = []                                                                                                  #|
    if selected.get("fitted_parameters_file"):                                                                                 #|
        parameter_candidates.append(results_dir / Path(selected["fitted_parameters_file"]).name)                               #|
    if not prefixed_selection:                                                                                                 #|
        parameter_candidates.extend([                                                                                          #|
            results_dir / f"{model_name}_parameters.csv",                                                                      #|
            results_dir / model_name / "fitted_parameters.csv",                                                                #|
        ])                                                                                                                     #|
    parameters_file = next((path for path in parameter_candidates if path.is_file()), None)                                    #|
    if parameters_file is None:                                                                                                #|
        parameters_file = find_saved_file(selected.get("fitted_parameters_file"), results_dir, project_dir)                    #|
    if selected.get("loaded_from"):                                                                                            #|
        parameters_file = find_saved_file(selected.get("fitted_parameters_file"), results_dir, project_dir)                    #|
                                                                                                                               #|
    return {                                                                                                                   #|
        "selection_file": selection_file.resolve(),                                                                            #|
        "selected": selected,                                                                                                  #|
        "model_file": model_file.resolve(),                                                                                    #|
        "fit_file": fit_file.resolve(),                                                                                        #|
        "parameters_file": parameters_file.resolve() if parameters_file else None,                                             #|
    }                                                                                                                          #|
#-------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------
# 2) Runs a script and reads its saved results.                                                     #|
                                                                                                    #|
# Runs one workflow step and keeps its detailed messages in a log.                                  #|
def run_script(project_dir, script_name, config_file, log_file, extra_args=None):                   #|
    project_dir = Path(project_dir).resolve()                                                       #|
    config_file = (project_dir / config_file).resolve()                                             #|
    log_file = (project_dir / log_file).resolve()                                                   #|
    log_file.parent.mkdir(parents=True, exist_ok=True)                                              #|
                                                                                                    #|
    command = [sys.executable, "-u", str(project_dir / script_name), "--config", str(config_file)]  #|
    if extra_args:                                                                                  #|
        command.extend(str(value) for value in extra_args)                                          #|
                                                                                                    #|
    # Saves plots without opening windows that would pause the workflow.                            #|
    environment = os.environ.copy()                                                                 #|
    environment["MPLBACKEND"] = "Agg"                                                               #|
    environment["PYTHONIOENCODING"] = "utf-8"                                                       #|
    print(f"  Detailed progress: {log_file}", flush=True)                                           #|
    try:                                                                                            #|
        with log_file.open("w", encoding="utf-8") as log:                                           #|
            subprocess.run(                                                                         #|
                command, cwd=str(project_dir), env=environment,                                     #|
                stdout=log, stderr=subprocess.STDOUT, check=True,                                   #|
            )                                                                                       #|
    except (subprocess.CalledProcessError, OSError) as error:                                       #|
        raise RuntimeError(                                                                         #|
            f"{script_name} did not complete. See the log:\n{log_file}\n{error}"                    #|
        ) from error                                                                                #|
#----------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Reads a saved error table and explains missing or incomplete results.             #|
def read_result_rows(filename):                                                     #|
    try:                                                                            #|
        with Path(filename).open("r", encoding="utf-8-sig", newline="") as handle:  #|
            return list(csv.DictReader(handle))                                     #|
    except (OSError, csv.Error) as error:                                           #|
        raise RuntimeError(f"Cannot read workflow results: {filename}") from error  #|
#------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
# Requires a usable error value before reporting a fold as completed.              #|
def read_error(row, name, filename):                                               #|
    try:                                                                           #|
        value = float(row[name])                                                   #|
    except (KeyError, TypeError, ValueError) as error:                             #|
        raise RuntimeError(f"Missing or invalid {name} in {filename}") from error  #|
    if not math.isfinite(value) or value < 0:                                      #|
        raise RuntimeError(f"Non-finite or negative {name} in {filename}")         #|
    return value                                                                   #|
#-----------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------
# 3) Refits the selected model with one condition omitted each time.                                                                        #|
                                                                                                                                            #|
def run_leave_one_out(project_dir, fit_config, condition_names, model_file, model_name=None):                                               #|
    project_dir = Path(project_dir).resolve()                                                                                               #|
    model_file = (project_dir / model_file).resolve()                                                                                       #|
    model_name = model_name or model_file.stem                                                                                              #|
    results_dir = (project_dir / fit_config.get("results_dir", "results")).resolve()                                                        #|
    figures_dir = (project_dir / fit_config.get("figures_dir", "figures")).resolve()                                                        #|
    base_prefix = fit_config.get("output_prefix", "")                                                                                       #|
    condition_names = list(condition_names)                                                                                                 #|
                                                                                                                                            #|
    # Checks the chosen pool before creating any fold files.                                                                                #|
    if (                                                                                                                                    #|
        len(condition_names) < 2                                                                                                            #|
        or any(not isinstance(name, str) or not name for name in condition_names)                                                           #|
        or len(set(condition_names)) != len(condition_names)                                                                                #|
    ):                                                                                                                                      #|
        raise ValueError("LOO requires at least two distinct condition names.")                                                             #|
    configured_names = fit_config.get("all_conditions", [])                                                                                 #|
    condition_settings = fit_config.get("conditions", {})                                                                                   #|
    for name in condition_names:                                                                                                            #|
        if name not in configured_names or name not in condition_settings:                                                                  #|
            raise ValueError(f"LOO condition {name!r} is missing from the fitting settings.")                                               #|
    if not model_file.is_file():                                                                                                            #|
        raise ValueError(f"The selected model file does not exist: {model_file}")                                                           #|
    for label, value in (("model_name", model_name), ("output_prefix", base_prefix)):                                                       #|
        if not isinstance(value, str) or ".." in value or any(character in value for character in '\\/:*?"<>|\x00'):                        #|
            raise ValueError(f"The {label} must be a simple filename label.")                                                               #|
                                                                                                                                            #|
    # Clears the old summary; a new one is written only after all folds succeed.                                                            #|
    results_dir.mkdir(parents=True, exist_ok=True)                                                                                          #|
    summary_csv = results_dir / f"{base_prefix}{model_name}_loo_summary.csv"                                                                #|
    summary_json = summary_csv.with_suffix(".json")                                                                                         #|
    for path in (summary_csv, summary_json):                                                                                                #|
        if path.is_file():                                                                                                                  #|
            path.unlink()                                                                                                                   #|
                                                                                                                                            #|
    rows = []                                                                                                                               #|
    for index, held_out in enumerate(condition_names, start=1):                                                                             #|
        print(f"LOO fold {index}/{len(condition_names)}: withholds {held_out}.", flush=True)                                                #|
        fold_prefix = f"{base_prefix}loo_{index:02d}_"                                                                                      #|
        config_file = results_dir / f"{fold_prefix}fitting_settings.json"                                                                   #|
        errors_file = results_dir / f"{fold_prefix}{model_name}_held_out_rmse.csv"                                                          #|
        comparison_file = results_dir / f"{fold_prefix}fit_scores.csv"                                                                      #|
        fit_file = results_dir / f"{fold_prefix}{model_name}_fit.npz"                                                                       #|
                                                                                                                                            #|
        # Keeps the Excel file intact; only the fitting condition list changes.                                                             #|
        fold_config = copy.deepcopy(fit_config)                                                                                             #|
        fitted_names = [name for name in condition_names if name != held_out]                                                               #|
        fold_config["fit_conditions"] = fitted_names                                                                                        #|
        fold_config["all_conditions"] = condition_names.copy()                                                                              #|
        fold_config["results_dir"] = str(results_dir)                                                                                       #|
        fold_config["figures_dir"] = str(figures_dir)                                                                                       #|
        fold_config["output_prefix"] = fold_prefix                                                                                          #|
        fold_config["save_model_selection"] = False                                                                                         #|
        config_file.write_text(json.dumps(fold_config, indent=2) + "\n", encoding="utf-8")                                                  #|
                                                                                                                                            #|
        # Removes old fold tables so a missing new result cannot look successful.                                                           #|
        for path in (errors_file, comparison_file, fit_file):                                                                               #|
            if path.is_file():                                                                                                              #|
                path.unlink()                                                                                                               #|
                                                                                                                                            #|
        # Forces a fresh fit of this one model using only the remaining conditions.                                                         #|
        run_script(                                                                                                                         #|
            project_dir, "Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py",                                                         #|
            config_file, results_dir / f"{fold_prefix}fit.log",                                                                             #|
            extra_args=["--mode", "fit", "--model", str(model_file), "--model-name", model_name],                                           #|
        )                                                                                                                                   #|
                                                                                                                                            #|
        # Confirms that the saved prediction belongs to this fold and this model.                                                           #|
        errors = read_result_rows(errors_file)                                                                                              #|
        if len(errors) != 1 or errors[0].get("condition") != held_out:                                                                      #|
            raise RuntimeError(f"Expected only held-out condition {held_out!r} in {errors_file}")                                           #|
        comparison = read_result_rows(comparison_file)                                                                                      #|
        if len(comparison) != 1 or comparison[0].get("model") != model_name:                                                                #|
            raise RuntimeError(f"Expected only selected model {model_name!r} in {comparison_file}")                                         #|
        if not fit_file.is_file():                                                                                                          #|
            raise RuntimeError(f"Cannot find the fitted parameters from this fold: {fit_file}")                                             #|
                                                                                                                                            #|
        # Records optimizer status separately from whether the script finished.                                                             #|
        success_text = str(comparison[0].get("optimizer_success", "")).strip().lower()                                                      #|
        optimizer_success = None                                                                                                            #|
        if success_text in ("true", "1"):                                                                                                   #|
            optimizer_success = True                                                                                                        #|
        elif success_text in ("false", "0"):                                                                                                #|
            optimizer_success = False                                                                                                       #|
        elif success_text:                                                                                                                  #|
            raise RuntimeError(f"Invalid optimizer_success in {comparison_file}")                                                           #|
                                                                                                                                            #|
        row = {                                                                                                                             #|
            "fold": index,                                                                                                                  #|
            "model": model_name,                                                                                                            #|
            "held_out_condition": held_out,                                                                                                 #|
            "fit_conditions": fitted_names,                                                                                                 #|
            "raw_rmse": read_error(errors[0], "raw_rmse", errors_file),                                                                     #|
            "normalized_rmse": read_error(errors[0], "normalized_rmse", errors_file),                                                       #|
            "training_normalized_rmse": read_error(comparison[0], "global_normalized_rmse", comparison_file),                               #|
            "optimizer_success": optimizer_success,                                                                                         #|
            "config_file": str(config_file),                                                                                                #|
            "results_dir": str(results_dir),                                                                                                #|
            "fit_results_file": str(fit_file),                                                                                              #|
        }                                                                                                                                   #|
        rows.append(row)                                                                                                                    #|
        status_note = " Optimizer did not report success." if optimizer_success is False else ""                                            #|
        print(f"  Held-out nRMSE = {row['normalized_rmse']:.6g}.{status_note}", flush=True)                                                 #|
                                                                                                                                            #|
    # Saves each held-out condition's score in the summary files.                                                                           #|
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:                                                                     #|
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))                                                                           #|
        writer.writeheader()                                                                                                                #|
        for row in rows:                                                                                                                    #|
            writer.writerow(dict(row, fit_conditions=json.dumps(row["fit_conditions"])))                                                    #|
    summary = {                                                                                                                             #|
        "model_name": model_name,                                                                                                           #|
        "model_file": str(model_file),                                                                                                      #|
        "validation_scope": "Selected-model LOO: the model identity is fixed before these folds; model selection is not cross-validated.",  #|
        "condition_pool": condition_names,                                                                                                  #|
        "number_of_folds": len(rows),                                                                                                       #|
        "folds": rows,                                                                                                                      #|
    }                                                                                                                                       #|
    summary_json.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")                                        #|
    print(f"Saved LOO results: {summary_csv}", flush=True)                                                                                  #|
    return summary_csv                                                                                                                      #|
#--------------------------------------------------------------------------------------------------------------------------------------------
