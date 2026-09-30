#---------------------------
# CODE 2: PERFORMANCE PREDICTION
#---------------------------
# Uses the model and fitted parameters from Code 1 to calculate
# model properties, then applies the chosen performance metric.
#
# INPUT
#   configs/<chosen_metric>.json                : chosen run settings
#   configs/model_interfaces.json                : shared model defaults
#   results/selected_model.json                  : points to model and fit files
#   metrics/<chosen_metric>.py                   : simulation requests and formula
#
# OUTPUT (uses <prefix><model>_<metric> as the filename)
#   figures/<name>.png                 : performance plot
#   results/<name>.npz                 : property and metric values
#   results/<name>.settings.json       : run settings used
#   results/<name>.model_interface.json: model defaults and run choices used
#   results/<name>*.csv                : convergence diagnostics
#
# WHAT THE CODE DOES
#   1. Loads Code 1's fitted model and this run's JSON settings.
#   2. Calculates the requested properties (properties/model_properties.py).
#   3. Applies the performance formula (metrics/<chosen_metric>.py).
#   4. Plots and saves the results.
#
# RUN THIS FILE
#   Chooses a run with --config configs/<chosen_metric>.json
#   Change species, conditions and sweeps in the chosen run JSON.
#
# Full file map and explanations: README.md, "File map and responsibilities".
#---------------------------


#------------------------------------------------------------------------------
from __future__ import annotations                                            #|
                                                                              #|
# Import a user-selected Python metric file.                                  #|
import importlib.util                                                         #|
# Select a configuration without editing this script.                         #|
import argparse                                                               #|
# Save a readable table for steady-state runs.                                #|
import csv                                                                    #|
# Read shared model defaults and explicit run-level overrides.                #|
from code2_configuration import load_model_interface                          #|
import sys                                                                    #|
# Read JSON workflow configuration files.                                     #|
import json                                                                   #|
# Build portable paths relative to this script.                               #|
from pathlib import Path                                                      #|
                                                                              #|
# Create and save performance plots.                                          #|
import matplotlib.pyplot as plt                                               #|
# Perform numerical array and grid calculations.                              #|
import numpy as np                                                            #|
                                                                              #|
# Simulate the selected SBML model with the solve_ivp adapter.                #|
from sbmltoodepy_solveivp_adapter import SBMLtoODEpySolveIVPModel             #|
# Built-in numerical steady-state checks, separate from the metric formula.   #|
from properties.steady_state_simulation import read_steady_state_settings     #|
from properties.model_properties import (                                     #|
    read_property_requests, validate_request_targets,                         #|
    load_property_functions, collect_properties, LIMIT_MODES,                 #|
)                                                                             #|
# Reads the sampling-time grid for time-course calculations.                  #|
from properties.time_course_simulation import read_time_values                #|
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
# 0) PATHS AND CONFIGURATION                                                                   #|
# 0a) Finds and reads this run's JSON settings.                                                #|
                                                                                               #|
# Finds this script's full path, then keeps its containing folder.                             #|
BASE_DIR = Path(__file__).resolve().parent                                                     #|
# Sets up the options that can be supplied when starting this script.                          #|
parser = argparse.ArgumentParser(description="Predict performance from a fitted SBML model.")  #|
# Requires --config so the chosen analysis is explicit when starting Code 2.                   #|
parser.add_argument(                                                                           #|
    "--config", required=True,                                                                 #|
    help="Absolute or project-relative JSON configuration path.",                              #|
)                                                                                              #|
# Reads the supplied options, using defaults for any that were omitted.                        #|
args = parser.parse_args()                                                                     #|
# Converts the chosen JSON filename into a path that Python can work with.                     #|
CONFIG_FILE = Path(args.config)                                                                #|
# Checks whether the filename needs the project folder added in front.                         #|
if not CONFIG_FILE.is_absolute():                                                              #|
    # Joins the folder and filename; '/' joins paths here, rather than dividing.               #|
    CONFIG_FILE = BASE_DIR / CONFIG_FILE                                                       #|
                                                                                               #|
# Opens the JSON for reading and closes it automatically when finished.                        #|
with open(CONFIG_FILE, "r", encoding="utf-8") as f:                                            #|
    # Stores its settings as a Python dictionary: names paired with values.                    #|
    config = json.load(f)                                                                      #|
#-----------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# 0b) Defines two helpers for later use: file paths and sweep values.            #|
                                                                                 #|
# Defines a reusable path helper; its calculation runs when it is called later.  #|
def resolve_path(path_value):                                                    #|
    """Resolve absolute or project-relative paths."""                            #|
    # Converts a filename from the JSON into a Python path.                      #|
    path = Path(path_value)                                                      #|
    # Keeps a full path, or adds the project folder to a relative path.          #|
    return path if path.is_absolute() else BASE_DIR / path                       #|
#---------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# Defines a helper that turns a supplied list or range into numerical values.               #|
def parse_values(value_spec):                                                               #|
    """Accept explicit values or a linear/logarithmic start/stop/points range."""           #|
    # A dictionary describes a range; a list already contains the chosen values.            #|
    if isinstance(value_spec, dict):                                                        #|
        # Reads the first value in the sweep as a number.                                   #|
        start = float(value_spec["start"])                                                  #|
        # Reads the last value in the sweep as a number.                                    #|
        stop = float(value_spec["stop"])                                                    #|
        # Reads how many values the sweep should contain.                                   #|
        points = value_spec["points"]                                                       #|
        # Rejects True/False, fractional counts and counts smaller than one.                #|
        if isinstance(points, bool) or int(points) != points or points < 1:                 #|
            raise ValueError("Range points must be a positive integer.")                    #|
        # Uses linear spacing if the JSON does not specify a spacing choice.                #|
        spacing = value_spec.get("spacing", "linear")                                       #|
        if spacing == "linear":                                                             #|
            # Creates values separated by equal differences, e.g. 1, 2, 3.                  #|
            return np.linspace(start, stop, int(points))                                    #|
        if spacing == "log":                                                                #|
            # Logarithmic spacing requires both endpoints to be above zero.                 #|
            if start <= 0 or stop <= 0:                                                     #|
                raise ValueError("Logarithmic control ranges require positive endpoints.")  #|
            # Creates values separated by equal ratios, e.g. 1, 10, 100.                    #|
            return np.geomspace(start, stop, int(points))                                   #|
        # Stops if the spacing choice is neither 'linear' nor 'log'.                        #|
        raise ValueError("Range spacing must be 'linear' or 'log'.")                        #|
    # Converts an explicit list, e.g. [1, 3, 10], into a numerical array.                   #|
    return np.asarray(value_spec, dtype=float)                                              #|
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# 0c) Reads which metric to use and the settings for its formula.                     #|
                                                                                      #|
# Reads the 'metric' block from the run JSON.                                         #|
METRIC_CONFIG = config["metric"]                                                      #|
# Reads the metric's name for reporting.                                              #|
METRIC_NAME = str(METRIC_CONFIG["name"])                                              #|
# Finds the Python file containing the metric's formula.                              #|
METRIC_FILE = resolve_path(METRIC_CONFIG["definition_file"])                          #|
# Reads formula-specific settings; {} means no settings were supplied.                #|
METRIC_SETTINGS = METRIC_CONFIG.get("settings", {})                                   #|
                                                                                      #|
# Checks that these settings are named values, rather than a list or number.          #|
if not isinstance(METRIC_SETTINGS, dict):                                             #|
    raise ValueError("metric.settings must be a JSON object.")                        #|
                                                                                      #|
# Checks that the selected metric file actually exists.                               #|
if not METRIC_FILE.is_file():                                                         #|
    raise FileNotFoundError(f"Metric definition file was not found:\n{METRIC_FILE}")  #|
#--------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# 0d) Loads that Python file so its functions can be used below.                         #|
                                                                                         #|
# Prepares instructions for importing the metric from its chosen file path.              #|
metric_spec = importlib.util.spec_from_file_location("performance_metric", METRIC_FILE)  #|
# Stops if Python cannot prepare a way to load that file.                                #|
if metric_spec is None or metric_spec.loader is None:                                    #|
    raise ImportError(f"Could not import metric file:\n{METRIC_FILE}")                   #|
                                                                                         #|
# Creates a module: a container for the metric file's functions and variables.           #|
metric_module = importlib.util.module_from_spec(metric_spec)                             #|
# Registers that container with Python under the name 'performance_metric'.              #|
sys.modules["performance_metric"] = metric_module                                        #|
# Runs the file's top-level code, making its functions and settings available.           #|
metric_spec.loader.exec_module(metric_module)                                            #|
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
# 0e) Reads optional property files and locates the metric function.                           #|
                                                                                               #|
# Reads any additional property definitions; the metric will request what it needs.            #|
PROPERTY_DEFINITIONS = config.get("property_definitions", {})                                  #|
                                                                                               #|
# Finds the metric's calculate_metric function; None means it was not found.                   #|
calculate_metric = getattr(metric_module, "calculate_metric", None)                            #|
# Checks that it is something Python can call as a function.                                   #|
if not callable(calculate_metric):                                                             #|
    raise ValueError(                                                                          #|
        f"{METRIC_FILE.name} must define calculate_metric(properties, experiment, settings)."  #|
    )                                                                                          #|
                                                                                               #|
# Reports the chosen metric's name and the Python file it came from.                           #|
print(f"Loaded metric         : {METRIC_NAME}")                                                #|
print(f"Metric definition     : {METRIC_FILE}")                                                #|
#-----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------
# 1) LOAD SELECTED MODEL AND BEST-FIT PARAMETERS                                                                            #|
                                                                                                                            #|
# Locate the model-selection record produced by Part 1.                                                                     #|
SELECTED_MODEL_FILE = resolve_path(config["selected_model_file"])                                                           #|
                                                                                                                            #|
# Require the model-selection record produced by Part 1.                                                                    #|
if not SELECTED_MODEL_FILE.exists():                                                                                        #|
    # Explain how to generate the missing selection record.                                                                 #|
    raise FileNotFoundError(                                                                                                #|
        "Part 1 selection record was not found:\n"                                                                          #|
        f"{SELECTED_MODEL_FILE}\n\n"                                                                                        #|
        "Run Code1_MasterSBMLTimecourseFitting first."                                                                      #|
    )                                                                                                                       #|
                                                                                                                            #|
# Load the selected mechanism and its fitted-parameter files.                                                               #|
with open(SELECTED_MODEL_FILE, "r", encoding="utf-8") as f:                                                                 #|
    # Convert the selection record into a Python dictionary.                                                                #|
    selected = json.load(f)                                                                                                 #|
                                                                                                                            #|
# Get the selected mechanism name for reporting and plot labels.                                                            #|
MODEL_NAME = selected["model_name"]                                                                                         #|
# Resolve the selected SBML model path.                                                                                     #|
MODEL_FILE = resolve_path(selected["model_file"])                                                                           #|
# Uses the saved SBML copy when available, so later model edits do not change this fit.                                     #|
if selected.get("archived_model_file"):                                                                                     #|
    archived_model = resolve_path(selected["archived_model_file"])                                                          #|
    if archived_model.is_file():                                                                                            #|
        MODEL_FILE = archived_model                                                                                         #|
# Resolve the saved best-fit parameter path.                                                                                #|
FIT_RESULTS_FILE = resolve_path(selected["fit_results_file"])                                                               #|
                                                                                                                            #|
# Matches a copied results folder to its saved model and fitted parameters.                                                 #|
if SELECTED_MODEL_FILE.name.endswith("selected_model.json"):                                                                #|
    local_archive = SELECTED_MODEL_FILE.with_suffix(".xml")                                                                 #|
    if local_archive.is_file():                                                                                             #|
        MODEL_FILE = local_archive                                                                                          #|
    local_fit = SELECTED_MODEL_FILE.parent / Path(selected["fit_results_file"]).name                                        #|
    if local_fit.is_file():                                                                                                 #|
        FIT_RESULTS_FILE = local_fit                                                                                        #|
                                                                                                                            #|
# A load pointer must keep using its chosen source files.                                                                   #|
if selected.get("loaded_from") and (not MODEL_FILE.is_file() or not FIT_RESULTS_FILE.is_file()):                            #|
    raise FileNotFoundError("A source file used by loaded_model.json is missing. Restore it or choose another saved fit.")  #|
                                                                                                                            #|
# Accepts the original project layout when a recorded path is unavailable.                                                  #|
if not MODEL_FILE.exists():                                                                                                 #|
    local_model = BASE_DIR / "models" / f"{MODEL_NAME}.xml"                                                                 #|
    archived_model = BASE_DIR / "results" / "selected_model.xml"                                                            #|
    if local_model.exists():                                                                                                #|
        MODEL_FILE = local_model                                                                                            #|
    elif archived_model.exists():                                                                                           #|
        MODEL_FILE = archived_model                                                                                         #|
    else:                                                                                                                   #|
        raise FileNotFoundError(f"Selected SBML file was not found:\n{MODEL_FILE}")                                         #|
                                                                                                                            #|
# Use the local fitted results when the original folder is unavailable.                                                     #|
if not FIT_RESULTS_FILE.exists():                                                                                           #|
    if SELECTED_MODEL_FILE.name.endswith("selected_model.json") and SELECTED_MODEL_FILE.name != "selected_model.json":      #|
        raise FileNotFoundError(f"The named selection's fitted results were not found:\n{FIT_RESULTS_FILE}")                #|
    local_fit = BASE_DIR / "results" / f"{MODEL_NAME}_fit.npz"                                                              #|
    # Also accepts results saved by the older folder layout.                                                                #|
    if not local_fit.exists():                                                                                              #|
        local_fit = BASE_DIR / "results" / MODEL_NAME / f"{MODEL_NAME}_results.npz"                                         #|
    if local_fit.exists():                                                                                                  #|
        FIT_RESULTS_FILE = local_fit                                                                                        #|
    else:                                                                                                                   #|
        raise FileNotFoundError(                                                                                            #|
            f"Selected fit-results file was not found:\n{FIT_RESULTS_FILE}\n"                                               #|
            f"Local fallback was also not found:\n{local_fit}"                                                              #|
        )                                                                                                                   #|
                                                                                                                            #|
# Restore the best-fit parameters from Part 1.                                                                              #|
saved_fit = np.load(FIT_RESULTS_FILE, allow_pickle=True)                                                                    #|
                                                                                                                            #|
# Extract the fitted values stored in log10 space.                                                                          #|
best_x = np.asarray(saved_fit["best_x"], dtype=float)                                                                       #|
# Extract the parameter names in their fitted order.                                                                        #|
fit_params = [str(name) for name in saved_fit["fit_params"]]                                                                #|
                                                                                                                            #|
# Part 1 stores fitted kinetic parameters in log10 space.                                                                   #|
# Convert fitted log-space parameters back to physical values.                                                              #|
base_params = {                                                                                                             #|
    parameter_name: float(10 ** best_x[index])                                                                              #|
    for index, parameter_name in enumerate(fit_params)                                                                      #|
}                                                                                                                           #|
                                                                                                                            #|
# Locate the folder used for generated simulator modules.                                                                   #|
GENERATED_MODELS_DIR = resolve_path(config.get("generated_models_dir", "generated_models"))                                 #|
                                                                                                                            #|
# Initialize a simulator for the selected SBML mechanism.                                                                   #|
rr = SBMLtoODEpySolveIVPModel(sbml_file=MODEL_FILE, generated_models_dir=GENERATED_MODELS_DIR)                              #|
                                                                                                                            #|
# Report the selected mechanism.                                                                                            #|
print(f"Loaded selected model : {MODEL_NAME}")                                                                              #|
# Report the SBML file used for simulation.                                                                                 #|
print(f"Loaded SBML           : {MODEL_FILE}")                                                                              #|
# Report the best-fit parameter file.                                                                                       #|
print(f"Loaded fit            : {FIT_RESULTS_FILE}")                                                                        #|
                                                                                                                            #|
# Introduce the fitted-parameter summary.                                                                                   #|
print("\nBest-fit parameters:")                                                                                             #|
# Print each fitted parameter in scientific notation.                                                                       #|
for parameter_name, value in base_params.items():                                                                           #|
    # Display the current parameter name and physical value.                                                                #|
    print(f"{parameter_name:<20} = {value:.8e}")                                                                            #|
#----------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
# 2) SELECT THE MODEL-SPECIFIC SIMULATION INTERFACE                                                             #|
                                                                                                                #|
# Shared defaults describe the selected model; the run chooses the analysis.                                    #|
MODEL_INTERFACE = load_model_interface(config, MODEL_NAME, BASE_DIR)                                            #|
MODEL_INTERFACES_FILE = resolve_path(config["model_interfaces_file"])                                           #|
print(f"Model interface file  : {MODEL_INTERFACES_FILE}")                                                       #|
                                                                                                                #|
# Require an analysis control for Code 2.                                                                       #|
CONTROL_TARGET = MODEL_INTERFACE.get("analysis_control")                                                        #|
                                                                                                                #|
# Stop clearly when a fitting model has no analysis control.                                                    #|
if CONTROL_TARGET is None:                                                                                      #|
    # Explain that the model may still be valid for Code 1.                                                     #|
    raise ValueError("The run configuration must specify analysis_control.")                                    #|
                                                                                                                #|
# Define the supported ways to control a simulation condition.                                                  #|
SUPPORTED_CONTROL_TYPES = {                                                                                     #|
    "species_initial_concentration",  # <--------- choice for value to enter in                                 #|
    "parameter",                      # <--------- json config file row "type:"                                 #|
}                                                                                                               #|
                                                                                                                #|
# Get whether the analysis control changes a species or parameter.                                              #|
CONTROL_TYPE = str(CONTROL_TARGET["type"])                                                                      #|
# Get the selected model's SBML ID for the analysis control.                                                    #|
CONTROL_ID = str(CONTROL_TARGET["id"])                                                                          #|
                                                                                                                #|
# Reject unsupported analysis-control types.                                                                    #|
if CONTROL_TYPE not in SUPPORTED_CONTROL_TYPES:                                                                 #|
    # Report the invalid value and the supported choices.                                                       #|
    raise ValueError(                                                                                           #|
        f"Unsupported analysis_control type {CONTROL_TYPE!r}. "                                                 #|
        f"Choose from {sorted(SUPPORTED_CONTROL_TYPES)}."                                                       #|
    )                                                                                                           #|
                                                                                                                #|
# Get the SBML species used as the measured output.                                                             #|
OBSERVABLE_ID = str(MODEL_INTERFACE["observable_id"])                                                           #|
# This is a display label; no concentration-unit conversion is performed.                                       #|
OBSERVABLE_UNITS = str(MODEL_INTERFACE.get("observable_units", "model concentration units"))                    #|
                                                                                                                #|
# Use shared initial concentrations, with explicit run overrides already applied.                               #|
FIXED_INITIAL_SPECIES = {                                                                                       #|
    str(species_id): float(value)                                                                               #|
    for species_id, value in MODEL_INTERFACE.get(                                                               #|
        "fixed_initial_species",                                                                                #|
        {},                                                                                                     #|
    ).items()                                                                                                   #|
}                                                                                                               #|
                                                                                                                #|
# Use shared fixed-parameter defaults, with explicit run overrides already applied.                             #|
FIXED_PARAMETERS = {                                                                                            #|
    str(parameter_id): float(value)                                                                             #|
    for parameter_id, value in MODEL_INTERFACE.get(                                                             #|
        "fixed_parameters",                                                                                     #|
        {},                                                                                                     #|
    ).items()                                                                                                   #|
}                                                                                                               #|
                                                                                                                #|
# Get the control values and display settings.                                                                  #|
CONTROL_SETTINGS = config["control_values"]                                                                     #|
# The control grid is experimental input, not a quantity obtained by solving ODEs.                              #|
CONTROL_VALUES = parse_values(CONTROL_SETTINGS["values"])                                                       #|
if CONTROL_VALUES.ndim != 1 or CONTROL_VALUES.size == 0 or not np.all(np.isfinite(CONTROL_VALUES)):             #|
    raise ValueError("control_values.values must be a nonempty one-dimensional list of finite values.")         #|
FITTED_CONTROL_VALUE = base_params.get(CONTROL_ID) if CONTROL_TYPE == "parameter" else None                     #|
                                                                                                                #|
# A general grid can contain one point or many; each metric checks any                                          #|
# additional requirements, such as positive or increasing parameter values.                                     #|
if CONTROL_SETTINGS.get("include_fitted_value", False):                                                         #|
    if FITTED_CONTROL_VALUE is None:                                                                            #|
        raise ValueError("include_fitted_value requires a fitted analysis-control parameter.")                  #|
    if not np.min(CONTROL_VALUES) <= FITTED_CONTROL_VALUE <= np.max(CONTROL_VALUES):                            #|
        raise ValueError("The control range must contain the fitted value when include_fitted_value is true.")  #|
    matches = np.isclose(CONTROL_VALUES, FITTED_CONTROL_VALUE, rtol=1e-10, atol=0.0)                            #|
    CONTROL_VALUES = np.sort(np.append(CONTROL_VALUES[~matches], FITTED_CONTROL_VALUE))                         #|
                                                                                                                #|
                                                                                                                #|
# Get the requested subset of control values to plot.                                                           #|
CONTROL_VALUES_TO_PLOT = np.asarray(                                                                            #|
    CONTROL_SETTINGS.get("values_to_plot", CONTROL_VALUES),                                                     #|
    dtype=float,                                                                                                #|
)                                                                                                               #|
                                                                                                                #|
# Get the display label for the analysis control.                                                               #|
CONTROL_LABEL = str(CONTROL_SETTINGS.get("label", CONTROL_ID))                                                  #|
# Get the display units for the analysis control.                                                               #|
CONTROL_UNITS = str(CONTROL_SETTINGS.get("units", ""))                                                          #|
#----------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------
# 3) TIME AND NUMERICAL SETTINGS                                                                        #|
                                                                                                        #|
# Explicit time_values_min can request one or several times. Otherwise use                              #|
# the existing time_end_min/time_points grid. Limit-only metrics need no grid.                          #|
t_minutes = read_time_values(config)                                                                    #|
T_END_MIN = float(t_minutes[-1]) if t_minutes.size else 0.0                                             #|
MODEL_TIME_UNITS_PER_MINUTE = float(config.get("model_time_units_per_minute", 60.0))                    #|
if not np.isfinite(MODEL_TIME_UNITS_PER_MINUTE) or MODEL_TIME_UNITS_PER_MINUTE <= 0:                    #|
    raise ValueError("model_time_units_per_minute must be positive.")                                   #|
t_model = t_minutes * MODEL_TIME_UNITS_PER_MINUTE                                                       #|
if not np.all(np.isfinite(t_model)):                                                                    #|
    raise ValueError("Requested times overflow after conversion to model units.")                       #|
# Sets relative solver accuracy; this applies to every metric.                                          #|
RTOL = float(config.get("rtol", 1e-6))                                                                  #|
# Sets absolute solver accuracy for small state values; this is not a metric threshold.                 #|
ATOL = float(config.get("atol", 1e-6))                                                                  #|
if not np.isfinite(RTOL) or RTOL <= 0 or not np.isfinite(ATOL) or ATOL <= 0:                            #|
    raise ValueError("rtol and atol must be finite positive values.")                                   #|
if CONTROL_VALUES_TO_PLOT.ndim != 1 or not np.all(np.isfinite(CONTROL_VALUES_TO_PLOT)):                 #|
    raise ValueError("control_values.values_to_plot must be a one-dimensional list of finite values.")  #|
                                                                                                        #|
# Both limit properties use these numerical settings, but each request names                            #|
# its own scientific convergence criterion. A metric may request both.                                  #|
LIMIT_SETTINGS = config.get("limit_settings", {})                                                       #|
if not isinstance(LIMIT_SETTINGS, dict):                                                                #|
    raise ValueError("limit_settings must be a JSON object.")                                           #|
                                                                                                        #|
# Get the optional secondary parameter to sweep.                                                        #|
SWEEP_PARAM = config.get("sweep_parameter")                                                             #|
# Get the optional secondary sweep values.                                                              #|
SWEEP_VALUES = config.get("sweep_values")                                                               #|
                                                                                                        #|
# Convert configured sweep values to a numeric array when provided.                                     #|
if SWEEP_VALUES is not None:                                                                            #|
    # Store sweep values as floating-point numbers.                                                     #|
    SWEEP_VALUES = np.asarray(SWEEP_VALUES, dtype=float)                                                #|
#--------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# 4) VERIFY REQUIRED SBML IDS                                                                        #|
                                                                                                     #|
# Collect the species IDs available in the selected model.                                           #|
species_ids = set(rr.species_ids)                                                                    #|
# Collect the global parameter IDs available in the selected model.                                  #|
parameter_ids = set(rr.parameter_ids)                                                                #|
                                                                                                     #|
# Start the required-species set with all fixed initial species.                                     #|
required_species = set(FIXED_INITIAL_SPECIES)                                                        #|
                                                                                                     #|
# Require the analysis-control species when the analysis control changes its initial value.          #|
if CONTROL_TYPE == "species_initial_concentration":                                                  #|
    # Add the analysis-control species to the validation set.                                        #|
    required_species.add(CONTROL_ID)                                                                 #|
                                                                                                     #|
# Identify required species missing from the selected model.                                         #|
missing_species = sorted(required_species.difference(species_ids))                                   #|
                                                                                                     #|
# Report all missing required species in one error.                                                  #|
if missing_species:                                                                                  #|
    # List every missing species in the error message.                                               #|
    raise ValueError(                                                                                #|
        "The selected SBML model is missing required species:\n"                                     #|
        + "\n".join(f"  - {species_id}" for species_id in missing_species)                           #|
    )                                                                                                #|
                                                                                                     #|
# Require the configured observable species.                                                         #|
if OBSERVABLE_ID not in species_ids:                                                                 #|
    # Stop when the model cannot provide the requested output.                                       #|
    raise ValueError(f"Observable {OBSERVABLE_ID!r} was not found in the selected SBML model.")      #|
                                                                                                     #|
# Require the analysis-control parameter when the control changes a parameter.                       #|
if CONTROL_TYPE == "parameter" and CONTROL_ID not in parameter_ids:                                  #|
    # Report the missing analysis-control parameter and available alternatives.                      #|
    raise ValueError(                                                                                #|
        f"Analysis-control parameter {CONTROL_ID!r} was not found.\n"                                #|
        f"Available global parameters: {sorted(parameter_ids)}"                                      #|
    )                                                                                                #|
                                                                                                     #|
# Identify fixed parameters missing from the selected model.                                         #|
missing_fixed_parameters = sorted(set(FIXED_PARAMETERS).difference(parameter_ids))                   #|
                                                                                                     #|
# Report all missing fixed parameters in one error.                                                  #|
if missing_fixed_parameters:                                                                         #|
    # List every missing parameter in the error message.                                             #|
    raise ValueError(                                                                                #|
        "The selected SBML model is missing fixed parameters:\n"                                     #|
        + "\n".join(f"  - {parameter_id}" for parameter_id in missing_fixed_parameters)              #|
    )                                                                                                #|
                                                                                                     #|
# Records whether the main control varies a parameter.                                               #|
controlled_parameters = {CONTROL_ID} if CONTROL_TYPE == "parameter" else set()                       #|
                                                                                                     #|
# Prevent fixed values from silently overriding a controlled parameter.                              #|
conflicting_fixed_parameters = sorted(controlled_parameters.intersection(FIXED_PARAMETERS))          #|
                                                                                                     #|
# Report interface conflicts before any simulation begins.                                           #|
if conflicting_fixed_parameters:                                                                     #|
    # List every parameter assigned two different roles.                                             #|
    raise ValueError(                                                                                #|
        "Controlled parameters cannot also appear in fixed_parameters:\n"                            #|
        + "\n".join(f"  - {parameter_id}" for parameter_id in conflicting_fixed_parameters)          #|
    )                                                                                                #|
                                                                                                     #|
# Require the optional secondary sweep parameter when configured.                                    #|
if SWEEP_PARAM is not None and SWEEP_PARAM not in parameter_ids:                                     #|
    # Report the missing sweep parameter and available alternatives.                                 #|
    raise ValueError(                                                                                #|
        f"Sweep parameter {SWEEP_PARAM!r} was not found.\n"                                          #|
        f"Available global parameters: {sorted(parameter_ids)}"                                      #|
    )                                                                                                #|
                                                                                                     #|
                                                                                                     #|
# A secondary sweep must not also be fixed or used as the main control.                              #|
if SWEEP_PARAM is not None:                                                                          #|
    if SWEEP_PARAM in FIXED_PARAMETERS or SWEEP_PARAM in controlled_parameters:                      #|
        raise ValueError("sweep_parameter cannot also be fixed or controlled by analysis_control.")  #|
    if SWEEP_VALUES is not None and (                                                                #|
        SWEEP_VALUES.ndim != 1 or not np.all(np.isfinite(SWEEP_VALUES))                              #|
    ):                                                                                               #|
        raise ValueError("sweep_values must be a one-dimensional list of finite values.")            #|
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------
# 5) EXPERIMENT DESCRIPTION AND METRIC RESULT CHECKS                                                                         #|
                                                                                                                             #|
# A plain dictionary describes the numerical experiment. The parameter sets                                                  #|
# below add each panel's fitted/swept values. All species still evolve by SBML.                                              #|
experiment = {                                                                                                               #|
    "model_name": MODEL_NAME,                                                                                                #|
    "observable_id": OBSERVABLE_ID,                                                                                          #|
    "output_units": OBSERVABLE_UNITS,                                                                                        #|
    "control_values": CONTROL_VALUES.copy(),                                                                                 #|
    "control_units": CONTROL_UNITS,                                                                                          #|
    "control_label": CONTROL_LABEL,                                                                                          #|
    "fitted_control_value": FITTED_CONTROL_VALUE,                                                                            #|
    "analysis_control": CONTROL_TARGET.copy(),                                                                               #|
    # Passes the optional model input mapping to metrics that use it.                                                        #|
    "input_control": MODEL_INTERFACE.get("input_control"),                                                                   #|
    "initial_values": FIXED_INITIAL_SPECIES.copy(),                                                                          #|
    "fixed_parameters": FIXED_PARAMETERS.copy(),                                                                             #|
    "time_minutes": t_minutes.copy(),                                                                                        #|
    "time_model": t_model.copy(),                                                                                            #|
    "model_time_units_per_minute": MODEL_TIME_UNITS_PER_MINUTE,                                                              #|
    "rtol": RTOL,                                                                                                            #|
    "atol": ATOL,                                                                                                            #|
    "limit_settings": dict(LIMIT_SETTINGS),                                                                                  #|
    "parameters": dict(base_params, **FIXED_PARAMETERS),                                                                     #|
}                                                                                                                            #|
                                                                                                                             #|
                                                                                                                             #|
# Lets the metric describe its simulations as a dictionary of named requests.                                                #|
PROPERTY_REQUESTS = read_property_requests(metric_module, experiment, METRIC_SETTINGS, PROPERTY_DEFINITIONS)                 #|
# Checks only the model targets actually used by those requests.                                                             #|
validate_request_targets(PROPERTY_REQUESTS, experiment, species_ids, rr.state_species_ids, parameter_ids, SWEEP_PARAM)       #|
# Loads each requested property calculation, regardless of which metric requested it.                                        #|
REQUIRED_PROPERTIES = list(dict.fromkeys(request["property"] for request in PROPERTY_REQUESTS.values()))                     #|
PROPERTY_FUNCTIONS = load_property_functions(REQUIRED_PROPERTIES, PROPERTY_DEFINITIONS, BASE_DIR)                            #|
TIME_REQUESTS = [name for name, request in PROPERTY_REQUESTS.items() if request["property"] in ("concentrations", "rates")]  #|
LIMIT_REQUESTS = [name for name, request in PROPERTY_REQUESTS.items() if request["property"] in LIMIT_MODES]                 #|
if TIME_REQUESTS and t_minutes.size == 0:                                                                                    #|
    raise ValueError("Requested concentrations/rates need time_values_min or time_end_min and time_points.")                 #|
for name in LIMIT_REQUESTS:                                                                                                  #|
    property_name = PROPERTY_REQUESTS[name]["property"]                                                                      #|
    read_steady_state_settings(dict(LIMIT_SETTINGS, convergence_mode=LIMIT_MODES[property_name]))                            #|
print("Requested calculations: " + ", ".join(PROPERTY_REQUESTS))                                                             #|
#-----------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# Catches mistakes in a new metric before attempting to plot it.                             #|
def validate_metric_result(result):                                                          #|
    """Check the small result dictionary shared by all metric files."""                      #|
    if not isinstance(result, dict):                                                         #|
        raise ValueError("calculate_metric must return a dictionary.")                       #|
                                                                                             #|
    required = {"kind", "values", "label", "units"}                                          #|
    missing = required.difference(result)                                                    #|
    if missing:                                                                              #|
        raise ValueError(f"Metric result is missing: {sorted(missing)}")                     #|
                                                                                             #|
    kind = result["kind"]                                                                    #|
    if kind == "time_series":                                                                #|
        if t_minutes.size == 0:                                                              #|
            raise ValueError("A time_series metric needs a configured sampling-time grid.")  #|
        expected_shape = (len(CONTROL_VALUES), len(t_model))                                 #|
        reference_shape = (len(t_model),)                                                    #|
    elif kind == "per_control":                                                              #|
        expected_shape = (len(CONTROL_VALUES),)                                              #|
        reference_shape = ()                                                                 #|
    else:                                                                                    #|
        raise ValueError("Metric kind must be 'time_series' or 'per_control'.")              #|
                                                                                             #|
    result["values"] = np.asarray(result["values"], dtype=float)                             #|
    if result["values"].shape != expected_shape:                                             #|
        raise ValueError(                                                                    #|
            f"Metric values have shape {result['values'].shape}; "                           #|
            f"expected {expected_shape} for {kind}."                                         #|
        )                                                                                    #|
    if np.any(np.isinf(result["values"])):                                                   #|
        raise ValueError("Metric values contain infinity. Use NaN for undefined values.")    #|
                                                                                             #|
    if not isinstance(result["label"], str) or not isinstance(result["units"], str):         #|
        raise ValueError("Metric label and units must be strings.")                          #|
                                                                                             #|
    reference_values = result.get("reference_values")                                        #|
    if reference_values is not None:                                                         #|
        reference_values = np.asarray(reference_values, dtype=float)                         #|
        if reference_values.shape != reference_shape:                                        #|
            raise ValueError(                                                                #|
                f"Metric reference has shape {reference_values.shape}; "                     #|
                f"expected {reference_shape}."                                               #|
            )                                                                                #|
        if np.any(np.isinf(reference_values)):                                               #|
            raise ValueError("Metric reference contains infinity. Use NaN instead.")         #|
        result["reference_values"] = reference_values                                        #|
                                                                                             #|
    # Optional labels and reference locations come from the metric itself.                   #|
    for key in ("condition_label", "reference_label"):                                       #|
        if key in result and not isinstance(result[key], str):                               #|
            raise ValueError(f"Metric {key} must be a string.")                              #|
    if result.get("reference_control_value") is not None:                                    #|
        result["reference_control_value"] = float(result["reference_control_value"])         #|
        if not np.isfinite(result["reference_control_value"]):                               #|
            raise ValueError("Metric reference_control_value must be finite.")               #|
    return result                                                                            #|
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# 6) RUN SIMULATIONS AND CALCULATE THE SELECTED METRIC                                       #|
                                                                                             #|
# Keep one parameter dictionary and panel label per sweep condition.                         #|
parameter_sets = []                                                                          #|
row_labels = []                                                                              #|
                                                                                             #|
if SWEEP_PARAM is None:                                                                      #|
    parameter_sets.append(base_params.copy())                                                #|
    row_labels.append("Best-fit parameters")                                                 #|
else:                                                                                        #|
    if SWEEP_VALUES is None or SWEEP_VALUES.size == 0:                                       #|
        raise ValueError("sweep_parameter is set but sweep_values is empty or null.")        #|
                                                                                             #|
    for value in SWEEP_VALUES:                                                               #|
        params = base_params.copy()                                                          #|
        params[SWEEP_PARAM] = float(value)                                                   #|
        parameter_sets.append(params)                                                        #|
        row_labels.append(f"{SWEEP_PARAM} = {value:.2e}")                                    #|
                                                                                             #|
# Keep the requested properties available for inspection and saving.                         #|
properties_all = []                                                                          #|
metric_results = []                                                                          #|
                                                                                             #|
for params, row_label in zip(parameter_sets, row_labels):                                    #|
    print(f"Computing {row_label}")                                                          #|
    panel_experiment = dict(experiment)                                                      #|
    panel_experiment["parameters"] = dict(params, **FIXED_PARAMETERS)                        #|
    properties = collect_properties(                                                         #|
        PROPERTY_REQUESTS, rr, panel_experiment, params,                                     #|
        PROPERTY_FUNCTIONS, PROPERTY_DEFINITIONS,                                            #|
    )                                                                                        #|
    result = calculate_metric(properties, panel_experiment, METRIC_SETTINGS)                 #|
    result = validate_metric_result(result)                                                  #|
    properties_all.append(properties)                                                        #|
    metric_results.append(result)                                                            #|
                                                                                             #|
# Shared axes require the same result type and units in every sweep panel.                   #|
first_result = metric_results[0]                                                             #|
for result in metric_results[1:]:                                                            #|
    for key in ("kind", "label", "units"):                                                   #|
        if result[key] != first_result[key]:                                                 #|
            raise ValueError(f"Metric {key!r} changed between parameter-sweep conditions.")  #|
#---------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# 7) PLOT AND SAVE THE PERFORMANCE RESULTS                                    #|
                                                                              #|
ncols = len(metric_results)                                                   #|
fig, axes = plt.subplots(                                                     #|
    1, ncols,                                                                 #|
    figsize=(5.2 * ncols, 4.8),                                               #|
    sharex=True,                                                              #|
    sharey=True,                                                              #|
    constrained_layout=True,                                                  #|
)                                                                             #|
if ncols == 1:                                                                #|
    axes = np.asarray([axes])                                                 #|
                                                                              #|
# Select the nearest simulated value for each requested time-course curve.    #|
control_indices = [                                                           #|
    int(np.argmin(np.abs(CONTROL_VALUES - value)))                            #|
    for value in CONTROL_VALUES_TO_PLOT                                       #|
]                                                                             #|
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def format_control_value(value):                                              #|
    """Format the analysis-control value for a plot legend."""                #|
    if CONTROL_UNITS:                                                         #|
        return f"{value:g} {CONTROL_UNITS}"                                   #|
    return f"{value:g}"                                                       #|
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
for panel_index, (result, ax) in enumerate(zip(metric_results, axes)):                                                 #|
    metric_values = result["values"]                                                                                   #|
    reference_values = result.get("reference_values")                                                                  #|
    reference_label = result.get("reference_label", "Reference")                                                       #|
    reference_control = result.get("reference_control_value")                                                          #|
                                                                                                                       #|
    if result["kind"] == "time_series":                                                                                #|
        # The metric supplies its own reference, including any undefined gaps.                                         #|
        if reference_values is not None:                                                                               #|
            ax.plot(                                                                                                   #|
                t_minutes, reference_values,                                                                           #|
                linestyle="--", linewidth=2.5, color="black",                                                          #|
                marker="o" if t_minutes.size == 1 else None,                                                           #|
                label=reference_label,                                                                                 #|
            )                                                                                                          #|
                                                                                                                       #|
        for control_index in control_indices:                                                                          #|
            control_value = CONTROL_VALUES[control_index]                                                              #|
            if reference_values is not None and reference_control is not None and control_value == reference_control:  #|
                continue                                                                                               #|
            ax.plot(                                                                                                   #|
                t_minutes, metric_values[control_index],                                                               #|
                linewidth=2.0, label=format_control_value(control_value),                                              #|
                marker="o" if t_minutes.size == 1 else None,                                                           #|
            )                                                                                                          #|
                                                                                                                       #|
        ax.set_xlabel("Time (min)")                                                                                    #|
        if T_END_MIN > 0:                                                                                              #|
            ax.set_xlim(0, T_END_MIN)                                                                                  #|
    else:                                                                                                              #|
        # Scalar metrics use the full control grid on the horizontal axis.                                             #|
        ax.plot(CONTROL_VALUES, metric_values, marker="o", linewidth=2.0)                                              #|
        if reference_values is not None:                                                                               #|
            ax.axhline(                                                                                                #|
                float(reference_values),                                                                               #|
                linestyle="--", linewidth=2.5, color="black",                                                          #|
                label=reference_label,                                                                                 #|
            )                                                                                                          #|
        control_axis_label = CONTROL_LABEL                                                                             #|
        if CONTROL_UNITS:                                                                                              #|
            control_axis_label += f" ({CONTROL_UNITS})"                                                                #|
        ax.set_xlabel(control_axis_label)                                                                              #|
        if CONTROL_VALUES.size > 1:                                                                                    #|
            ax.set_xlim(np.min(CONTROL_VALUES), np.max(CONTROL_VALUES))                                                #|
            if FITTED_CONTROL_VALUE is not None and CONTROL_VALUES[0] <= FITTED_CONTROL_VALUE <= CONTROL_VALUES[-1]:   #|
                ax.axvline(                                                                                            #|
                    FITTED_CONTROL_VALUE, linestyle=":", color="grey",                                                 #|
                    label="Code 1 fitted value",                                                                       #|
                )                                                                                                      #|
            if isinstance(CONTROL_SETTINGS["values"], dict) and CONTROL_SETTINGS["values"].get("spacing") == "log":    #|
                ax.set_xscale("log")                                                                                   #|
            if not np.any(np.isfinite(metric_values)):                                                                 #|
                ax.text(                                                                                               #|
                    0.5, 0.5, "No valid metric values\nCheck properties and metric settings",                          #|
                    transform=ax.transAxes, ha="center", va="center",                                                  #|
                )                                                                                                      #|
                                                                                                                       #|
                                                                                                                       #|
    ax.set_title(row_labels[panel_index])                                                                              #|
    ax.grid(True, alpha=0.3)                                                                                           #|
                                                                                                                       #|
# Avoid magnifying solver round-off into an apparent trend in a constant metric.                                       #|
if first_result["kind"] == "per_control":                                                                              #|
    finite_values = np.concatenate([                                                                                   #|
        result["values"][np.isfinite(result["values"])] for result in metric_results                                   #|
    ])                                                                                                                 #|
    if finite_values.size and first_result["units"] == "" and np.max(np.abs(finite_values)) < 1e-12:                   #|
        # Give tiny dimensionless values a readable near-zero scale.                                                   #|
        # This changes only the display; the saved numbers retain full precision.                                      #|
        axes[0].set_ylim(-0.05, 0.05)                                                                                  #|
    elif finite_values.size and np.allclose(finite_values, finite_values[0], rtol=1e-6, atol=0.0):                     #|
        padding = 0.05 * abs(finite_values[0]) if finite_values[0] != 0 else 0.05                                      #|
        axes[0].set_ylim(finite_values.min() - padding, finite_values.max() + padding)                                 #|
                                                                                                                       #|
metric_axis_label = first_result["label"]                                                                              #|
if first_result["units"]:                                                                                              #|
    metric_axis_label += f" ({first_result['units']})"                                                                 #|
axes[0].set_ylabel(metric_axis_label)                                                                                  #|
                                                                                                                       #|
# Scalar metrics without a reference do not need a legend.                                                             #|
handles, labels = axes[-1].get_legend_handles_labels()                                                                 #|
if handles:                                                                                                            #|
    axes[-1].legend(                                                                                                   #|
        title=CONTROL_LABEL if first_result["kind"] == "time_series" else None,                                        #|
        loc="upper right",                                                                                             #|
    )                                                                                                                  #|
                                                                                                                       #|
# Uses an optional condition description supplied by the selected metric.                                              #|
condition_text = first_result.get("condition_label", "")                                                               #|
title = f"{first_result['label']} | {MODEL_NAME}"                                                                      #|
if condition_text:                                                                                                     #|
    title += " | " + condition_text                                                                                    #|
fig.suptitle(title)                                                                                                    #|
#-----------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------
# Saves the plot in figures/ and its numbers and settings in results/.                                                               #|
FIGURES_DIR = resolve_path(config.get("figures_dir", "figures"))                                                                     #|
RESULTS_DIR = resolve_path(config.get("results_dir", "results"))                                                                     #|
OUTPUT_PREFIX = config.get("output_prefix", "")                                                                                      #|
if not isinstance(OUTPUT_PREFIX, str) or any(                                                                                        #|
    ord(character) < 32 or character in '<>:"/\\|?*' for character in OUTPUT_PREFIX                                                  #|
):                                                                                                                                   #|
    raise ValueError("output_prefix must be a filename prefix, for example experiment_A_.")                                          #|
OUTPUT_NAME = f"{OUTPUT_PREFIX}{MODEL_NAME}_{METRIC_FILE.stem}"                                                                      #|
# An explicitly supplied image filename still takes precedence.                                                                      #|
if config.get("output_file"):                                                                                                        #|
    OUTPUT_FILE = resolve_path(config["output_file"])                                                                                #|
else:                                                                                                                                #|
    OUTPUT_FILE = FIGURES_DIR / f"{OUTPUT_NAME}.png"                                                                                 #|
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)                                                                                #|
RESULTS_DIR.mkdir(parents=True, exist_ok=True)                                                                                       #|
fig.savefig(OUTPUT_FILE, dpi=400, bbox_inches="tight")                                                                               #|
print(f"Saved figure -> {OUTPUT_FILE}")                                                                                              #|
                                                                                                                                     #|
# Uses the image's filename to match it to its saved numerical results.                                                              #|
RESULTS_FILE = RESULTS_DIR / f"{OUTPUT_FILE.stem}.npz"                                                                               #|
saved_values = {                                                                                                                     #|
    "time_minutes": t_minutes,                                                                                                       #|
    "control_values": CONTROL_VALUES,                                                                                                #|
    "metric_values": np.stack([result["values"] for result in metric_results]),                                                      #|
    "metric_name": np.asarray(METRIC_NAME),                                                                                          #|
    "metric_kind": np.asarray(first_result["kind"]),                                                                                 #|
    "metric_label": np.asarray(first_result["label"]),                                                                               #|
    "metric_units": np.asarray(first_result["units"]),                                                                               #|
    "row_labels": np.asarray(row_labels),                                                                                            #|
}                                                                                                                                    #|
for index, result in enumerate(metric_results):                                                                                      #|
    if result.get("reference_values") is not None:                                                                                   #|
        saved_values[f"reference_values_{index}"] = result["reference_values"]                                                       #|
    if result.get("reference_control_value") is not None:                                                                            #|
        saved_values[f"reference_control_value_{index}"] = np.asarray(result["reference_control_value"])                             #|
# Save each requested property under its own name. Mixed requests therefore                                                          #|
# retain independent convergence flags and diagnostics. No pickle is needed.                                                         #|
saved_values["request_names"] = np.asarray(list(PROPERTY_REQUESTS))                                                                  #|
saved_values["requested_properties"] = np.asarray([request["property"] for request in PROPERTY_REQUESTS.values()])                   #|
for name in PROPERTY_REQUESTS:                                                                                                       #|
    fields = set(properties_all[0][name])                                                                                            #|
    for properties in properties_all[1:]:                                                                                            #|
        if set(properties[name]) != fields:                                                                                          #|
            raise ValueError(f"Property {name!r} changed result fields between sweep panels.")                                       #|
    for field in sorted(fields):                                                                                                     #|
        output_key = f"property_{name}_{field}"                                                                                      #|
        if output_key in saved_values:                                                                                               #|
            raise ValueError(f"Property output key {output_key!r} is ambiguous. Rename the imported property or diagnostic field.")  #|
        saved_values[output_key] = np.stack([p[name][field] for p in properties_all])                                                #|
                                                                                                                                     #|
if FITTED_CONTROL_VALUE is not None:                                                                                                 #|
    saved_values["fitted_control_value"] = np.asarray(FITTED_CONTROL_VALUE)                                                          #|
np.savez(RESULTS_FILE, **saved_values)                                                                                               #|
                                                                                                                                     #|
# One readable diagnostic table per requested limit. For a single limit the                                                          #|
# existing output.csv name is preserved. A time-series metric has no single                                                          #|
# metric value per control, so its table leaves that column empty.                                                                   #|
for name in LIMIT_REQUESTS:                                                                                                          #|
    property_name = PROPERTY_REQUESTS[name]["property"]                                                                              #|
    CSV_FILE = RESULTS_DIR / (OUTPUT_FILE.stem + (".csv" if len(LIMIT_REQUESTS) == 1 else f".{name}.csv"))                           #|
    with open(CSV_FILE, "w", newline="", encoding="utf-8") as f:                                                                     #|
        writer = csv.writer(f)                                                                                                       #|
        estimate_column = "xss" if property_name == "steady_state" else "x_plateau"                                                  #|
        writer.writerow([                                                                                                            #|
            "condition", "request", "parameter", "parameter_value", "species", estimate_column,                                      #|
            "converged", "metric_value", "final_output", "final_time_min",                                                           #|
            "max_abs_rate_per_model_time", "max_scaled_change", "worst_rate_species", "status",                                      #|
            "convergence_mode", "full_system_converged", "output_max_abs_rate_per_model_time",                                       #|
            "output_max_scaled_change", "plateau_start_time_min", "plateau_target_time_min",                                         #|
            "plateau_max_scaled_change",                                                                                             #|
        ])                                                                                                                           #|
        for row_label, properties, result in zip(row_labels, properties_all, metric_results):                                        #|
            limit = properties[name]                                                                                                 #|
            same_grid = np.array_equal(limit["control_values"], CONTROL_VALUES)                                                      #|
            for index, control_value in enumerate(limit["control_values"]):                                                          #|
                writer.writerow([                                                                                                    #|
                    row_label, name, CONTROL_ID, control_value, OBSERVABLE_ID,                                                       #|
                    limit["values"][index], limit["converged"][index],                                                               #|
                    result["values"][index] if result["kind"] == "per_control" and same_grid else "",                                #|
                    limit["final_output"][index], limit["time_minutes"][index],                                                      #|
                    limit["max_abs_rate"][index], limit["max_scaled_change"][index],                                                 #|
                    limit["worst_rate_species"][index], limit["status"][index],                                                      #|
                    LIMIT_MODES[property_name], limit["full_system_converged"][index],                                               #|
                    limit["output_max_abs_rate"][index], limit["output_max_scaled_change"][index],                                   #|
                    limit["plateau_start_time_minutes"][index], limit["plateau_target_time_minutes"][index],                         #|
                    limit["plateau_max_scaled_change"][index],                                                                       #|
                ])                                                                                                                   #|
    print(f"Saved convergence table -> {CSV_FILE}")                                                                                  #|
                                                                                                                                     #|
SETTINGS_FILE = RESULTS_DIR / f"{OUTPUT_FILE.stem}.settings.json"                                                                    #|
with open(SETTINGS_FILE, "w", encoding="utf-8") as f:                                                                                #|
    json.dump(config, f, indent=2)                                                                                                   #|
    f.write("\n")                                                                                                                    #|
                                                                                                                                     #|
print(f"Saved metric values -> {RESULTS_FILE}")                                                                                      #|
print(f"Saved run settings  -> {SETTINGS_FILE}")                                                                                     #|
                                                                                                                                     #|
# Shared defaults may change later. Save the resolved interface used for this run.                                                   #|
INTERFACE_SNAPSHOT_FILE = RESULTS_DIR / f"{OUTPUT_FILE.stem}.model_interface.json"                                                   #|
with open(INTERFACE_SNAPSHOT_FILE, "w", encoding="utf-8") as f:                                                                      #|
    json.dump({                                                                                                                      #|
        "model_name": MODEL_NAME,                                                                                                    #|
        "model_interfaces_file": str(MODEL_INTERFACES_FILE),                                                                         #|
        "resolved_interface": MODEL_INTERFACE,                                                                                       #|
        "requested_properties": REQUIRED_PROPERTIES,                                                                                 #|
        "property_requests": PROPERTY_REQUESTS,                                                                                      #|
    }, f, indent=2)                                                                                                                  #|
    f.write("\n")                                                                                                                    #|
print(f"Saved model interface -> {INTERFACE_SNAPSHOT_FILE}")                                                                         #|
                                                                                                                                     #|
plt.show()                                                                                                                           #|
#-------------------------------------------------------------------------------------------------------------------------------------
