# ============================================================
# MASTER SBML TIMECOURSE FITTING, VALIDATION, AND MODEL COMPARISON
# ============================================================
#
# PURPOSE
# -------
# Fit candidate SBML models to experimental time-course data.
#
# The workflow combines:
#
#     experimental time-course data
#                  +
#     candidate SBML kinetic models
#                  +
#     fitted-parameter bounds
#                  ↓
#     fitted parameters, held-out predictions, and model ranking
#
# Candidate models are loaded from SBML files in the configured models folder.
#
# CONDITIONS
# ----------
# Conditions in fit_conditions are used to estimate model parameters.
#
# Conditions in all_conditions but not fit_conditions are held out from fitting.
# Their trajectories are predicted using parameters fitted to fit_conditions.
#
# If both lists contain the same conditions, no held-out prediction is performed.
#
# This script evaluates the configured split once.
# It does not automatically repeat every possible leave-one-out split.
#
# WHAT THE SCRIPT DOES
# --------------------
# 1. Reads the shared experimental time-course data.
# 2. Builds the fitted and held-out experimental datasets.
# 3. Finds every candidate SBML model in the models folder.
# 4. Validates each model's required species and parameters.
# 5. Fits each model using only fit_conditions.
# 6. Predicts conditions excluded from fitting.
# 7. Calculates fitting and held-out prediction errors.
# 8. Saves parameters, errors, predictions, diagnostics, and figures.
# 9. Ranks the models using the configured selection metric.
# 10. Records the selected model and fitted-parameter files for Code 2.
#
# MODEL SELECTION
# ---------------
# global_normalized_rmse ranks models using fitted-condition error.
#
# held_out_normalized_rmse ranks models using held-out prediction error.
#
# FIT AND LOAD MODES
# ------------------
# MODE = "fit" estimates new parameters using the current fit_conditions.
#
# MODE = "load" reuses previously saved fitted parameters.
#
# Use fit mode whenever fit_conditions changes.
# Loading parameters fitted to different conditions would invalidate validation.
#
# USER SETTINGS
# -------------
# Change MODE to either "fit" or "load".
#
# Set CONFIG_FILE to the shared JSON configuration.
#
# Define datasets, models, parameters, bounds, solver settings,
# output locations, and the model-selection metric in the JSON file.
#
# ============================================================

# ============================================================
# 1) PACKAGES
# ============================================================

from __future__ import annotations

import json
import math
import shutil
import sys
import xml.etree.ElementTree as ET

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.optimize import least_squares
from scipy.stats import qmc


# ============================================================
# 2) EXECUTION MODE
# ============================================================
# "fit"  : perform all Sobol-start optimizations again.
# "load" : skip optimization and reuse the saved .npz file for each model.

MODE = "fit"


# ============================================================
# 3) PATHS AND CONFIGURATION
# ============================================================
# The JSON file contains all the user-specified settings for this run.
# It tells this script where the data and models live, which conditions
# are fitted, which SBML species are initialized, which parameters move,
# and what bounds and solver settings should be used.
#
# For Windows users:
# BASE_DIR makes all paths relative to this script, 
# so the full project folder can be moved to Windows from Mac without rewriting paths.

BASE_DIR = Path(__file__).resolve().parent

CONFIG_FILE = (
    BASE_DIR
    / "configs"
    / "fit_all_models.json"
)

# Open the JSON configuration file and read it into a Python dictionary.
# The dictionary values are then used in the rest of the script.

with open(CONFIG_FILE, "r", encoding="utf-8") as f:
    config = json.load(f)

RUN_NAME = config["run_name"]

MODELS_DIR = BASE_DIR / config["models_dir"]
MODEL_PATTERN = config.get("model_pattern", "*.xml")

EXCEL_FILE = BASE_DIR / config["excel_file"]
EXCEL_SHEET = config["excel_sheet"]

FIT_CONDITIONS = list(config["fit_conditions"])
ALL_CONDITIONS = list(config["all_conditions"])
CONDITION_INFO = config["conditions"]

DEFAULT_FIT_PARAMS = list(config["fit_params"])

MODEL_FIT_PARAMS = {
    str(model_name): [
        str(parameter_name)
        for parameter_name in parameter_names
    ]
    for model_name, parameter_names in config.get(
        "model_fit_params",
        {}
    ).items()
}

# Central values serve two roles:
#   1. they are the default values for parameters that are not fitted;
#   2. they provide a complete parameter dictionary before fitted values
#      overwrite the relevant entries.
ALL_PARAM_INFO = {
    p: {"central": float(v)}
    for p, v in config["all_param_central"].items()
}

CUSTOM_PARAM_BOUNDS = {
    p: tuple(map(float, bounds))
    for p, bounds in config["custom_param_bounds"].items()
}

MODEL_PARAM_CENTRAL = {
    str(model_name): {
        str(parameter_name): float(value)
        for parameter_name, value in parameter_values.items()
    }
    for model_name, parameter_values in config.get(
        "model_param_central",
        {}
    ).items()
}

MODEL_PARAM_BOUNDS = {
    str(model_name): {
        str(parameter_name): tuple(
            map(float, bounds)
        )
        for parameter_name, bounds in parameter_bounds.items()
    }
    for model_name, parameter_bounds in config.get(
        "model_param_bounds",
        {}
    ).items()
}

MODEL_FIXED_PARAMETERS = {
    str(model_name): {
        str(parameter_name): float(value)
        for parameter_name, value in parameter_values.items()
    }
    for model_name, parameter_values in config.get(
        "model_fixed_parameters",
        {}
    ).items()
}

# Most models use the shared species IDs from the configuration. Models with
# different SBML naming conventions can map each shared ID to their local ID.
MODEL_SPECIES_ALIASES = {
    str(model_name): {
        str(shared_species_id): str(model_species_id)
        for shared_species_id, model_species_id in species_aliases.items()
    }
    for model_name, species_aliases in config.get(
        "model_species_aliases",
        {},
    ).items()
}

# Most models use the shared observable. A model may provide its local output ID.
MODEL_OBSERVABLE_IDS = {
    str(model_name): str(observable_id)
    for model_name, observable_id in config.get(
        "model_observable_ids",
        {},
    ).items()
}

N_STARTS_REQUESTED = int(config["n_starts"])
RANDOM_SEED = int(config["random_seed"])
MAX_NFEV = int(config.get("max_nfev", 1000))

RESULTS_DIR = BASE_DIR / config["results_dir"]
FIGURES_DIR = BASE_DIR / config["figures_dir"]

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

OBSERVABLE_ID = str(config["observable_id"])

DATA_SIGNAL_MULTIPLIER = float(
    config.get("data_signal_multiplier", 1.0)
)

# Experimental time is stored in minutes. The SBML may use seconds.
# Multiplying experimental minutes by this number converts the experimental
# time grid into the model's time units before simulation.
MODEL_TIME_UNITS_PER_MINUTE = float(
    config.get("model_time_units_per_minute", 60.0)
)

RTOL = float(config.get("rtol", 1e-6))
ATOL = float(config.get("atol", 1e-6))

row_end_map_config = config.get("row_end_map") or {}
row_end_map = {
    str(condition_name): int(row_end)
    for condition_name, row_end in row_end_map_config.items()
}

CUTOFF_FRAC = float(config.get("cutoff_frac", 0.95))
CUTOFF_HOLD_POINTS = int(config.get("cutoff_hold_points", 20))
CUTOFF_DRIFT_TOL = float(config.get("cutoff_drift_tol", 0.03))
USE_TRANSIENT_CUTOFF = bool(config.get("use_transient_cutoff", True))


# ============================================================
# 4) SBMLtoODEpy + solve_ivp IMPORT
# ============================================================
# SBMLtoODEpy converts each SBML file into a generated Python model. 
# It needs a file named sbmltoodepy_solveivp_adapter.py in the project folder.
# This turns SBML to ODEs of the form dx/dt = f(x, p, t),
# and then uses solve_ivp to integrate the ODEs.

from sbmltoodepy_solveivp_adapter import SBMLtoODEpySolveIVPModel

GENERATED_MODELS_DIR = BASE_DIR / config.get(
    "generated_models_dir",
    "generated_models",
)

# ============================================================
# 5) FIND CANDIDATE SBML MODELS
# ============================================================
# Every matching XML file is treated as one candidate model.
# The filename stem is used as the model name. 
# Candidate files should have unique filenames even when their internal SBML IDs are equal.

MODEL_FILES = sorted(
    MODELS_DIR.glob(MODEL_PATTERN),
    key=lambda path: path.name.lower(),
)

if len(MODEL_FILES) == 0:
    raise FileNotFoundError(
        f"No SBML files matching {MODEL_PATTERN!r} were found in:\n"
        f"{MODELS_DIR}"
    )

print("RUN SUMMARY")
print("=" * 72)
print(f"Run name             : {RUN_NAME}")
print(f"Mode                 : {MODE}")
print(f"Configuration        : {CONFIG_FILE}")
print(f"Experimental file    : {EXCEL_FILE}")
print(f"Experimental sheet   : {EXCEL_SHEET}")
print(f"Observable           : {OBSERVABLE_ID}")
print(f"Fit conditions       : {FIT_CONDITIONS}")
print(f"Default fit parameters: {DEFAULT_FIT_PARAMS}")
if MODEL_FIT_PARAMS:
    print(f"Model-specific fits  : {MODEL_FIT_PARAMS}")
print(f"Requested starts     : {N_STARTS_REQUESTED}")
print(f"Candidate SBML files : {[p.name for p in MODEL_FILES]}")
print("=" * 72)


# ============================================================
# 6) LOAD EXPERIMENTAL DATA
# ============================================================
# We are assuming the data is stored in an excel workbook with one sheet containing all the conditions. 
# The spreadsheet is read without assuming a standard single-row header.
# First row lists condition names, second row lists measurement names, 
# Data starts from the third row.
# Time is assumed to be in the first column, and the rest of the columns contain measurements for each condition.

raw_data = pd.read_excel(
    EXCEL_FILE,
    sheet_name=EXCEL_SHEET,
    header=None,
    engine="openpyxl",
)

condition_header_row = int(config.get("condition_header_row", 0))
measurement_header_row = int(config.get("measurement_header_row", 1))
data_start_row = int(config.get("data_start_row", 2))
time_column_index = int(config.get("time_column_index", 0))

condition_headers = {
    str(raw_data.iloc[condition_header_row, column_index]).strip(): column_index
    for column_index in range(raw_data.shape[1])
    if (
        column_index != time_column_index
        and pd.notna(raw_data.iloc[condition_header_row, column_index])
    )
}

time_full = pd.to_numeric(
    raw_data.iloc[data_start_row:, time_column_index],
    errors="coerce",
).to_numpy(dtype=float)

valid_time = np.isfinite(time_full)
time_full = time_full[valid_time]


def read_condition_column(condition_name):
    # Purpose of function:
    # Read one experimental trajectory from Excel.
    # The data and time columns are cleaned independently of the model. 
    # The signal is then multiplied by DATA_SIGNAL_MULTIPLIER.

    if condition_name not in CONDITION_INFO:
        raise KeyError(
            f"Condition {condition_name!r} is missing from config['conditions']."
        )

    data_column_name = str(
        CONDITION_INFO[condition_name]["data_column"]
    ).strip()

    if data_column_name not in condition_headers:
        raise KeyError(
            f"Excel column {data_column_name!r} was not found.\n"
            f"Available condition headers: {list(condition_headers)}"
        )

    column_index = condition_headers[data_column_name]

    signal_full = pd.to_numeric(
        raw_data.iloc[data_start_row:, column_index],
        errors="coerce",
    ).to_numpy(dtype=float)

    signal_full = signal_full[valid_time]

    valid_signal = np.isfinite(signal_full)

    t = time_full[valid_signal]
    x = signal_full[valid_signal] * DATA_SIGNAL_MULTIPLIER

    return np.asarray(t, dtype=float), np.asarray(x, dtype=float)


# ============================================================
# 7) AUTOMATIC TRANSIENT / PLATEAU CUTOFF
# ============================================================
# The fit is intended to focus on the informative transient dynamics rather
# than allowing a long plateau to dominate the residual simply because the
# plateau contains many time points. The automatic rule estimates where the
# signal has reached and remained near its final level.

def estimate_cutoff_index(
    t,
    x,
    frac=CUTOFF_FRAC,
    hold_points=CUTOFF_HOLD_POINTS,
    drift_tol=CUTOFF_DRIFT_TOL,
):
    # Purpose of function:
    # Estimate the first point at which the trajectory has effectively settled onto its plateau.

    # A cutoff is accepted only when the final portion of the trace is itself sufficiently stable. 
    # The signal must then remain above a chosen fraction of that final level for several consecutive points. 
    # If those checks fail, the full curve is retained rather than forcing an unreliable cutoff.

    if len(x) < max(hold_points + 2, 10):
        return len(x)

    tail_start = int(0.9 * len(x))
    tail = np.asarray(x[tail_start:], dtype=float)

    if len(tail) < 2:
        return len(x)

    final_level = float(np.mean(tail))

    if (not np.isfinite(final_level)) or abs(final_level) < 1e-12:
        return len(x)

    tail_drift = float((tail[-1] - tail[0]) / final_level)

    if abs(tail_drift) > drift_tol:
        return len(x)

    threshold = frac * final_level

    for i in range(len(x) - hold_points):
        if np.all(x[i:i + hold_points] >= threshold):
            return i

    return len(x)


# ============================================================
# 8) BUILD DATASETS
# ============================================================
# Each dataset dictionary joins together the measured trajectory and the
# condition-specific SBML initial values. This is the bridge between an Excel
# column such as "25 nM I1" and the species values that SBMLtoODEpy/solve_ivp must set.

def build_dataset(condition_name):
    # Purpose of function:
    # Build a dataset dictionary for one experimental condition. 
    # The dictionary contains the time and signal arrays, the initial species values, and the plateau region that is excluded from fitting.
    # Note:
    # User can manually supply a row_end_map to override the automatic cutoff. 
    # This preserves the option to make a deliberate condition-by-condition choice when the automated plateau detector is not appropriate.

    t_all, x_all = read_condition_column(condition_name)

    if condition_name in row_end_map:
        row_end = row_end_map[condition_name]
    elif USE_TRANSIENT_CUTOFF:
        row_end = estimate_cutoff_index(t_all, x_all)
    else:
        row_end = len(x_all)

    row_end = max(2, min(int(row_end), len(x_all)))

    return {
        "name": condition_name,
        "t_exp": t_all[:row_end],
        "x_exp": x_all[:row_end],
        "t_full": t_all,
        "x_full": x_all,
        "t_plateau": t_all[row_end:],
        "x_plateau": x_all[row_end:],
        "initial_values": {
            str(species_id): float(value)
            for species_id, value in CONDITION_INFO[
                condition_name
            ]["initial_values"].items()
        },
    }


datasets = [
    build_dataset(condition_name)
    for condition_name in FIT_CONDITIONS
]

# Find every condition that is not used for parameter fitting.
HELD_OUT_CONDITIONS = [
    condition_name
    for condition_name in ALL_CONDITIONS
    if condition_name not in FIT_CONDITIONS
]

# Build prediction-only datasets for all conditions omitted from fitting.
HELD_OUT_DATASETS = [
    build_dataset(condition_name)
    for condition_name in HELD_OUT_CONDITIONS
]

print("\nDATASET SUMMARY")
print("=" * 72)

# Report all conditions used for parameter fitting.
for d in datasets:
    initial_values_text = ", ".join(
        f"{species}={value:g}"
        for species, value in d["initial_values"].items()
    )

    print(
        f"{d['name']:<22} | "
        f"type=fit      | "
        f"points={len(d['t_exp']):<4} | "
        f"{initial_values_text}"
    )

# Report all conditions reserved for prediction.
for d in HELD_OUT_DATASETS:
    initial_values_text = ", ".join(
        f"{species}={value:g}"
        for species, value in d["initial_values"].items()
    )

    print(
        f"{d['name']:<22} | "
        f"type=held-out | "
        f"points={len(d['t_exp']):<4} | "
        f"{initial_values_text}"
    )

# Report when every condition is included in fitting.
if not HELD_OUT_DATASETS:
    print("Held-out conditions: None")

# ============================================================
# 9) VISUAL CHECK OF THE FITTING WINDOW
# ============================================================
# This figure is deliberately produced before any model is fitted. 
# It lets us see exactly which experimental points contribute to the objective function
# and which late points were excluded as plateau.

n_datasets = len(datasets)
ncols = min(3, n_datasets)
nrows = int(np.ceil(n_datasets / ncols))

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(5.0 * ncols, 4.0 * nrows),
    squeeze=False,
)

for ax, d in zip(axes.flat, datasets):
    ax.plot(
        d["t_exp"],
        d["x_exp"],
        "o",
        markersize=2.5,
        label="fitted region",
    )

    if len(d["x_plateau"]) > 0:
        ax.plot(
            d["t_plateau"],
            d["x_plateau"],
            "o",
            markersize=2.5,
            alpha=0.45,
            label="excluded region",
        )

    ax.set_title(d["name"])
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Reacted reporter (nM)")
    ax.grid(True, alpha=0.3)

for ax in axes.flat[n_datasets:]:
    ax.remove()

handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper right")
fig.suptitle("Experimental fitting windows")
fig.tight_layout()

fit_window_file = (
    FIGURES_DIR
    / f"{RUN_NAME}_fitting_windows.png"
)

fig.savefig(
    fit_window_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)

print(f"\nSaved fitting-window plot: {fit_window_file}")


# ============================================================
# 10) MODEL-SPECIFIC PARAMETER SETUP
# ============================================================
# Most candidate mechanisms share a common parameter set. Future models may
# contain one or more additional leak parameters. This block starts from the
# shared settings and then applies any model-specific additions or overrides.

def get_model_fit_setup(model_name):
    # Purpose of function:
    # Assemble the parameter rules for one candidate model.

    # Shared values keep comparisons consistent. Model-specific entries allow a
    # genuine leak model to add a parameter such as klkg without forcing that
    # parameter to exist in every other SBML file.

    fit_params = list(
        MODEL_FIT_PARAMS.get(
            model_name,
            DEFAULT_FIT_PARAMS,
        )
    )

    central_values = {
        p: float(info["central"])
        for p, info in ALL_PARAM_INFO.items()
    }

    central_values.update(
        MODEL_PARAM_CENTRAL.get(
            model_name,
            {}
        )
    )

    parameter_bounds = {
        p: tuple(bounds)
        for p, bounds in CUSTOM_PARAM_BOUNDS.items()
    }

    parameter_bounds.update(
        MODEL_PARAM_BOUNDS.get(
            model_name,
            {}
        )
    )

    missing_central = [
        p
        for p in fit_params
        if p not in central_values
    ]

    missing_bounds = [
        p
        for p in fit_params
        if p not in parameter_bounds
    ]

    if missing_central or missing_bounds:
        messages = []

        if missing_central:
            messages.append(
                "missing central values: "
                + ", ".join(missing_central)
            )

        if missing_bounds:
            messages.append(
                "missing bounds: "
                + ", ".join(missing_bounds)
            )

        raise ValueError(
            f"{model_name}: "
            + "; ".join(messages)
        )

    for p in fit_params:
        lower, upper = parameter_bounds[p]

        if not (
            lower > 0
            and upper > lower
        ):
            raise ValueError(
                f"{model_name}: log10 fitting requires "
                f"0 < lower < upper for {p}; "
                f"received ({lower}, {upper})."
            )

    return (
        fit_params,
        central_values,
        parameter_bounds,
    )


def get_model_species_setup(model_name):
    # Map the shared experimental species interface onto one model's SBML IDs.
    species_aliases = dict(
        MODEL_SPECIES_ALIASES.get(
            model_name,
            {},
        )
    )

    observable_id = MODEL_OBSERVABLE_IDS.get(
        model_name,
        OBSERVABLE_ID,
    )

    return species_aliases, observable_id


# ============================================================
# 11) SBML VALIDATION HELPERS
# ============================================================
# Before an expensive optimization begins, the script checks that the SBML
# contains every species needed for the experimental conditions, the chosen
# observable, and every configured fitted or fixed parameter.

def strip_namespace(tag):
    return tag.split("}", 1)[-1]


def inspect_sbml(model_file):
    # ElementTree gives a lightweight inventory before SBMLtoODEpy generates the executable Python ODE model.
    tree = ET.parse(model_file)

    species_ids = []
    global_parameter_ids = []
    reaction_ids = []

    for element in tree.getroot().iter():
        element_type = strip_namespace(element.tag)
        element_id = element.attrib.get("id")

        if not element_id:
            continue

        if element_type == "species":
            species_ids.append(element_id)

        elif element_type == "parameter":
            global_parameter_ids.append(element_id)

        elif element_type == "reaction":
            reaction_ids.append(element_id)

    return {
        "species_ids": species_ids,
        "global_parameter_ids": global_parameter_ids,
        "reaction_ids": reaction_ids,
    }


def validate_model(model_file):
    # Purpose of function:
    # Validate that the SBML file contains all species and parameters needed for the fitting process.

    model_name = model_file.stem

    (
        model_fit_params,
        _,
        _,
    ) = get_model_fit_setup(
        model_name
    )

    fixed_parameters = MODEL_FIXED_PARAMETERS.get(
        model_name,
        {},
    )

    inventory = inspect_sbml(model_file)

    species_aliases, observable_id = get_model_species_setup(
        model_name
    )

    required_species = {observable_id}

    # Validate species required for both fitting and held-out prediction.
    for d in datasets + HELD_OUT_DATASETS:
        required_species.update(
            species_aliases.get(species_id, species_id)
            for species_id in d["initial_values"]
        )

    missing_species = sorted(
        required_species.difference(inventory["species_ids"])
    )

    required_parameters = set(model_fit_params)
    required_parameters.update(fixed_parameters)

    missing_params = sorted(
        parameter_name
        for parameter_name in required_parameters
        if parameter_name not in inventory["global_parameter_ids"]
    )

    if missing_species or missing_params:
        messages = []

        if missing_species:
            messages.append(
                "missing species: " + ", ".join(missing_species)
            )

        if missing_params:
            messages.append(
                "missing configured global parameters: "
                + ", ".join(missing_params)
            )

        raise ValueError(
            f"{model_file.name}: " + "; ".join(messages)
        )

    try:
        model = SBMLtoODEpySolveIVPModel(
            sbml_file=model_file,
            generated_models_dir=GENERATED_MODELS_DIR,
        )
    except Exception as exc:
        raise RuntimeError(
            f"SBMLtoODEpy could not generate/import {model_file.name}: {exc}"
        ) from exc

    # Attach the model-specific interface so every downstream simulation uses
    # the correct SBML species and observable without changing the datasets.
    model.configured_species_aliases = species_aliases
    model.configured_observable_id = observable_id

    return model, inventory


# ============================================================
# 12) SIMULATION HELPERS
# ============================================================
# The SBML file is converted once into a Python class. 
# For every experimental condition, the adapter creates a fresh instance, 
# applies that condition's initial species concentrations and the current trial parameters, 
# and asks solve_ivp/LSODA for the trajectory at the exact experimental times.


def predict(
    rr,
    dataset,
    params,
    fixed_parameters,
):
    # Purpose of function:
    # Simulate one condition and return the configured observable.

    # The variable name ``rr`` is retained so the rest of the established fitting code changes as little as possible. 
    # It now holds an SBMLtoODEpySolveIVPModel adapter rather than a SBMLtoODEpy/solve_ivp object.

    t_eval_min = np.asarray(dataset["t_exp"], dtype=float)
    t_eval_model = t_eval_min * MODEL_TIME_UNITS_PER_MINUTE

    if t_eval_model[-1] <= t_eval_model[0]:
        return np.full_like(t_eval_min, 1e9, dtype=float)

    try:
        # Translate shared condition IDs into this model's local SBML IDs.
        model_initial_values = {
            rr.configured_species_aliases.get(species_id, species_id): value
            for species_id, value in dataset["initial_values"].items()
        }

        # Ignore shared defaults that do not belong to this model. Fitted and
        # fixed parameters are still validated before optimization begins.
        model_parameters = {
            parameter_id: value
            for parameter_id, value in params.items()
            if parameter_id in rr.parameter_ids
        }

        prediction = rr.simulate(
            initial_values=model_initial_values,
            parameters=model_parameters,
            fixed_parameters=fixed_parameters,
            t_eval=t_eval_model,
            observable_id=rr.configured_observable_id,
            method="LSODA",
            rtol=RTOL,
            atol=ATOL,
        )
    except Exception as exc:
        # A finite penalty lets least_squares reject a bad parameter region
        # without terminating the entire Sobol multistart fit.
        
        # The first failure is printed so a systematic model or programming error
        # is not silently mistaken for an unsuccessful parameter combination.
        # Later failures remain quiet to avoid printing thousands of messages.
        if not hasattr(predict, "first_error_reported"):
            print(
                "\nFirst simulation failure: "
                f"{type(exc).__name__}: {exc}"
            )
            predict.first_error_reported = True

        return np.full_like(t_eval_min, 1e9, dtype=float)

    prediction = np.asarray(prediction, dtype=float)

    if prediction.shape != t_eval_min.shape or np.any(~np.isfinite(prediction)):
        return np.full_like(t_eval_min, 1e9, dtype=float)

    return prediction


# ============================================================
# 13) RESIDUAL AND RMSE FUNCTIONS
# ============================================================
# The optimizer sees normalized residuals rather than raw concentration errors.
# For one curve:
#
#     normalized residual = (model - experiment)
#                           --------------------
#                           max(experiment) * sqrt(N)
#
# Dividing by max(experiment) prevents large-amplitude curves from dominating.
# Dividing by sqrt(N) gives each trajectory approximately equal total weight
# even when curves contain different numbers of fitted time points.

def compute_normalized_residual(
    rr,
    dataset,
    params,
    fixed_parameters,
):
    x_model = predict(
        rr=rr,
        dataset=dataset,
        params=params,
        fixed_parameters=fixed_parameters,
    )

    res = x_model - dataset["x_exp"]

    scale = (
        np.max(dataset["x_exp"])
        if np.max(dataset["x_exp"]) > 0
        else 1.0
    )

    points = len(dataset["x_exp"])

    return res / scale / np.sqrt(points)


def compute_rmse(
    rr,
    dataset,
    params,
    fixed_parameters,
):
    x_model = predict(
        rr=rr,
        dataset=dataset,
        params=params,
        fixed_parameters=fixed_parameters,
    )

    return float(
        np.sqrt(
            np.mean(
                (
                    x_model
                    - dataset["x_exp"]
                ) ** 2
            )
        )
    )


def compute_normalized_rmse(
    rr,
    dataset,
    params,
    fixed_parameters,
):
    rmse = compute_rmse(
        rr=rr,
        dataset=dataset,
        params=params,
        fixed_parameters=fixed_parameters,
    )

    scale = (
        np.max(dataset["x_exp"])
        if np.max(dataset["x_exp"]) > 0
        else 1.0
    )

    return float(rmse / scale)


def compute_global_rmse(
    rr,
    dataset_list,
    params,
    fixed_parameters,
):
    all_residuals = np.concatenate([
        compute_normalized_residual(
            rr=rr,
            dataset=d,
            params=params,
            fixed_parameters=fixed_parameters,
        )
        for d in dataset_list
    ])

    return float(
        np.sqrt(
            np.mean(all_residuals ** 2)
        )
    )


# ============================================================
# 14) FIT OR LOAD ONE MODEL
# ============================================================
# Parameters are optimized in log10 space. This is appropriate for positive
# kinetic parameters that span orders of magnitude and guarantees that trial
# values remain positive. The saved best_x therefore contains log10 values;
# physical parameter values are recovered with 10 ** best_x.

def fit_or_load_model(
    model_file,
    rr,
):
    model_name = model_file.stem

    (
        model_fit_params,
        model_central_values,
        model_param_bounds,
    ) = get_model_fit_setup(
        model_name
    )

    # Convert physical bounds into log10 space once. 
    # Every Sobol point and every optimizer step below uses this transformed coordinate system.
    lower_bounds = np.asarray([
        np.log10(
            model_param_bounds[p][0]
        )
        for p in model_fit_params
    ])

    upper_bounds = np.asarray([
        np.log10(
            model_param_bounds[p][1]
        )
        for p in model_fit_params
    ])

    model_results_dir = (
        RESULTS_DIR
        / model_name
    )

    model_figures_dir = (
        FIGURES_DIR
        / model_name
    )

    model_results_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    model_figures_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    fit_results_file = (
        model_results_dir
        / f"{model_name}_results.npz"
    )

    fixed_parameters = MODEL_FIXED_PARAMETERS.get(
        model_name,
        {},
    )

    def residuals(log_params):
        # Begin with a complete parameter dictionary. 
        # Only parameters listed for fitting are overwritten by the current optimizer trial.
        params = dict(
            model_central_values
        )

        for i, p in enumerate(
            model_fit_params
        ):
            params[p] = 10 ** log_params[i]

        return np.concatenate([
            compute_normalized_residual(
                rr=rr,
                dataset=d,
                params=params,
                fixed_parameters=fixed_parameters,
            )
            for d in datasets
        ])

    if MODE == "fit":
        # Sobol sequences are generated most evenly in powers of two.
        # The requested number of starts is therefore rounded up to 2**m.
        m = int(
            np.ceil(
                np.log2(
                    N_STARTS_REQUESTED
                )
            )
        )

        n_starts = 2 ** m

        print(
            f"\nRequested {N_STARTS_REQUESTED} starts for {model_name}."
            f"\nUsing {n_starts} Sobol starts.\n"
        )

        # Each Sobol point is one widely distributed initial guess in the
        # multidimensional bounded parameter space. least_squares then refines
        # that point locally. The lowest-cost refined solution is retained.
        sampler = qmc.Sobol(
            d=len(model_fit_params),
            scramble=True,
            seed=RANDOM_SEED,
        )

        U = sampler.random_base2(m=m)

        x0_all = (
            lower_bounds
            + U
            * (
                upper_bounds
                - lower_bounds
            )
        )

        best = None

        for i, x0 in enumerate(x0_all):
            result = least_squares(
                residuals,
                x0,
                bounds=(
                    lower_bounds,
                    upper_bounds,
                ),
                method=str(
                    config.get(
                        "optimizer_method",
                        "trf",
                    )
                ),
                max_nfev=MAX_NFEV,
            )

            print(
                f"{model_name} | "
                f"Run {i + 1}/{n_starts} | "
                f"cost = {result.cost:.4e} | "
                f"success = {result.success}"
            )

            if (
                best is None
                or result.cost < best.cost
            ):
                best = result

        if best is None:
            raise RuntimeError(
                f"No optimization result was produced for {model_name}."
            )

        best_x = np.asarray(
            best.x,
            dtype=float,
        )

        best_cost = float(best.cost)

        optimizer_success = bool(best.success)
        optimizer_message = str(best.message)

        np.savez(
            fit_results_file,
            best_x=best_x,
            best_cost=best_cost,
            fit_params=np.asarray(
                model_fit_params,
                dtype=object,
            ),
            model_file=str(model_file),
            optimizer_success=optimizer_success,
            optimizer_message=optimizer_message,
        )

        print(
            f"\nSaved fit results -> {fit_results_file}"
        )

    elif MODE == "load":
        # Loading is deliberately strict: 
        # the saved fitted-parameter order must match the current configuration because best_x is position-dependent.
        if not fit_results_file.exists():
            raise FileNotFoundError(
                f"Saved fit does not exist for {model_name}:\n"
                f"{fit_results_file}"
            )

        saved = np.load(
            fit_results_file,
            allow_pickle=True,
        )

        best_x = np.asarray(
            saved["best_x"],
            dtype=float,
        )

        best_cost = float(
            saved["best_cost"]
        )

        saved_fit_params = [
            str(parameter_name)
            for parameter_name in saved["fit_params"]
        ]

        if saved_fit_params != model_fit_params:
            raise ValueError(
                f"Saved fitted-parameter list for {model_name} "
                f"does not match the current config.\n"
                f"Saved: {saved_fit_params}\n"
                f"Config: {model_fit_params}"
            )

        optimizer_success = (
            bool(saved["optimizer_success"])
            if "optimizer_success" in saved.files
            else True
        )

        optimizer_message = (
            str(saved["optimizer_message"])
            if "optimizer_message" in saved.files
            else "Loaded result"
        )

        print(
            f"\nLoaded fit results -> {fit_results_file}"
        )

    else:
        raise ValueError(
            "MODE must be either 'fit' or 'load'."
        )

    best_fit_params = dict(
        model_central_values
    )

    for i, p in enumerate(
        model_fit_params
    ):
        best_fit_params[p] = 10 ** best_x[i]

    return {
        "model_name": model_name,
        "model_file": model_file,
        "rr": rr,
        "species_aliases": dict(rr.configured_species_aliases),
        "observable_id": rr.configured_observable_id,
        "fixed_parameters": fixed_parameters,
        "fit_params": model_fit_params,
        "param_bounds": model_param_bounds,
        "lower_bounds": lower_bounds,
        "upper_bounds": upper_bounds,
        "best_x": best_x,
        "best_cost": best_cost,
        "best_fit_params": best_fit_params,
        "optimizer_success": optimizer_success,
        "optimizer_message": optimizer_message,
        "fit_results_file": fit_results_file,
        "results_dir": model_results_dir,
        "figures_dir": model_figures_dir,
    }


# ============================================================
# 15) FIT EVERY MODEL
# ============================================================
# This is the generalized model-comparison loop. 
# The data and fitting rules are held constant while only the SBML mechanism changes.

model_results = []

for model_index, model_file in enumerate(
    MODEL_FILES,
    start=1,
):
    print("\n" + "=" * 72)
    print(
        f"MODEL {model_index}/{len(MODEL_FILES)}: "
        f"{model_file.name}"
    )
    print("=" * 72)

    rr, inventory = validate_model(
        model_file
    )

    print(
        f"Imported successfully | "
        f"species={len(inventory['species_ids'])} | "
        f"parameters={len(inventory['global_parameter_ids'])} | "
        f"reactions={len(inventory['reaction_ids'])}"
    )

    result = fit_or_load_model(
        model_file=model_file,
        rr=rr,
    )

    model_results.append(result)


# ============================================================
# 16) CALCULATE MODEL-SPECIFIC RMSE AND SAVE OUTPUTS
# ============================================================
# Optimization returns the best parameter vector and cost. 
# This block returns: 
# fitted trajectories
# per-curve RMSE
# global normalized RMSE
# parameter tables
# boundary checks
# figures
# human-readable summary for each candidate model.

# Note:
# This global normalised RMSE score is intended only for comparing models fitted to exactly
# the same experimental conditions, time points, fitting windows, data scaling,
# and residual definition. It should not be compared with results from runs
# that use different datasets or different numbers of fitted points.

comparison_rows = []

for result in model_results:
    model_name = result["model_name"]
    rr = result["rr"]
    observable_id = result["observable_id"]
    params = result["best_fit_params"]
    fixed_parameters = result["fixed_parameters"]

    per_curve_rows = []

    # Store held-out error statistics for this model.
    held_out_rows = []

    # Store held-out trajectories for later plotting and export.
    held_out_predictions = {}

    # Store held-out predictions in a long-form output table.
    held_out_prediction_rows = []

    print("\n" + "-" * 72)
    print(f"FIT SUMMARY: {model_name}")
    print("-" * 72)

    print("\nBest-fit parameters:")
    for p, value in params.items():
        print(f"{p:<15} = {value:.8e}")

    if fixed_parameters:
        print("\nModel-specific fixed parameters:")
        for p, value in fixed_parameters.items():
            print(f"{p:<15} = {value:.8e}")

    fitted_predictions = {}

    for d in datasets:
        prediction = predict(
            rr=rr,
            dataset=d,
            params=params,
            fixed_parameters=fixed_parameters,
        )

        fitted_predictions[
            d["name"]
        ] = prediction

        raw_rmse = compute_rmse(
            rr=rr,
            dataset=d,
            params=params,
            fixed_parameters=fixed_parameters,
        )

        normalized_rmse = compute_normalized_rmse(
            rr=rr,
            dataset=d,
            params=params,
            fixed_parameters=fixed_parameters,
        )

        per_curve_rows.append({
            "condition": d["name"],
            "raw_rmse": raw_rmse,
            "normalized_rmse": normalized_rmse,
            "n_fit_points": len(d["x_exp"]),
        })

        print(
            f"{d['name']:<24} | "
            f"RMSE={raw_rmse:.6e} | "
            f"normalized RMSE={normalized_rmse:.6e}"
        )

    # Evaluate conditions that were not used during fitting.
    if HELD_OUT_DATASETS:
        print("\nHELD-OUT PREDICTIONS")
        print("-" * 72)

    # Predict every held-out condition using the fitted parameters.
    for d in HELD_OUT_DATASETS:
        prediction = predict(
            rr=rr,
            dataset=d,
            params=params,
            fixed_parameters=fixed_parameters,
        )

        # Save the held-out trajectory for later export.
        held_out_predictions[
            d["name"]
        ] = prediction

        # Calculate raw held-out prediction error.
        raw_rmse = float(
            np.sqrt(
                np.mean(
                    (
                        prediction
                        - d["x_exp"]
                    ) ** 2
                )
            )
        )

        # Get the experimental scale used for normalization.
        scale = (
            np.max(d["x_exp"])
            if np.max(d["x_exp"]) > 0
            else 1.0
        )

        # Calculate normalized held-out prediction error.
        normalized_rmse = float(
            raw_rmse
            / scale
        )

        # Save the error statistics for this held-out condition.
        held_out_rows.append({
            "condition": d["name"],
            "raw_rmse": raw_rmse,
            "normalized_rmse": normalized_rmse,
            "n_prediction_points": len(
                d["x_exp"]
            ),
        })

        # Save every held-out trajectory point.
        for (
            time_value,
            experimental_value,
            predicted_value,
        ) in zip(
            d["t_exp"],
            d["x_exp"],
            prediction,
        ):
            held_out_prediction_rows.append({
                "condition": d["name"],
                "time_min": float(
                    time_value
                ),
                "experimental_value": float(
                    experimental_value
                ),
                "predicted_value": float(
                    predicted_value
                ),
            })

        # Report the held-out error.
        print(
            f"{d['name']:<24} | "
            f"RMSE={raw_rmse:.6e} | "
            f"normalized RMSE="
            f"{normalized_rmse:.6e}"
        )

    # Combine held-out errors with equal condition weighting.
    if held_out_rows:
        held_out_normalized_rmse = float(
            np.sqrt(
                np.mean([
                    row[
                        "normalized_rmse"
                    ] ** 2
                    for row in held_out_rows
                ])
            )
        )

    # Use a missing value when no conditions are held out.
    else:
        held_out_normalized_rmse = np.nan

    # Report the combined validation score.
    if held_out_rows:
        print(
            "\nCombined held-out normalized RMSE = "
            f"{held_out_normalized_rmse:.8e}"
        )

    global_normalized_rmse = compute_global_rmse(
        rr=rr,
        dataset_list=datasets,
        params=params,
        fixed_parameters=fixed_parameters,
    )

    print(
        f"\nGlobal normalized RMSE = "
        f"{global_normalized_rmse:.8e}"
    )

    # --------------------------------------------------------
    # Plot model fits
    # --------------------------------------------------------

    ncols = min(
        3,
        len(datasets),
    )

    nrows = int(
        np.ceil(
            len(datasets)
            / ncols
        )
    )

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(
            5.0 * ncols,
            4.0 * nrows,
        ),
        squeeze=False,
        sharey=True,
    )

    for ax, d in zip(
        axes.flat,
        datasets,
    ):
        ax.plot(
            d["t_exp"],
            d["x_exp"],
            "o",
            markersize=2.5,
            label="experiment",
        )

        ax.plot(
            d["t_exp"],
            fitted_predictions[d["name"]],
            linewidth=2.0,
            label="fit",
        )

        row = next(
            row
            for row in per_curve_rows
            if row["condition"] == d["name"]
        )

        ax.set_title(
            f"{d['name']}\n"
            f"nRMSE={row['normalized_rmse']:.3e}"
        )

        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Reacted reporter (nM)")
        ax.grid(True, alpha=0.3)

    for ax in axes.flat[len(datasets):]:
        ax.remove()

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper right",
    )

    fig.suptitle(
        f"{model_name}\n"
        f"Global normalized RMSE = "
        f"{global_normalized_rmse:.3e}"
    )

    fig.tight_layout()

    fit_plot_file = (
        result["figures_dir"]
        / f"{model_name}_fit.png"
    )

    fig.savefig(
        fit_plot_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Plot held-out predictions
    # --------------------------------------------------------

    # Create a held-out figure only when prediction datasets exist.
    if HELD_OUT_DATASETS:
        ncols_held_out = min(
            3,
            len(HELD_OUT_DATASETS),
        )

        nrows_held_out = int(
            np.ceil(
                len(HELD_OUT_DATASETS)
                / ncols_held_out
            )
        )

        # Create one panel for each held-out condition.
        fig_held_out, axes_held_out = plt.subplots(
            nrows_held_out,
            ncols_held_out,
            figsize=(
                5.0 * ncols_held_out,
                4.0 * nrows_held_out,
            ),
            squeeze=False,
            sharey=True,
        )

        # Plot experimental data and prediction for each condition.
        for ax, d in zip(
            axes_held_out.flat,
            HELD_OUT_DATASETS,
        ):
            ax.plot(
                d["t_exp"],
                d["x_exp"],
                "o",
                markersize=2.5,
                label="experiment",
            )

            ax.plot(
                d["t_exp"],
                held_out_predictions[
                    d["name"]
                ],
                linewidth=2.0,
                label="prediction",
            )

            # Find the error statistics for this condition.
            row = next(
                row
                for row in held_out_rows
                if row["condition"]
                == d["name"]
            )

            ax.set_title(
                f"{d['name']}\n"
                f"Held-out nRMSE="
                f"{row['normalized_rmse']:.3e}"
            )

            ax.set_xlabel("Time (min)")
            ax.set_ylabel(
                "Reacted reporter (nM)"
            )
            ax.grid(True, alpha=0.3)

        # Remove unused panels.
        for ax in axes_held_out.flat[
            len(HELD_OUT_DATASETS):
        ]:
            ax.remove()

        # Add the shared legend.
        handles, labels = (
            axes_held_out.flat[0]
            .get_legend_handles_labels()
        )

        fig_held_out.legend(
            handles,
            labels,
            loc="upper right",
        )

        # Add the model and validation score.
        fig_held_out.suptitle(
            f"{model_name}\n"
            f"Held-out normalized RMSE = "
            f"{held_out_normalized_rmse:.3e}"
        )

        fig_held_out.tight_layout()

        # Define the held-out figure path.
        held_out_plot_file = (
            result["figures_dir"]
            / f"{model_name}_held_out_prediction.png"
        )

        # Save the held-out prediction figure.
        fig_held_out.savefig(
            held_out_plot_file,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(fig_held_out)

        print(
            "Held-out prediction plot saved -> "
            f"{held_out_plot_file}"
        )

    # --------------------------------------------------------
    # Boundary proximity check
    # --------------------------------------------------------
    # Note: This part is a diagnostic only. It does not change the fit or the model ranking.

    # A parameter close to a search bound may indicate that the allowed range is too narrow, 
    # the parameter is weakly identifiable, or the model would prefer a value outside the assumed regime. 

    boundary_rows = []

    model_fit_params = result["fit_params"]
    model_param_bounds = result["param_bounds"]
    model_lower_bounds = result["lower_bounds"]
    model_upper_bounds = result["upper_bounds"]

    for i, p in enumerate(
        model_fit_params
    ):
        log_value = result["best_x"][i]
        log_lower = model_lower_bounds[i]
        log_upper = model_upper_bounds[i]

        parameter_range = (
            log_upper
            - log_lower
        )

        distance_to_lower = (
            log_value
            - log_lower
        ) / parameter_range

        distance_to_upper = (
            log_upper
            - log_value
        ) / parameter_range

        if distance_to_lower < 0.05:
            status = "near lower bound"
        elif distance_to_upper < 0.05:
            status = "near upper bound"
        else:
            status = "not within 5% of bound"

        boundary_rows.append({
            "parameter": p,
            "value": params[p],
            "lower": model_param_bounds[p][0],
            "upper": model_param_bounds[p][1],
            "distance_to_lower": distance_to_lower,
            "distance_to_upper": distance_to_upper,
            "status": status,
        })

    # --------------------------------------------------------
    # Save tables
    # --------------------------------------------------------

    parameter_rows = [
        {
            "parameter": p,
            "value": params[p],
            "fitted": p in model_fit_params,
        }
        for p in params
    ]

    for p, value in fixed_parameters.items():
        parameter_rows.append({
            "parameter": p,
            "value": value,
            "fitted": False,
        })

    # Save fitted parameters.
    pd.DataFrame(
        parameter_rows
    ).to_csv(
        result["results_dir"]
        / "fitted_parameters.csv",
        index=False,
    )

    # Save errors for conditions used during fitting.
    pd.DataFrame(
        per_curve_rows
    ).to_csv(
        result["results_dir"]
        / "per_curve_rmse.csv",
        index=False,
    )

    # Save held-out results when prediction conditions exist.
    if held_out_rows:
        pd.DataFrame(
            held_out_rows
        ).to_csv(
            result["results_dir"]
            / "held_out_rmse.csv",
            index=False,
        )

        # Save every point in the held-out trajectories.
        pd.DataFrame(
            held_out_prediction_rows
        ).to_csv(
            result["results_dir"]
            / "held_out_predictions.csv",
            index=False,
        )

    # Save the parameter-boundary diagnostic.
    pd.DataFrame(
        boundary_rows
    ).to_csv(
        result["results_dir"]
        / "boundary_check.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Save text summary
    # --------------------------------------------------------

    summary_file = (
        result["results_dir"]
        / f"{model_name}_summary.txt"
    )

    with open(
        summary_file,
        "w",
        encoding="utf-8",
    ) as f:
        f.write("RUN SUMMARY\n")
        f.write("=" * 60 + "\n\n")

        f.write(
            f"Run name       : {RUN_NAME}\n"
        )
        f.write(
            f"Model          : {model_name}\n"
        )
        f.write(
            f"SBML file      : {result['model_file']}\n"
        )
        f.write(
            f"Mode           : {MODE}\n"
        )
        f.write(
            f"Fit conditions : {FIT_CONDITIONS}\n"
        )
        f.write(
            f"Held-out conditions: "
            f"{HELD_OUT_CONDITIONS}\n"
        )
        f.write(
            f"Fit parameters : {model_fit_params}\n"
        )
        f.write(
            f"N starts       : {N_STARTS_REQUESTED}\n"
        )
        f.write(
            f"Random seed    : {RANDOM_SEED}\n"
        )
        f.write(
            f"Observable     : {observable_id}\n"
        )
        f.write(
            f"Time conversion: "
            f"{MODEL_TIME_UNITS_PER_MINUTE} "
            f"model units per minute\n\n"
        )

        f.write("BEST-FIT PARAMETERS\n")
        f.write("-" * 60 + "\n")

        for p, value in params.items():
            f.write(
                f"{p:<20} = {value:.8e}\n"
            )

        if fixed_parameters:
            f.write(
                "\nMODEL-SPECIFIC FIXED PARAMETERS\n"
            )
            f.write("-" * 60 + "\n")

            for p, value in fixed_parameters.items():
                f.write(
                    f"{p:<20} = {value:.8e}\n"
                )

        f.write("\nFIT STATISTICS\n")
        f.write("-" * 60 + "\n")
        f.write(
            f"Best cost (0.5 * SSE)  = "
            f"{result['best_cost']:.8e}\n"
        )
        f.write(
            f"Global normalized RMSE = "
            f"{global_normalized_rmse:.8e}\n"
        )
        f.write(
            f"Optimizer success      = "
            f"{result['optimizer_success']}\n"
        )
        f.write(
            f"Optimizer message      = "
            f"{result['optimizer_message']}\n"
        )

        f.write("\nRMSE PER DATASET\n")
        f.write("-" * 60 + "\n")

        for row in per_curve_rows:
            f.write(
                f"{row['condition']:<24} | "
                f"RMSE={row['raw_rmse']:.8e} | "
                f"nRMSE={row['normalized_rmse']:.8e}\n"
            )

        f.write(
            "\nHELD-OUT PREDICTION RMSE\n"
        )
        f.write("-" * 60 + "\n")

        # Record every held-out prediction error.
        for row in held_out_rows:
            f.write(
                f"{row['condition']:<24} | "
                f"RMSE={row['raw_rmse']:.8e} | "
                f"nRMSE="
                f"{row['normalized_rmse']:.8e}\n"
            )

        # Record the combined held-out score.
        if held_out_rows:
            f.write(
                "Combined held-out nRMSE = "
                f"{held_out_normalized_rmse:.8e}\n"
            )
        else:
            f.write("None\n")

        f.write("\nBOUNDARY PROXIMITY\n")
        f.write("-" * 60 + "\n")

        for row in boundary_rows:
            f.write(
                f"{row['parameter']:<15} | "
                f"value={row['value']:.8e} | "
                f"{row['status']}\n"
            )

    result[
        "global_normalized_rmse"
    ] = global_normalized_rmse

    result[
        "per_curve_rows"
    ] = per_curve_rows

    # Store held-out results with this model.
    result[
        "held_out_normalized_rmse"
    ] = held_out_normalized_rmse

    result[
        "held_out_rows"
    ] = held_out_rows

    comparison_rows.append({
        "model": model_name,
        "model_file": result[
            "model_file"
        ].name,
        "number_of_fitted_parameters": len(
            model_fit_params
        ),
        "best_cost": result[
            "best_cost"
        ],
        "global_normalized_rmse": (
            global_normalized_rmse
        ),
        "held_out_normalized_rmse": (
            held_out_normalized_rmse
        ),
        "optimizer_success": result[
            "optimizer_success"
        ],
    })


# ============================================================
# 17) COMPARE AND SELECT MODELS
# ============================================================
# All candidates have now been evaluated using the same data and residual definition. 

comparison = pd.DataFrame(
    comparison_rows
)

selection_metric = str(
    config.get(
        "selection_metric",
        "global_normalized_rmse",
    )
)

if selection_metric not in comparison.columns:
    raise KeyError(
        f"Selection metric {selection_metric!r} "
        f"is not in comparison table."
    )

# Require prediction-only conditions for validation selection.
if (
    selection_metric
    == "held_out_normalized_rmse"
    and not HELD_OUT_DATASETS
):
    raise ValueError(
        "selection_metric is "
        "'held_out_normalized_rmse', "
        "but fit_conditions contains "
        "every condition."
    )

comparison = comparison.sort_values(
    by=[
        selection_metric,
        "model",
    ],
    ascending=[
        True,
        True,
    ],
).reset_index(drop=True)

comparison["rank"] = (
    np.arange(len(comparison))
    + 1
)

comparison["selected"] = (
    comparison["rank"]
    == 1
)

comparison_file = (
    RESULTS_DIR
    / "model_comparison.csv"
)

comparison.to_csv(
    comparison_file,
    index=False,
)

selected_row = comparison.iloc[0]
selected_model_name = str(
    selected_row["model"]
)

selected_result = next(
    result
    for result in model_results
    if result["model_name"]
    == selected_model_name
)

selected_model_record = {
    "run_name": RUN_NAME,
    "selection_metric": selection_metric,
    "model_name": selected_model_name,
    "model_file": str(
        selected_result["model_file"]
    ),
    "fit_results_file": str(
        selected_result["fit_results_file"]
    ),
    "fitted_parameters_file": str(
        selected_result["results_dir"]
        / "fitted_parameters.csv"
    ),
    "global_normalized_rmse": float(
        selected_row[
            "global_normalized_rmse"
        ]
    ),
    "held_out_normalized_rmse": (
        float(
            selected_row[
                "held_out_normalized_rmse"
            ]
        )
        if pd.notna(
            selected_row[
                "held_out_normalized_rmse"
            ]
        )
        else None
    ),
    "observable_id": selected_result["observable_id"],
    "species_aliases": selected_result["species_aliases"],
    "model_time_units_per_minute": (
        MODEL_TIME_UNITS_PER_MINUTE
    ),
}

selected_model_json = (
    RESULTS_DIR
    / "selected_model.json"
)

with open(
    selected_model_json,
    "w",
    encoding="utf-8",
) as f:
    json.dump(
        selected_model_record,
        f,
        indent=2,
    )

selected_model_copy = (
    RESULTS_DIR
    / "selected_model.xml"
)

shutil.copy2(
    selected_result["model_file"],
    selected_model_copy,
)

print("\n" + "=" * 72)
print("MODEL COMPARISON")
print("=" * 72)
print(
    comparison.to_string(
        index=False
    )
)

print(
    f"\nSelected model: "
    f"{selected_model_name}"
)

print(
    f"Comparison saved -> "
    f"{comparison_file}"
)

print(
    f"Selection record saved -> "
    f"{selected_model_json}"
)

print(
    f"Selected SBML copied -> "
    f"{selected_model_copy}"
)
