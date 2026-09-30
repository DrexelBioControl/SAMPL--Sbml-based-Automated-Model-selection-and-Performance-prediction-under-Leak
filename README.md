# `SAMPLE-CRN`: 
# `S`emi `A`utomated `M`odel selection and `P`erformance prediction for `L`eaky `E`xperimental `C`hemical `R`eaction `N`etworks
## _A Python workflow for model selection, model validation, and performance prediction of leaky biochemical circuits._

## Start the workflow

From the project folder, run:

```powershell
python .\Run_SAMPL.py
```

The launcher first asks which fitting JSON to use:

```text
Fitting settings filename [fit_all_models.json]:
```

Type a filename from `configs/`, such as `my_fitting_settings.json`.
The `.json` extension is optional. A full path or a project-relative path
also works. Press Enter to use `fit_all_models.json`.

It then asks whether to **fit** or **load**:

| Mode | What the launcher does |
|---|---|
| `fit` | Detects the Excel file and candidate models, then fits and compares the models. |
| `load` | Asks for a previous results folder and reuses its selected model and fitted parameters. |

It then asks whether to perform LOO and performance prediction. For prediction,
it lists the metric files in `metrics/`, followed by the matching JSON files in
`configs/`. Type a filename or its number, or press Enter to accept the displayed
default. The defaults are fit, no LOO, yes prediction, and amplification.

In **fit mode**, model identification fits every condition in `all_conditions`:
currently all six conditions. The main fitting JSON only needs this list; the
launcher creates `fit_conditions` in the saved settings copy. If the ranking
setting was `held_out_normalized_rmse`, the launcher uses
`global_normalized_rmse` because identification has no held-out condition.
Other ranking choices are retained.

In **load mode**, enter the results folder containing the previous selection:

```text
results
```

An absolute path also works. If the selection has a prefix, give the JSON
filename itself, for example `results/experiment_A_selected_model.json`.
Older workflow run folders and their `identification/results` folders are
still accepted. The saved SBML and fitted-parameter NPZ must remain available.

The launcher reads these existing files without copying them or comparing
the candidates again. It saves `results/loaded_model.json`, a small record
pointing to the files it reused. Keep those source files available. A fit
loaded from an older run may have used a subset of conditions.

**LOO always performs fresh fits in either mode.** It refits only the selected
model, using the current fitting JSON's Excel data, conditions and parameter
settings. With the supplied six conditions, it fits five and predicts the
omitted condition, repeated six times. Use the same experimental setup when
loading a previous model and requesting LOO.

Performance prediction uses the fit from model identification, or the loaded
fit, never a LOO fold's fit. Load mode without LOO does not need the original
Excel file or candidate-model folder. The metric JSON must use species and
parameter IDs compatible with the selected model.

The supplied prediction examples use `ROL`, `F_temp` and/or `basal_frac` from
`model1_activegate_fuel`. If another candidate wins, those examples need the
appropriate IDs for that model. Code 2 reports incompatible settings in
the prediction log; the identification and LOO results remain saved.

To supply the fitting JSON directly and skip its filename prompt, or check
the configured input paths without running the calculations:

```powershell
python .\Run_SAMPL.py --fit-config configs/my_fitting_settings.json
python .\Run_SAMPL.py --fit-config configs/fit_all_models.json --check
```

`--check` checks the chosen fitting inputs, then exits before the fit/load
and metric prompts. With no `--fit-config`, it first asks for the filename.
The fitting JSON supplies the shared output folders in either mode.

The launcher uses the existing fitting and metric JSON files; it needs no
separate workflow configuration. All calculations use the same three output
folders. Filenames identify the model, metric or omitted condition:

```text
figures/
    fitting_windows.png
    <model>_fit.png
    loo_01_<model>_held_out_prediction.png
    <model>_amplification.png
    <model>_sensitivity.png

results/
    selected_model.json, selected_model.xml
    <model>_fit.npz, <model>_parameters.csv
    model_comparison.csv
    loo_01_<model>_fit.npz, loo_01_<model>_held_out_rmse.csv
    <model>_loo_summary.csv, <model>_loo_summary.json
    <model>_amplification.npz, <model>_sensitivity.npz
    ...saved settings, diagnostics, logs and workflow_summary.json

generated_models/
    ...Python equations generated from the SBML files
```

LOO and prediction files appear when those steps are requested. `loo_01_`
means the first condition in `all_conditions` was omitted; its saved settings
also record the condition names. Each fold has separate files and does not
replace the identification fit or announce another selected model.

**Running the same calculation again replaces files with the same names.**
It does not delete other outputs. To keep another version, add a prefix to
the fitting JSON before launching:

```json
"output_prefix": "experiment_A_"
```

The prefix is added exactly as written, including the final underscore.
For example, the fit becomes `experiment_A_<model>_fit.npz` and its selection
becomes `experiment_A_selected_model.json`. The launcher uses the same prefix
for LOO and prediction. Omit it, or use `""`, for the usual filenames.
`run_name` is a descriptive label; it does not change the filenames.

The Excel file and source JSON files in `configs/` stay editable. The launcher
saves the fitting settings only when it performs a fresh identification.
Loading a fit leaves that saved fitting-settings file unchanged. Each LOO
fit and prediction saves its own settings and log in `results/`. Figures
are saved without opening plot windows.

Generated equations are shared between runs. The adapter checks the SBML
contents and regenerates the Python model when those contents change.

## Contents

- [Start the workflow](#start-the-workflow)
- [Workflow and file map](#workflow-and-file-map)
- [Important terms](#important-terms)
- [Data file structure](#data-file-structure)
- [Code 1: fit and compare models](#code-1-fit-and-compare-models)
- [Code 2: predict circuit performance](#code-2-predict-circuit-performance)
- [To perform a leave-one-out (LOO) model validation](#to-perform-a-leave-one-out-loo-model-validation)
- [Advanced settings and troubleshooting](#advanced-settings-and-troubleshooting)

## Workflow and file map

`Run_SAMPL.py` asks for the run choices and carries out these steps:

```text
Data + candidate SBML models + fitting settings
                        |
                        v
                  Fit mode: Code 1
              Fits and compares models
                        |
                        +<--- Load mode: previous selected model + fitted parameters
                        |
             Selected model + fitted parameters ---> Optional fresh LOO fits
                        |
                        v
Run JSON -----------> Code 2 <--------- Shared model defaults
                        |
                        +--> Metric file: build_requests(...)
                        |    Specifies properties and their simulation conditions
                        |
                        +--> properties/model_properties.py
                        |    Prepares each condition and calls its calculation:
                        |      time_course_simulation.py -> concentrations or rates
                        |      steady_state_simulation.py -> steady state or plateau
                        |      a user-defined property file -> another quantity
                        |
                        +--> Metric file: calculate_metric(...)
                        |    Applies the performance formula to those properties
                        |
                        v
              Figure + numerical results + settings
```

A **model property** is a quantity calculated from the model, such as a
concentration, a rate, or a verified limiting output. A **performance metric**
specifies which quantities it needs and combines them mathematically.

For example, amplification needs rates under four sets of conditions:
ON, OFF, reference ON, and reference OFF. Sensitivity needs the limiting
output over a parameter grid. These choices belong to their metric files.
Code 2 runs the requests it receives without requiring every metric to use
ON/OFF conditions or a reference simulation.

| File or group | What it contains or does |
|---|---|
| [`Run_SAMPL.py`](Run_SAMPL.py) | The file to launch the complete workflow: asks for fit/load, LOO and prediction choices, then runs the requested steps. |
| [`workflow_steps.py`](workflow_steps.py) | Runs Code 1/2 in the same Python environment and repeats the selected-model LOO fits. |
| `models/*.xml` | Model equations, species, parameters, compartments, and initial values. |
| `data/` | Experimental time courses. |
| [`configs/fit_all_models.json`](configs/fit_all_models.json) | Code 1's data selection, fitting conditions, model mappings, and fitting settings. |
| [`Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py`](Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py) | Fits and compares models, then saves the selected model and parameter estimates. |
| `results/selected_model.json` and `results/<model>_fit.npz` | Code 1's model selection and fitted parameters, read by Code 2. |
| `results/loaded_model.json` | Points to the existing model and fitted parameters when using the launcher in load mode. |
| [`Code2_MasterSBMLPrediction_GeneralInterfaces.py`](Code2_MasterSBMLPrediction_GeneralInterfaces.py) | Loads the chosen run, executes property requests, calls the metric, and plots and saves results. This is the Code 2 file to run. |
| [`configs/model_interfaces.json`](configs/model_interfaces.json) | Shared prediction defaults for each model: initial concentrations, optional fixed parameters, display units, and an optional input mapping. |
| [`configs/amplification_prediction.json`](configs/amplification_prediction.json) | Amplification run choices, including metric-specific ON/OFF and reference settings. |
| [`configs/sensitivity_prediction.json`](configs/sensitivity_prediction.json) | Sensitivity run choices, including species, parameter grid, initial conditions, and convergence settings. |
| [`code2_configuration.py`](code2_configuration.py) | Combines the selected shared defaults with this run's overrides. |
| [`properties/model_properties.py`](properties/model_properties.py) | Reads requests, prepares conditions, calls the calculations, and collects the results across each sweep. |
| [`properties/time_course_simulation.py`](properties/time_course_simulation.py) | Calculates concentrations and rates at the requested times for one prepared condition. |
| [`properties/steady_state_simulation.py`](properties/steady_state_simulation.py) | Integrates one prepared condition over successive time windows and checks for a steady state or output plateau. |
| [`sbmltoodepy_solveivp_adapter.py`](sbmltoodepy_solveivp_adapter.py) | Connects the SBML-generated equations to numerical integration. |
| [`metrics/amplification.py`](metrics/amplification.py), [`metrics/sensitivity.py`](metrics/sensitivity.py) | Define each metric's requests, formula, validity rules, and labels. |
| [`metrics/MetricTemplate.py`](metrics/MetricTemplate.py) | A simple example that requests concentrations and returns the last requested value. |
| [`properties/PropertyTemplate.py`](properties/PropertyTemplate.py) | An example of an imported property: the initial observable concentration at each control value. |
| `tests/` | Checks for metric formulas, requests, configuration, convergence, and workflow setup. |
| `generated_models/` | Python models generated from SBML. Change the source SBML to change the model. |
| `figures/` | Shared folder for fitting, LOO and performance figures. |
| `results/` | Shared folder for numerical results, selection records, saved settings and logs. |

### Reading the property files

```text
properties/
    __init__.py                 Marks this as a Python package
    model_properties.py         Coordinates requests and sweeps
    time_course_simulation.py   Calculates concentrations and rates
    steady_state_simulation.py  Calculates steady states and output plateaus
    PropertyTemplate.py        Shows how to add another calculation
```

Start with `collect_properties(...)` in `properties/model_properties.py`.
It sends each request to `collect_time_properties(...)`,
`collect_limit_property(...)`, or the user's imported property function.
The first two prepare one condition at a time and call the appropriate
simulation file. Read that simulation file next if you want to follow its
numerical calculation; the coordinator collects its answers into arrays.

Concentrations and rates stay together because they use the same trajectory.
Steady state and output plateau stay together because they use the same
integration loop with different convergence checks. All are ordinary
functions; these files are imported automatically by Code 2.

### Before and after

| Part | Original Code 2 | Current Code 2 |
|---|---|---|
| Main file | Amplification-specific name and calculation | General prediction entry point |
| Simulation conditions | ON/OFF and reference conditions built into Code 2 | Each metric builds its own named requests |
| Performance formula | Amplification formula inside Code 2 | Formula in the chosen metric file |
| Quantities available | Rates used for amplification | Concentrations, rates, verified limits, or imported properties |
| Property calculations | Simulation and sweep logic mixed together | Coordinator and two simulation files together in `properties/` |
| Model defaults | Repeated in run configurations | Shared model-interface JSON |
| Numerical settings | Integration settings in Code 2 | General settings read from the run JSON |

The amplification equation and numerical choices are unchanged. Code 1's
fitting calculations are unchanged; the launcher supplies its settings and
repeats it for automatic LOO. The shared model interfaces are used by Code 2;
Code 1 continues to use its fitting JSON.

### Where to make a change

| What you want to change | Where it belongs |
|---|---|
| Reactions or model equations | The SBML file; refit with Code 1 when needed |
| Shared prediction defaults for a model | `configs/model_interfaces.json` |
| Observable, main control grid, initial conditions, or parameter overrides | The run JSON |
| Sampling times, solver tolerances, or convergence settings | The run JSON |
| Metric-specific choices, such as ON/OFF values or a rate floor | The run JSON's `metric.settings` |
| Which simulations a metric needs | Its `build_requests(...)` function |
| How the requested quantities are combined | Its `calculate_metric(...)` function |
| A new type of model property | A property Python file, registered under `property_definitions` |

The Python files start with a short purpose and input/output outline.
Local comments explain the implementation as it proceeds. This README
contains the fuller map and configuration examples.

Set up the environment once:

```powershell
conda env create -f environment.yml
conda activate SBMLtoODEpyWorkflow
```

## Important terms

| Term | Meaning |
|---|---|
| Fitting condition | A complete experimental curve used to estimate parameters. |
| Held-out condition | A complete curve predicted without using it during fitting. |
| Shared species name | A model-independent experimental role used in Code 1's conditions. |
| Species alias | Translation from a shared name to an exact model-specific SBML ID. |
| Observable | The species whose concentration or rate is analysed. |
| Model property | A quantity calculated from the model, such as a rate or limiting output. |
| Performance metric | A mathematical expression that uses model properties to quantify performance. |
| Analysis control | The species initial concentration or parameter varied along the main Code 2 grid. |
| Property request | A named calculation: which property to calculate, under which conditions and control values. |
| Request alias | A name chosen by the metric for a request, such as `on_rates` or `response`. |
| Secondary parameter sweep | Repeats the main-grid analysis at another parameter value, usually in a separate panel. |
| ON/OFF input | Amplification's two input conditions, defined in its metric settings. |
| Amplification reference | The control value used for amplification's baseline rates; separate from the OFF input. |

## Data file structure

Code 1 assumes that the experimental data are stored in one Excel worksheet with the following structure:

| Excel row | Column A | Column B | Column C | Additional columns |
|---|---|---|---|---|
| 1 | Blank | Condition 1 name | Condition 2 name | ... |
| 2 | Time label and units | Measurement label | Measurement label | ... |
| 3 | Numeric time values | Condition 1 measurements | Condition 2 measurements | ... |
| ... | ... | ... | ... | ... |

The first column contains the time points shared by all conditions. Each remaining column contains one complete experimental time course.

For example:

|  | Condition 1 | Condition 2 |
|---|---:|---:|
| Time (min) | Signal | Signal |
| 0 | 0.000 | 0.000 |
| 5 | 0.125 | 0.080 |
| 10 | 0.260 | 0.170 |

Each `data_column` value in the Code 1 JSON must exactly match a condition name in the first Excel row. Measurements are multiplied by `data_signal_multiplier` after loading; use `1.0` when the spreadsheet values are already in the model's required units. Replicate columns are not automatically averaged.


## Code 1: fit and compare models

Edit the existing [`configs/fit_all_models.json`](configs/fit_all_models.json). It is the complete template. Work through the sections below in order.

### 1. Choose models and data

Set:

- `models_dir` and `model_pattern` for the candidate XML files;
- `excel_file` and `excel_sheet` for the measurements; and
- `run_name` to identify the analysis.

Each model's configuration name must equal its XML filename without `.xml`.

### 2. Define conditions and translate species names

```json
{
  "all_conditions": ["condition_high", "condition_low"],
  "conditions": {
    "condition_high": {
      "data_column": "High condition",
      "initial_values": {
        "varied_species": 50.0,
        "fixed_species": 25.0
      }
    },
    "condition_low": {
      "data_column": "Low condition",
      "initial_values": {
        "varied_species": 10.0,
        "fixed_species": 25.0
      }
    }
  },
  "model_species_aliases": {
    "model_1": {
      "varied_species": "model1_species_A",
      "fixed_species": "model1_species_B"
    },
    "model_2": {
      "varied_species": "model2_species_X",
      "fixed_species": "model2_species_Y"
    }
  },
  "observable_id": "shared_output",
  "model_observable_ids": {
    "model_2": "model2_output"
  }
}
```

How to read this example:

- Each condition is an entire data curve, not one data point.
- `all_conditions` lists the curves to analyse. Code 1 fits them all when `fit_conditions` is omitted.
- `data_column` identifies that curve in the spreadsheet.
- `initial_values` lists only SBML species that must be set at time zero.
- `varied_species` changes between conditions; `fixed_species` does not.
- The left side of an alias must match a name in `initial_values`.
- The right side must exactly match the model's SBML species ID.
- The observable is configured separately from `initial_values`.

For `condition_low`, Code 1 sets the mapped varied species to `10` and the mapped fixed species to `25` in both models.

Use biologically meaningful shared names when possible, such as `input_template`, `gate_template`, or `substrate_pool`. Other internal species should retain their values from the SBML file.

Every shared initial-value name must exist in every candidate model, either directly or through an alias. If a species exists in only one candidate, normally leave its value in that model's SBML file.

### 3. Define fitted parameters

Use shared settings when models use the same parameter IDs. Add model-specific entries only for exceptions.

```json
{
  "fit_params": ["shared_rate"],
  "all_param_central": {
    "shared_rate": 0.001
  },
  "custom_param_bounds": {
    "shared_rate": [0.000001, 0.1]
  },
  "model_fit_params": {
    "model_2": ["model2_rate"]
  },
  "model_param_central": {
    "model_2": {
      "model2_rate": 0.002
    }
  },
  "model_param_bounds": {
    "model_2": {
      "model2_rate": [0.000001, 0.1]
    }
  }
}
```

| Field | Meaning |
|---|---|
| `fit_params` | Parameters fitted in models that use the shared IDs. |
| `all_param_central` | Starting or default values for shared parameters. |
| `custom_param_bounds` | Lower and upper bounds for shared parameters. |
| `model_fit_params` | Replacement fit list for a named model. |
| `model_param_central` | Starting or default values for that model. |
| `model_param_bounds` | Bounds for that model. |
| `model_fixed_parameters` | Optional parameter values that are set but not fitted. |

All fitted bounds must satisfy `0 < lower < upper` because fitting is performed in log space.

### 4. Choose validation, ranking, and fitting effort

```json
{
  "selection_metric": "global_normalized_rmse",
  "n_starts": 8,
  "random_seed": 20260806,
  "max_nfev": 400
}
```

| Setting | Meaning |
|---|---|
| `selection_metric` | Use `global_normalized_rmse` when fitting all conditions. Use `held_out_normalized_rmse` only when some conditions are held out. |
| `n_starts` | Number of starting parameter sets tried per model. More starts may find a better solution but take longer. |
| `random_seed` | Makes the starting sets reproducible. Keep it fixed when comparing models. |
| `max_nfev` | Maximum model evaluations per start. Larger values allow more fitting work but take longer. |

A condition in `all_conditions` but outside an explicitly supplied
`fit_conditions` list is held out. If `fit_conditions` is omitted, there are
no held-out conditions.

Code 1 also includes the fitting diagnostics from ver1:

- **Sobol points:** a fresh fit saves `<model>_sobolsamples.png` in `figures/`.
  Each row shows the starting values, the best fit, and the parameter bounds.
  Load mode skips this plot because it does not generate new starting points.
- **Cost landscapes:** set `"save_cost_landscapes": true` in the fitting JSON
  to save `<model>_cost_landscapes.png`. This works in Code 1's fit and load modes.
  Each curve varies one parameter while holding the others at their best-fit values.
  It uses 50 points, with a log scan from 1000 times below to 1000 times above the
  best fit. `basal_frac` and `misfold_frac` instead use a linear scan from 0 to 1.
  These diagnostic scans can go beyond the fitting bounds, which are shaded.
  The vertical axis shows the change in normalized SSE from the lowest sampled
  value on that curve. This is twice the change in `least_squares` cost.
- **Boundary checks:** `basal_frac` and `misfold_frac` use linear distances to
  their fitting bounds; other parameters use log10 distances. The CSV and text
  summary state which scale was used. Within 5% of the range counts as near a bound.

The two fraction names are listed in `LINEAR_SCAN_RANGES` and
`LINEAR_BOUNDARY_PARAMS` in Code 1's section 16. Add another fraction name there
if needed. These settings affect diagnostics only: Sobol sampling and fitting
still use log10 values for every fitted parameter, as in ver1. The cost curves
are slices with other parameters fixed, not confidence intervals or a formal
identifiability test.

The optional output prefix also applies to these figures. Automatic LOO uses
the same diagnostic settings for each fold and adds its `loo_01_`, etc. prefix.
Cost scans add simulations, so `save_cost_landscapes` can stay `false` for a
quicker run. They use only that run's fitting conditions, not its held-out data.

### 5. Run Code 1 by itself

To fit all conditions without the launcher:

1. List the conditions in `all_conditions` in `configs/fit_all_models.json`.
2. Leave out `fit_conditions` and use `"selection_metric": "global_normalized_rmse"`.
3. Run this from the project folder:

```powershell
conda activate SBMLtoODEpyWorkflow
python .\Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py --config configs/fit_all_models.json --mode fit
```

Code 1 fits and compares the candidate models. It does not automatically repeat
LOO folds or run performance prediction. Results and figures go to the JSON's
`results_dir` and `figures_dir` (normally `results/` and `figures/`).

To fit only some conditions, add `fit_conditions` to the JSON. For example,
with the supplied six conditions, this holds out `fig2d_2p5nM_I1`:

```json
"fit_conditions": [
  "fig2d_50nM_I0",
  "fig2d_50nM_I1",
  "fig2d_25nM_I1",
  "fig2d_12p5nM_I1",
  "fig2d_6p25nM_I1"
]
```

Keep all six in `all_conditions`. Code 1 fits the five listed above and predicts
the remaining curve. Use `held_out_normalized_rmse` if ranking models by that
held-out prediction. Run the same command with `--mode fit` after changing the
conditions. To return to fitting all conditions, remove `fit_conditions` and
set the ranking back to `global_normalized_rmse`.

Other ways to run Code 1:

```powershell
python .\Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py
python .\Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py --config configs/fit_all_models.json --mode load --model models/model1_activegate_fuel.xml
```

With no options, Code 1 uses `MODE` and `CONFIG_FILE` near the start of the file.
`--model` restricts the candidate files and can be repeated. Relative paths
start at the project folder. Use `load` only for results from the same fitting
setup.

Code 1's `--mode load` reads the saved fits for the requested candidate models
and compares them again using the current settings. The launcher's **load**
choice simply reuses a previously selected model and its saved parameters.

The terminal shows the run choices, fitted and held-out conditions, progress
through the Sobol starts, model scores, and output folders. Detailed parameters,
per-condition errors, optimizer messages, and boundary checks stay in
`results/<model>_summary.txt` and the corresponding CSV files. Simulation
failures and unsuccessful optimizer results are still reported.

Code 1 still checks the required data columns, SBML species and parameter IDs,
positive fitting bounds, valid simulation outputs, and saved parameter order.
These checks protect which data and parameter values enter the calculation.

Check these outputs, whether running Code 1 directly or through the launcher:

- `results/model_comparison.csv`: model scores and ranks;
- `results/selected_model.json`: selected model and fitted-results path;
- `results/<model>_parameters.csv`: readable parameter values;
- `results/<model>_held_out_rmse.csv`: error when conditions are held out;
- `figures/<model>_fit.png`: fitted curves; and
- `figures/<model>_held_out_prediction.png`: curves when conditions are held out.

The machine-readable fitted parameters are stored in `results/<model>_fit.npz`
and loaded by Code 2. An optional `output_prefix` is added to these filenames,
including `selected_model.json`, its saved XML and `model_comparison.csv`.

Direct Code 1 load mode reads `results/<prefix><model>_fit.npz`, using the
results folder and optional prefix from the fitting JSON. It reads the optimizer
status saved in that file. It no longer searches the older per-model subfolders
or assumes success when the saved optimizer status is missing. Files saved by
the current Code 1 already contain this status. New tables and plots go to the
shared folders. The launcher's separate loading behaviour is unchanged.

## Code 2: predict circuit performance

### 1. Choose the run settings and run Code 2

From the project directory:

```powershell
conda activate SBMLtoODEpyWorkflow
python .\Code2_MasterSBMLPrediction_GeneralInterfaces.py --config configs/amplification_prediction.json
python .\Code2_MasterSBMLPrediction_GeneralInterfaces.py --config configs/sensitivity_prediction.json
```

Run the command for the metric you want. `--config` is required: the main
script does not assume amplification. Helper files such as
`code2_configuration.py` are imported automatically.

Code 2 reads the selected model and fitted parameters from Code 1's saved
results. It does not refit the model. If the project was copied and the
selected-model record points to a missing original folder, Code 2 tries
the corresponding local model and fitted results. `selected_model_file` in
the prediction JSON chooses the record to read. For a prefixed fit, set it
to that fit's `results/experiment_A_selected_model.json`; for a previous
launcher load, it can point to `results/loaded_model.json`.

### 2. Separate shared model defaults from run choices

The shared `configs/model_interfaces.json` contains one entry per model:

```json
{
  "model1_activegate_fuel": {
    "input_control": {"type": "species_initial_concentration", "id": "IN_temp"},
    "default_initial_species": {"RSD_temp": 25, "DRL": 500, "F_temp": 0},
    "default_observable_units": "nM"
  }
}
```

`input_control` is an optional model mapping. Amplification uses it to
translate its ON/OFF values into the appropriate species or parameter.
A metric that does not compare input conditions need not use this mapping.

The run JSON chooses what is analysed. For the sensitivity example:

```json
"model_interfaces_file": "configs/model_interfaces.json",
"observable_id": "ROL",
"analysis_control": {"type": "parameter", "id": "basal_frac"},
"initial_species_overrides": {"IN_temp": 2.5},
"fixed_parameters": {}
```

The amplification example instead varies the initial concentration of
`F_temp`. Both runs use the same shared model defaults.

The code applies initial concentrations in this order: shared defaults,
run overrides, request overrides, then the request's current control value.
An initial concentration sets the starting state; it does not clamp a species
during integration. Unspecified species retain their SBML initial values.

Code 1 supplies the baseline parameter values. Shared `default_fixed_parameters`
and run `fixed_parameters` provide explicit overrides, with run values taking
precedence. The optional secondary sweep and each request determine the
conditions simulated for that panel. Avoid assigning the same parameter
competing roles as fixed, swept, or request-controlled.

`observable_units` overrides `default_observable_units` as a display label;
it does not convert concentrations. File paths may be absolute or relative
to the directory containing Code 2. Species and parameter IDs must match
the selected SBML model; the supplied examples use `model1_activegate_fuel`.

### 3. Define the main grid and numerical settings

For example, amplification's main grid is:

```json
"control_values": {
  "values": {"start": 0, "stop": 50, "points": 51},
  "values_to_plot": [1, 3, 5, 10, 30, 50],
  "label": "Fuel template",
  "units": "nM"
}
```

`values` can also be a list such as `[0, 5, 10, 50]`. Time-series plots
choose the nearest simulated value for each `values_to_plot` entry.
Metrics returning one value per control plot the full grid.

For a parameter control, `include_fitted_value: true` adds the exact Code 1
fitted value to the grid; it must lie within the requested range. For a
per-control plot, Code 2 marks that fitted value with a vertical line when
it lies in the plotted range. This marker is separate from a metric's
reference curve or reference simulation.

`time_end_min` and `time_points` define an evenly spaced sampling grid from
zero. Alternatively, provide explicit times:

```json
"time_values_min": [0, 10, 30, 60]
```

Explicit times take precedence and must be nonnegative and strictly
increasing. `[60]` requests only the result at 60 minutes, while integration
still starts at time zero. `[0]` requests the initial concentration or rate.
Limit-only runs do not require this sampling grid.

`model_time_units_per_minute` converts displayed minutes to model time;
use `60` when the model uses seconds. `rtol` and `atol` are the ODE solver's
relative and absolute error tolerances. They apply to numerical integration
for any metric; they are not amplification settings.

Optional settings repeat the entire analysis at other parameter values:

```json
"sweep_parameter": "basal_frac",
"sweep_values": [0.01, 0.03, 0.1]
```

Use `null` for both to disable this secondary sweep. Each condition starts
from a fresh copy of the configured initial state, not the previous
control value's final state.

### 4. Follow a named request through the calculation

Every metric exposes two ordinary functions:

```python
def build_requests(experiment, settings):
    return {"output": {"property": "concentrations"}}

def calculate_metric(properties, experiment, settings):
    values = properties["output"]["values"][:, -1]
    return {
        "kind": "per_control",
        "values": values,
        "label": "Final requested concentration",
        "units": experiment["output_units"]
    }
```

In this example, `output` is the request alias and `concentrations` is the
property to calculate. The property result is returned under the same alias.
The last requested concentration is an endpoint; it is not automatically
a steady state. This example is also available in `metrics/MetricTemplate.py`.

Code 2 builds the request dictionary once for the run. A secondary sweep
changes the parameter baseline for each panel while reusing those requests.

| Request field | Meaning |
|---|---|
| `property` | `concentrations`, `rates`, `steady_state`, `output_plateau`, or a registered custom property |
| `initial_values` | Optional species initial-concentration overrides for this request |
| `parameters` | Optional parameter overrides for this request |
| `control_values` | Optional explicit grid; otherwise uses the run's main grid |

For example, a metric could request a second trajectory with
`"initial_values": {"IN_temp": 0}`. A reference request can use
`"control_values": [0]`. The main script does not assign either condition
a scientific meaning; the metric defines how to use it.

The property result is a flat dictionary:

| Property | `values` shape and meaning |
|---|---|
| `concentrations` | `(N, T)`: observable concentrations at the requested times |
| `rates` | `(N, T)`: observable ODE derivative at the requested times |
| `steady_state` | `(N,)`: verified output after every dynamic species settles |
| `output_plateau` | `(N,)`: verified output after the observable passes the extended-horizon check |

`N` is the number of control values in **that request**, and `T` is the
number of requested times. Each result also records `units` and
`control_values`; limiting properties include convergence diagnostics.
An unverified limiting value is `NaN`, with its endpoint stored separately.

Rates are concentration per **model time unit**. Multiply by
`experiment["model_time_units_per_minute"]` if rates per minute are needed.

The three dictionaries passed to `calculate_metric` have different jobs:

| Dictionary | Contents |
|---|---|
| `properties` | The calculated results, indexed by the metric's request aliases |
| `experiment` | Observable, main grid, initial state, parameters, sampling times, model mappings, and numerical settings |
| `settings` | Only the chosen metric's `metric.settings` from the run JSON |

The calculation returns `kind` (`time_series` or `per_control`), `values`
with shape `(C, T)` or `(C,)`, a vertical-axis `label`, and `units`.
Here `C` is the run's main grid length. A metric can additionally return
`reference_values`, `reference_control_value`, `reference_label`, or
`condition_label` to describe its own reference and plotted conditions.
The supplied metric files show these optional fields in use.

### 5. Amplification defines its own ON/OFF experiment

The amplification run selects:

```json
"metric": {
  "name": "amplification",
  "definition_file": "metrics/amplification.py",
  "settings": {
    "rate_floor": 1e-8,
    "input_values": {"on_value": 2.5, "off_value": 0, "label": "Input", "units": "nM"},
    "reference_value": 0
  }
}
```

Its `build_requests(...)` uses the model's input mapping and these settings
to create four requests:

| Request alias | Input condition | Control grid |
|---|---|---|
| `on_rates` | ON input | Main grid |
| `off_rates` | OFF input | Main grid |
| `reference_on` | ON input | One reference value |
| `reference_off` | OFF input | One reference value |

All four request `rates`. Code 2 calculates them, then amplification's
`calculate_metric(...)` applies the original equation:

```text
A(t,a) = [v_ON(t,a) * v_OFF(t,a_ref)] / [v_OFF(t,a) * v_ON(t,a_ref)]
```

All four absolute rates must be at least `rate_floor`, in concentration per
model time unit. Undefined values stay `NaN`. The reference line uses the
same validity rule and therefore has gaps where the ratio is undefined.

ON/OFF values and the reference value are amplification settings. A different
metric is free to request one condition, several conditions, or different
control grids using the same general engine.

### 6. Sensitivity requests one limiting response

The supplied example chooses `ROL`, varies `basal_frac`, and explicitly
sets `IN_temp` to 2.5 nM in `initial_species_overrides`. It does not request
ON/OFF comparisons or reference simulations.

```json
"metric": {
  "name": "steady_state_sensitivity",
  "definition_file": "metrics/sensitivity.py",
  "settings": {
    "response_property": "output_plateau",
    "concentration_floor": 1e-8,
    "absolute": false
  }
}
```

`build_requests(...)` returns one alias, `response`, requesting the selected
limiting property over the main grid. The formula then reads
`properties["response"]["values"]` and its parameter grid:

```text
S(theta) = (d xss / d theta) * (theta / xss)
```

For an output plateau, `xss` in this expression is the accepted plateau value.
The plot label records which meaning was requested. The metric uses
`np.gradient(xss, theta, edge_order=2)`, followed by multiplication by
`theta / xss`. The result is signed and dimensionless unless `absolute` is
enabled.

The parameter grid needs at least three distinct, positive, increasing
values. A list or a range is supported, including logarithmic spacing:

```json
"values": {"start": 0.005, "stop": 0.05, "points": 41, "spacing": "log"}
```

The derivative is taken against physical parameter values. Failed or
too-small responses become `NaN`, as do derivatives whose three-point
stencil touches an invalid response. Refine the grid and tighten numerical
tolerances when checking accuracy; tiny parameter spacing can amplify
simulation errors. Code 2 adds the vertical marker for Code 1's fitted
parameter value; sensitivity does not return a reference line or require
a `reference_value` setting.

**Choose the convergence criterion:**

- `steady_state`: every dynamic species must have sufficiently small rates
  and concentration changes for consecutive check windows.
- `output_plateau`: only the observable must pass those local checks, then
  continue passing them over an extended simulation horizon. Other dynamic
  species can still change.

Both integrate the complete coupled ODE system. The property name selects
the criterion; these general numerical settings control the checks:

```json
"limit_settings": {
  "max_time_min": 5000,
  "check_interval_min": 100,
  "consecutive_checks": 3,
  "samples_per_check": 9,
  "plateau_extension_factor": 2,
  "rate_tolerance": 1e-8,
  "change_atol": 1e-6,
  "change_rtol": 1e-6
}
```

| Setting | Meaning |
|---|---|
| `max_time_min` | Total integration cap, including plateau verification, in minutes |
| `check_interval_min` | Length of each check window, in minutes |
| `consecutive_checks` | Required consecutive passing windows, at least 2 |
| `samples_per_check` | Sampled points per window, at least 3 |
| `plateau_extension_factor` | Multiplier for the candidate plateau time, greater than 1 |
| `rate_tolerance` | Absolute rate limit, concentration per model time unit |
| `change_atol`, `change_rtol` | Allowed concentration range: `change_atol + change_rtol * max(abs(concentration))` |

A plateau candidate at 600 minutes with factor 2 must remain valid through
at least 1200 minutes. Verification also requires at least
`consecutive_checks` extra full windows. The cumulative sampled output range
must remain within tolerance throughout the extension, which catches gradual
drift and excursions even if the endpoint returns to its earlier value.
A failed candidate is discarded, and integration can continue looking for
another before the time cap.

Nonfinite or significantly negative concentrations anywhere in the state
invalidate either mode. A finite-horizon check provides numerical evidence
for a persistent limit; it does not prove an infinite-time limit or find
every steady-state branch. Increasing the cap and extension factor checks
whether the result persists longer. The limiting solver converts its
concentration `atol` to species amounts using compartment volumes.

**The supplied ROL / basal_frac example is a useful zero-sensitivity check.**
With positive RSD_temp and k_txn, `model1_activegate_fuel` has:

```text
d(uRSDg + RSDg + I_RSDg + F_RSDg)/dt = k_txn * RSD_temp > 0
```

The full system cannot reach a finite steady state under these conditions.
ROL can nevertheless plateau near 500 nM as reporter is exhausted. Its
plateau sensitivity to basal_frac is then approximately zero. This describes
the limiting output, not its transient response or a full-system equilibrium.

### 7. Add another metric or property

For a new metric, copy `metrics/MetricTemplate.py` and edit its two functions:

1. `build_requests(...)` specifies the required properties and conditions.
2. `calculate_metric(...)` combines the returned arrays and supplies labels.

Point the run JSON's `metric.definition_file` to the new file. Select species,
grids, and initial conditions in the run JSON; put formula-specific options
under `metric.settings`.

If the built-in properties do not supply a required quantity, copy
`properties/PropertyTemplate.py`. Its function has this interface:

```python
def calculate_property(simulator, experiment, parameters, settings):
    return {"values": values, "units": experiment["output_units"]}
```

The property receives the request's conditions and control grid. Return a
flat dictionary with one row per requested control, a string `units`, and
any diagnostic arrays. Use `NaN` for undefined values. Fields need consistent
shapes between secondary-sweep panels so they can be saved together.

Register the file in the run JSON:

```json
"property_definitions": {
  "initial_output": {
    "definition_file": "properties/PropertyTemplate.py",
    "settings": {}
  }
}
```

Then a metric can request
`{"initial": {"property": "initial_output"}}` and read
`properties["initial"]["values"]`. Registration alone does not run the
property. Built-in names cannot be overwritten. If a custom calculation
needs several steps, write those steps in its function or ordinary helper
functions; registering properties does not create dependencies between them.

### 8. Read the saved results

By default, Code 2 saves the figure in `figures/` and its numerical results
in `results/`. It names them `<prefix><model>_<metric>`, where `<metric>` is
the metric Python filename without `.py`. The supplied examples therefore
produce `<model>_amplification` and `<model>_sensitivity`.

When running Code 2 by itself, its prediction JSON controls `figures_dir`,
`results_dir` and the optional `output_prefix`. When using the launcher,
these shared folders and the prefix come from the fitting JSON.

An explicit `output_file`, such as `figures/my_amplification.png`, chooses
the figure path and basename instead. The launcher keeps its basename and
puts the figure in the shared `figures_dir`. The numerical files still go
to `results_dir`, using that same basename:

| Output | Contents |
|---|---|
| `.png` in `figures/` | Performance plot, unless `output_file` chooses another location |
| `.npz` in `results/` | Metric values, property arrays, grids, request names, labels, and units |
| `.settings.json` | The run configuration |
| `.model_interface.json` | Resolved model interface and exact named property requests |
| `.csv` | Convergence diagnostics when one limiting request is present |
| `.<alias>.csv` | Separate diagnostics for each alias when several limiting requests are present |

All non-figure files in this table go to `results_dir`. The launcher also
saves `<model>_<metric>.log` there, with the same optional prefix or chosen
basename. This log contains Code 2's detailed messages.

The launcher temporarily passes Code 2 a `.pending_settings.json` file.
After prediction succeeds, the final `.settings.json` is saved and the pending
file is removed. If prediction stops, the pending file and log show what was
attempted; the saved settings beside an older result stay unchanged.

The NPZ `metric_values` shape is `(panels, controls, times)` for a time-series
metric or `(panels, controls)` for one value per control. Each property field
is saved under `property_<alias>_<field>` with a leading panel dimension.
For example, `property_on_rates_values` stores amplification's ON rates;
`property_response_values` stores sensitivity's accepted limiting output.

`request_names` records the aliases and `requested_properties` records their
property types in the same order. Optional metric references are recorded
as `reference_values_0`, `reference_control_value_0`, and so on for each panel
that returns them. There are no duplicate legacy `steady_state_*` aliases.

Limiting-property tables record their actual request grid, observable,
verified values, endpoint values, convergence flags, stopping time, and
diagnostics. The value column is `xss` for full-system convergence or
`x_plateau` for the output-plateau criterion. A scalar metric value is included
when it corresponds to the request's grid; a time-series metric has no single
scalar for that column.

Useful limiting-property fields include:

| Field | Meaning |
|---|---|
| `values`, `converged`, `status` | Accepted output or NaN, acceptance flag, and reason |
| `convergence_mode` | `full_system` or `output_plateau` |
| `time_minutes`, `final_output` | Stopping time and final simulated output, including failed runs |
| `max_abs_rate`, `max_scaled_change` | Whole-system diagnostics |
| `output_max_abs_rate`, `output_max_scaled_change` | Observable diagnostics |
| `full_system_converged` | Whether the whole system also passed its criterion |
| `plateau_start_time_minutes`, `plateau_target_time_minutes` | Candidate and confirmation times |
| `plateau_max_scaled_change` | Cumulative output variation during extension relative to tolerance |
| `worst_rate_species` | Species with the largest absolute rate |
| `final_concentrations`, `state_species_ids` | Full-state endpoint and species order |

The endpoint is never substituted for an unverified limiting value.
`output_plateau_verified` indicates a verified plateau;
`plateau_extension_incomplete` means the time cap prevented confirmation.
The sensitivity derivative `dxss_dtheta` is also returned in memory by the
metric; extra metric fields are not automatically saved to NPZ.

Use a different prefix or an explicit output filename to keep several analyses.
Running again with the same filename replaces those outputs. To reproduce a run, retain its
settings and resolved-interface snapshot along with the SBML, metric files,
and fitted-parameter results.

### Updating older prediction files

Both supplied configurations and metrics are already updated.

The property helpers now live in `properties/`. Custom files that import
`build_condition` should use `from properties.model_properties import build_condition`.
Imports of the limiting-value solver should use `properties.steady_state_simulation`.
This folder change leaves run commands, JSON settings, metric requests,
and saved output formats unchanged.

| Previous arrangement | Current arrangement |
|---|---|
| Amplification-specific Code 2 filename | `Code2_MasterSBMLPrediction_GeneralInterfaces.py` |
| Running without a configuration argument | Required `--config <run-settings.json>` |
| Run-level `input_values` | Amplification's `metric.settings.input_values`; other metrics set their own conditions |
| `control_values.reference_value` | Amplification's `metric.settings.reference_value`; Code 2 marks the fitted parameter directly on sensitivity plots |
| Sensitivity's `limit_settings.input_condition` | Explicit starting conditions, e.g. `initial_species_overrides.IN_temp: 2.5` |
| `REQUIRED_PROPERTIES` or `get_required_properties(settings)` | `build_requests(experiment, settings)` returning named requests |
| `properties["rates"]["on"]` | `properties["on_rates"]["values"]` |
| `properties["rates"]["reference_on"]` | `properties["reference_on"]["values"][0]` |
| `properties[response_property]["values"]` | `properties["response"]["values"]` |
| Property NPZ fields named by type | Fields named by request alias |
| Duplicate `steady_state_*` NPZ fields | Use `property_<alias>_<field>` |

Earlier migrations still apply: shared defaults belong in the separate
model-interface file; observable and analysis control belong in the run JSON;
`enabled` and `simulation_mode` are no longer used. Old `steady_state` settings
belong under `limit_settings` without `convergence_mode`. Sensitivity chooses
the criterion through `metric.settings.response_property`.

### Check the calculations

```powershell
python -m unittest discover -s tests -v
```

Tests cover amplification and sensitivity formulas, request conditions and
grids, property sampling, convergence checks, imported properties, and model
configuration. The supplied amplification calculation retains the original
four-rate formula and rate-floor rules.

## To perform a leave-one-out (LOO) model validation

Run `Run_SAMPL.py` and answer **yes** to LOO. It withholds one complete
experimental condition at a time by changing `fit_conditions` in separate
settings copies. It never deletes Excel rows or columns.

Only the selected model is refitted, including when the launcher is in load
mode. LOO uses the current fitting JSON's `all_conditions`, Excel data and
parameter settings. With the current six conditions, there are six fits,
each using five conditions and predicting the sixth.

The model choice stays fixed throughout LOO. These scores evaluate refitting
that selected model; they do not cross-validate the complete model-selection
procedure. Prediction afterwards still uses the identification fit or loaded
fit, not a fold's fitted parameters.

Read `results/<model>_loo_summary.csv` or `results/<model>_loo_summary.json`
for each condition's prediction error. No average across folds is calculated.
Each fold saves its own plots, parameters, settings
and optimizer status with a prefix such as `loo_01_`. An `output_prefix`
is added before these names, for example `experiment_A_loo_01_`.

To set up individual splits manually with Code 1 instead:

1. List every experimental condition in `all_conditions`.
2. Add an explicit `fit_conditions` list for each run, leaving out one condition while keeping it in `all_conditions`.
3. Use `--mode fit` so the parameters are estimated using only the remaining conditions.
4. Give each split a different `output_prefix`, such as `loo_without_condition1_`, to keep its files in the shared folders. Set `save_model_selection` to `false` to save the fit and its scores without writing another selection record.
5. Run Code 1 with `--config <split settings>` and `--model <selected SBML file>`, then inspect its held-out prediction plot and `<prefix><model>_held_out_rmse.csv`. When the file is a saved `selected_model.xml`, also give `--model-name <original model name>` so its model-specific settings are retained.
6. Repeat until every condition has been held out once.


## Advanced settings and troubleshooting

- Use `--mode fit` after changing data, models, parameters, or fitting conditions.
- Use `--mode load` only for results from the exact same setup.
- Keep the current Excel-layout, cutoff, and solver settings unless the data layout or units change.
- A held-out plot exists only when at least one condition is outside `fit_conditions` and the model is refitted.
- After an XML changes, refit with Code 1. The adapter checks its contents and refreshes the matching generated Python model automatically.
- A parameter at its fitting bound may be weakly identified or may require scientifically justified bounds.
- An undefined amplification value means an ON or OFF rate fell below `metric.settings.rate_floor`.
- A constant-species warning means the workflow attempted to change a species marked constant in SBML; check whether the configured value already matches the SBML value.
- The adapter supports the SBML features used here, but not arbitrary assignment rules, rate rules, algebraic rules, or events.
