# SBML Fitting and Amplification Prediction

This project fits candidate SBML models to time-course data, compares their predictions, and uses the selected fitted model to predict amplification.

## Contents

- [Workflow](#workflow)
- [Important terms](#important-terms)
- [Data file structure](#data-file-structure)
- [Code 1: fit and compare models](#code-1-fit-and-compare-models)
- [Code 2: predict amplification](#code-2-predict-amplification)
- [To perform a leave-one-out (LOO) model validation](#to-perform-a-leave-one-out-loo-model-validation)
- [Advanced settings and troubleshooting](#advanced-settings-and-troubleshooting)

## Workflow

```text
Data + candidate SBML models
              |
              v
   Code 1: fit and compare
              |
 selected model + fitted parameters
              |
              v
 Code 2: predict amplification
```

| Needed item | Purpose |
|---|---|
| `models/*.xml` | Candidate reaction mechanisms. |
| `data/` | Experimental time courses. |
| [`configs/fit_all_models.json`](configs/fit_all_models.json) | Settings for Code 1. |
| [`Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py`](Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py) | Fits and compares the models. |
| `results/selected_model.json` | Passes the selected model and fitted parameters to Code 2. |
| [`configs/amplification_prediction_general.json`](configs/amplification_prediction_general.json) | Settings for Code 2. |
| [`Code2_MasterSBMLAmplificationPrediction_GeneralInterfaces.py`](Code2_MasterSBMLAmplificationPrediction_GeneralInterfaces.py) | Predicts amplification. |
| `sbmltoodepy_solveivp_adapter.py` | Simulates the SBML equations. |

Set up the environment once:

```powershell
conda env create -f environment.yml
conda activate SBMLtoODEpyWorkflow
```

## Important terms

| Term | Meaning |
|---|---|
| Condition | One complete experimental curve containing many time points. |
| Fitting condition | A curve used to estimate parameters. |
| Held-out condition | A complete curve predicted without using it during fitting. |
| `initial_values` | SBML species values that Code 1 sets at time zero. These are not automatically the input or observable. |
| Shared species name | A model-independent experimental role used in the conditions. |
| Species alias | Translation from a shared name to an exact model-specific SBML ID. |
| Observable | The simulated species compared with the measured data. |
| ON/OFF input | The two input values compared by Code 2. |
| Amplification reference | The amplification-control value used as the unamplified baseline. It is different from the OFF input. |

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
  "fit_conditions": ["condition_high"],
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
  "selection_metric": "held_out_normalized_rmse",
  "n_starts": 8,
  "random_seed": 20260806,
  "max_nfev": 400
}
```

| Setting | Meaning |
|---|---|
| `selection_metric` | Use `held_out_normalized_rmse` for held-out prediction ranking or `global_normalized_rmse` for fitted-curve ranking. |
| `n_starts` | Number of starting parameter sets tried per model. More starts may find a better solution but take longer. |
| `random_seed` | Makes the starting sets reproducible. Keep it fixed when comparing models. |
| `max_nfev` | Maximum model evaluations per start. Larger values allow more fitting work but take longer. |

A condition in `all_conditions` but not `fit_conditions` is held out. One run evaluates one split. Full leave-one-out validation requires repeated runs with a different condition omitted each time.

### 5. Run Code 1

Set `MODE = "fit"` near the top of Code 1, then run:

```powershell
python .\Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py
```

Check these outputs:

- `results/model_comparison.csv`: model scores and ranks;
- `results/selected_model.json`: selected model and fitted-results path;
- `results/<model>/fitted_parameters.csv`: readable parameter values;
- `results/<model>/held_out_rmse.csv`: held-out error;
- `figures/<model>/<model>_fit.png`: fitted curves; and
- `figures/<model>/<model>_held_out_prediction.png`: held-out curves.

The machine-readable fitted parameters are stored in `results/<model>/<model>_results.npz` and loaded by Code 2.

## Code 2: predict amplification

Code 2 simulates the selected fitted model at ON and OFF input values. It repeats these simulations across amplification-control values and compares them with the unamplified reference.

Edit [`configs/amplification_prediction_general.json`](configs/amplification_prediction_general.json).

### 1. Describe each model's interface

```json
{
  "enabled": true,
  "selected_model_file": "results/selected_model.json",
  "model_interfaces": {
    "example_amplifier": {
      "input_control": {
        "type": "species_initial_concentration",
        "id": "input_species"
      },
      "observable_id": "output_species",
      "amplification_control": {
        "type": "species_initial_concentration",
        "id": "amplification_species"
      },
      "fixed_initial_species": {
        "fixed_species": 500.0
      },
      "fixed_parameters": {}
    }
  }
}
```

Set `enabled` to `true` to run Code 2. 'enabled' is a safety switch if you wish to automatically run Code 2 after Code 1.
Each `model_interfaces` key must match an XML filename without `.xml`.

| Field | Meaning |
|---|---|
| `input_control` | Species or parameter changed between ON and OFF. |
| `observable_id` | Output species whose rate is measured. |
| `amplification_control` | Species or parameter varied to test amplification. |
| `fixed_initial_species` | Model-specific species values set at time zero. |
| `fixed_parameters` | Optional parameter overrides. |

The control `type` may be `species_initial_concentration` or `parameter`. All `id` values must exactly match the selected model's SBML IDs.

Use `"amplification_control": null` for a model without an amplification variable. Code 1 can fit that model, but Code 2 cannot predict amplification for it.

### 2. Define the prediction experiment

```json
{
  "input_values": {
    "on_value": 50.0,
    "off_value": 0.0,
    "label": "Input",
    "units": "nM"
  },
  "amplification_values": {
    "reference_value": 0.0,
    "values": {
      "start": 0.0,
      "stop": 50.0,
      "points": 51
    },
    "values_to_plot": [1.0, 5.0, 10.0, 50.0],
    "label": "Amplification control",
    "units": "nM"
  },
  "time_end_min": 200.0,
  "time_points": 1000,
  "model_time_units_per_minute": 60.0,
  "output_file": "figures/amplification_prediction.png"
}
```

- `on_value` and `off_value` define the two input simulations.
- `reference_value` is the baseline amplification-control value; it is not the OFF input.
- `values` defines all amplification values to simulate.
- `values_to_plot` selects the displayed curves.
- `model_time_units_per_minute` converts minutes into the SBML time units. Use `60.0` for a model expressed in seconds.

`values` may also be a direct list such as `[0.0, 5.0, 10.0, 50.0]`.

### 3. Optionally sweep another parameter

```json
{
  "sweep_parameter": null,
  "sweep_values": null
}
```

JSON `null` disables the optional secondary parameter sweep. To vary another model parameter during prediction without refitting, set `sweep_parameter` to its exact SBML parameter ID and provide a list of values to sweep through in `sweep_values`.
For example:

```json
{
"sweep_parameter": "leak_rate",
"sweep_values": [0.000001, 0.00001, 0.0001]
}
```

### 4. Run Code 2

```powershell
python .\Code2_MasterSBMLAmplificationPrediction_GeneralInterfaces.py
```

The figure is saved to the configured `output_file`, currently `figures/amplification_prediction_general_interface.png`.

## To perform a leave-one-out (LOO) model validation

1. List every experimental condition in `all_conditions`.
2. For each run, omit one condition from `fit_conditions` while keeping it in `all_conditions`.
3. Set `MODE = "fit"` so the parameters are estimated using only the remaining conditions.
4. Give each run different `results_dir` and `figures_dir` values so its results are not overwritten. A different `run_name` may also be used to identify the run.
5. Run Code 1 and inspect the held-out prediction plot and `held_out_rmse.csv`.
6. Repeat until every condition has been held out once.


## Advanced settings and troubleshooting

- Use `MODE = "fit"` after changing data, models, parameters, or fitting conditions.
- Use `MODE = "load"` only for results from the exact same setup.
- Keep the current Excel-layout, cutoff, and solver settings unless the data layout or units change.
- A held-out plot exists only when at least one condition is outside `fit_conditions` and the model is refitted.
- If an XML changes, delete only its matching file from `generated_models/` and rerun Code 1.
- A parameter at its fitting bound may be weakly identified or may require scientifically justified bounds.
- An undefined amplification value means an ON or OFF rate fell below `rate_floor`.
- A constant-species warning means the workflow attempted to change a species marked constant in SBML; check whether the configured value already matches the SBML value.
- The adapter supports the SBML features used here, but not arbitrary assignment rules, rate rules, algebraic rules, or events.
