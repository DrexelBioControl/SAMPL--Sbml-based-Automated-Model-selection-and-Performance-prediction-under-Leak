# ============================================================
# MASTER SBML AMPLIFICATION PREDICTION
# ============================================================
#
# PURPOSE
# -------
# Predict how an amplification control changes ON/OFF output
# discrimination for the model selected by Code 1.
#
# Code 2 does not fit parameters or compare models.
# It uses the selected model and fitted parameters produced by Code 1.
#
# WORKFLOW
# --------
# Code 1 produces:
#
#     selected model
#           +
#     fitted parameter values
#           ↓
#     selected_model.json
#
# Code 2 then combines:
#
#     selected model and fitted parameters
#                  +
#     model-specific simulation interface
#                  +
#     amplification-control values
#                  ↓
#     predicted amplification time courses
#
# MODEL INTERFACES
# ----------------
# Each candidate model has an interface in the JSON configuration.
#
# The interface identifies:
#
#     input_control
#     observable_id
#     amplification_control
#     optional fixed initial species
#     optional fixed parameters
#
# Models may use completely different SBML species and parameter IDs.
# Gates, reporters, and other fixed species are included only when required.
#
# A fitting model may define amplification_control as null.
# Such a model can be used by Code 1 but cannot be simulated by Code 2.
#
# INPUT AND AMPLIFICATION CONTROLS
# --------------------------------
# An input or amplification control may be:
#
#     a species initial concentration
#     or
#     an SBML model parameter
#
# The JSON specifies the control type and its model-specific SBML ID.
#
# AMPLIFICATION DEFINITION
# ------------------------
# For each amplification-control value a, Code 2 calculates:
#
#                   v_ON(t,a) * v_OFF(t,a_ref)
#     A(t,a) = -----------------------------------
#                   v_OFF(t,a) * v_ON(t,a_ref)
#
# a_ref is the configured unamplified reference value.
#
# A value of 1 means no change relative to the reference condition.
# A value above 1 means improved ON/OFF rate discrimination.
# A value below 1 means reduced ON/OFF rate discrimination.
#
# Relevant output rates must exceed rate_floor.
# Amplification is undefined where a required rate is too small.
#
# OPTIONAL PARAMETER SWEEP
# ------------------------
# Code 2 can optionally vary one additional model parameter.
#
# All other parameters remain at their fitted values.
# This is a sensitivity analysis and does not refit the model.
#
# Set sweep_parameter and sweep_values to null to disable this feature.
#
# WHAT THE SCRIPT DOES
# --------------------
# 1. Loads the model-selection record produced by Code 1.
# 2. Loads the selected SBML model.
# 3. Restores its fitted parameter values.
# 4. Selects the model-specific simulation interface.
# 5. Validates all required species and parameter IDs.
# 6. Simulates the ON and OFF reference conditions.
# 7. Simulates every configured amplification-control value.
# 8. Calculates normalized amplification time courses.
# 9. Optionally repeats the analysis for a parameter sweep.
# 10. Saves and displays the amplification figure.
#
# USER SETTINGS
# -------------
# User-configurable settings are stored in:
#
#     configs/amplification_prediction_general.json
#
# Add a model_interfaces entry when supporting a new SBML model.
#
# ============================================================

# Defer evaluation of type annotations for compatibility.
from __future__ import annotations

# Read JSON workflow configuration files.
import json
# Build portable paths relative to this script.
from pathlib import Path

# Create and save amplification plots.
import matplotlib.pyplot as plt
# Perform numerical array and grid calculations.
import numpy as np

# Simulate the selected SBML model with the solve_ivp adapter.
from sbmltoodepy_solveivp_adapter import SBMLtoODEpySolveIVPModel


# ============================================================
# 0) PATHS AND CONFIGURATION
# ============================================================

# Locate the project directory containing this script.
BASE_DIR = Path(__file__).resolve().parent
# Locate the amplification-prediction configuration file.
CONFIG_FILE = BASE_DIR / "configs" / "amplification_prediction_general.json"

# Load the amplification-prediction settings.
with open(CONFIG_FILE, "r", encoding="utf-8") as f:
    # Convert the JSON settings into a Python dictionary.
    config = json.load(f)

# Exit cleanly when this workflow is disabled in the configuration.
if not bool(config.get("enabled", False)):
    # Tell the user which configuration disabled the workflow.
    print(f"Amplification prediction is disabled in:\n{CONFIG_FILE}")
    # Stop before loading or simulating a model.
    raise SystemExit(0)


# Convert an absolute or project-relative value into a usable path.
def resolve_path(path_value):
    """Resolve absolute or project-relative paths."""
    # Interpret the configured value as a filesystem path.
    path = Path(path_value)
    # Keep absolute paths or anchor relative paths to the project directory.
    return path if path.is_absolute() else BASE_DIR / path


# Convert a configured list or range object into a numeric array.
def parse_values(value_spec):
    """Accept either a list of values or a start/stop/points range object."""
    # Expand range dictionaries into evenly spaced values.
    if isinstance(value_spec, dict):
        # Generate the requested start-to-stop sampling grid.
        return np.linspace(
            float(value_spec["start"]),
            float(value_spec["stop"]),
            int(value_spec["points"]),
        )

    # Convert explicit value lists directly to floating-point arrays.
    return np.asarray(value_spec, dtype=float)


# ============================================================
# 1) LOAD SELECTED MODEL AND BEST-FIT PARAMETERS
# ============================================================

# Locate the model-selection record produced by Part 1.
SELECTED_MODEL_FILE = resolve_path(config["selected_model_file"])

# Require the model-selection record produced by Part 1.
if not SELECTED_MODEL_FILE.exists():
    # Explain how to generate the missing selection record.
    raise FileNotFoundError(
        "Part 1 selection record was not found:\n"
        f"{SELECTED_MODEL_FILE}\n\n"
        "Run Code1_MasterSBMLTimecourseFitting first."
    )

# Load the selected mechanism and its fitted-parameter files.
with open(SELECTED_MODEL_FILE, "r", encoding="utf-8") as f:
    # Convert the selection record into a Python dictionary.
    selected = json.load(f)

# Get the selected mechanism name for reporting and plot labels.
MODEL_NAME = selected["model_name"]
# Resolve the selected SBML model path.
MODEL_FILE = resolve_path(selected["model_file"])
# Resolve the saved best-fit parameter path.
FIT_RESULTS_FILE = resolve_path(selected["fit_results_file"])

# Use the archived selected model if its original SBML path is unavailable.
if not MODEL_FILE.exists():
    # Locate the archived copy of the selected SBML model.
    fallback_model = BASE_DIR / "results" / "selected_model.xml"
    # Substitute the archived model when it is available.
    if fallback_model.exists():
        # Use the archived model as the simulation source.
        MODEL_FILE = fallback_model
    else:
        # Stop when neither selected-model path is valid.
        raise FileNotFoundError(f"Selected SBML file was not found:\n{MODEL_FILE}")

# Require the best-fit parameter file produced by Part 1.
if not FIT_RESULTS_FILE.exists():
    # Stop because simulations require the fitted parameter values.
    raise FileNotFoundError(
        f"Selected fit-results file was not found:\n{FIT_RESULTS_FILE}"
    )

# Restore the best-fit parameters from Part 1.
saved_fit = np.load(FIT_RESULTS_FILE, allow_pickle=True)

# Extract the fitted values stored in log10 space.
best_x = np.asarray(saved_fit["best_x"], dtype=float)
# Extract the parameter names in their fitted order.
fit_params = [str(name) for name in saved_fit["fit_params"]]

# Part 1 stores fitted kinetic parameters in log10 space.
# Convert fitted log-space parameters back to physical values.
base_params = {
    parameter_name: float(10 ** best_x[index])
    for index, parameter_name in enumerate(fit_params)
}

# Locate the folder used for generated simulator modules.
GENERATED_MODELS_DIR = resolve_path(
    config.get("generated_models_dir", "generated_models")
)

# Initialize a simulator for the selected SBML mechanism.
rr = SBMLtoODEpySolveIVPModel(
    sbml_file=MODEL_FILE,
    generated_models_dir=GENERATED_MODELS_DIR,
)

# Report the selected mechanism.
print(f"Loaded selected model : {MODEL_NAME}")
# Report the SBML file used for simulation.
print(f"Loaded SBML           : {MODEL_FILE}")
# Report the best-fit parameter file.
print(f"Loaded fit            : {FIT_RESULTS_FILE}")

# Introduce the fitted-parameter summary.
print("\nBest-fit parameters:")
# Print each fitted parameter in scientific notation.
for parameter_name, value in base_params.items():
    # Display the current parameter name and physical value.
    print(f"{parameter_name:<20} = {value:.8e}")


# ============================================================
# 2) SELECT THE MODEL-SPECIFIC SIMULATION INTERFACE
# ============================================================

# Get the per-model interface definitions from the configuration file.
MODEL_INTERFACES = config["model_interfaces"]

# Require an interface for the model selected by Code 1.
if MODEL_NAME not in MODEL_INTERFACES:
    # List the configured alternatives to make spelling errors easy to find.
    raise KeyError(
        f"No model interface was configured for {MODEL_NAME!r}.\n"
        f"Configured models: {sorted(MODEL_INTERFACES)}"
    )

# Use only the interface belonging to the selected model.
MODEL_INTERFACE = MODEL_INTERFACES[MODEL_NAME]

# Require the selected model to define its ON/OFF input control.
if "input_control" not in MODEL_INTERFACE:
    # Stop because Code 2 cannot construct the ON and OFF conditions.
    raise KeyError(
        f"The interface for {MODEL_NAME!r} has no input_control."
    )

# Get the selected model's input-control target.
INPUT_TARGET = MODEL_INTERFACE["input_control"]

# Require an amplification control for Code 2.
AMPLIFICATION_TARGET = MODEL_INTERFACE.get("amplification_control")

# Stop clearly when a fitting model has no amplification mechanism.
if AMPLIFICATION_TARGET is None:
    # Explain that the model may still be valid for Code 1.
    raise ValueError(
        f"The selected model {MODEL_NAME!r} has no amplification_control. "
        "It can be fitted by Code 1 but cannot be used for Code 2."
    )

# Define the supported ways to control a simulation condition.
SUPPORTED_CONTROL_TYPES = {
    "species_initial_concentration",  # <--------- choice for value to enter in
    "parameter",                      # <--------- json config file row "type:"
}

# Get whether the ON/OFF input changes a species or parameter.
INPUT_TYPE = str(INPUT_TARGET["type"])
# Get the selected model's SBML ID for the ON/OFF input.
INPUT_ID = str(INPUT_TARGET["id"])

# Reject unsupported input-control types.
if INPUT_TYPE not in SUPPORTED_CONTROL_TYPES:
    # Report the invalid value and the supported choices.
    raise ValueError(
        f"Unsupported input_control type {INPUT_TYPE!r}. "
        f"Choose from {sorted(SUPPORTED_CONTROL_TYPES)}."
    )

# Get whether amplification changes a species or parameter.
AMP_TYPE = str(AMPLIFICATION_TARGET["type"])
# Get the selected model's SBML ID for the amplification control.
AMP_ID = str(AMPLIFICATION_TARGET["id"])

# Reject unsupported amplification-control types.
if AMP_TYPE not in SUPPORTED_CONTROL_TYPES:
    # Report the invalid value and the supported choices.
    raise ValueError(
        f"Unsupported amplification_control type {AMP_TYPE!r}. "
        f"Choose from {sorted(SUPPORTED_CONTROL_TYPES)}."
    )

# Get the SBML species used as the measured output.
OBSERVABLE_ID = str(MODEL_INTERFACE["observable_id"])

# Get species concentrations held fixed for the selected model.
FIXED_INITIAL_SPECIES = {
    str(species_id): float(value)
    for species_id, value in MODEL_INTERFACE.get(
        "fixed_initial_species",
        {},
    ).items()
}

# Get model parameters held fixed during every Code 2 simulation.
FIXED_PARAMETERS = {
    str(parameter_id): float(value)
    for parameter_id, value in MODEL_INTERFACE.get(
        "fixed_parameters",
        {},
    ).items()
}

# Get the shared ON/OFF input values and display settings.
INPUT_VALUES = config["input_values"]
# Get the configured input value for the ON condition.
INPUT_ON = float(INPUT_VALUES["on_value"])
# Get the configured input value for the OFF condition.
INPUT_OFF = float(INPUT_VALUES["off_value"])

# Get the amplification values and display settings.
AMPLIFICATION_VALUES = config["amplification_values"]
# Get the unamplified reference value.
AMP_REFERENCE = float(AMPLIFICATION_VALUES["reference_value"])
# Get all amplification-control values to simulate.
AMP_VALUES = parse_values(AMPLIFICATION_VALUES["values"])

# Get the requested subset of amplification values to plot.
AMP_VALUES_TO_PLOT = np.asarray(
    AMPLIFICATION_VALUES.get(
        "values_to_plot",
        AMP_VALUES,
    ),
    dtype=float,
)

# Get the display label for the amplification control.
AMP_LABEL = str(AMPLIFICATION_VALUES.get("label", AMP_ID))
# Get the display units for the amplification control.
AMP_UNITS = str(AMPLIFICATION_VALUES.get("units", ""))


# ============================================================
# 3) TIME AND NUMERICAL SETTINGS
# ============================================================

# Get the simulation end time in minutes.
T_END_MIN = float(config["time_end_min"])
# Get the number of requested simulation time points.
N_TIMEPOINTS = int(config["time_points"])

# Get the conversion from minutes to the model's time units.
MODEL_TIME_UNITS_PER_MINUTE = float(
    config.get("model_time_units_per_minute", 60.0)
)

# Create the output time grid in minutes.
t_minutes = np.linspace(0.0, T_END_MIN, N_TIMEPOINTS)
# Convert the output grid to the model's time units.
t_model = t_minutes * MODEL_TIME_UNITS_PER_MINUTE

# Get the minimum rate magnitude permitted in amplification ratios.
RATE_FLOOR = float(config.get("rate_floor", 1e-8))

# Get the relative integration tolerance.
RTOL = float(config.get("rtol", 1e-6))
# Get the absolute integration tolerance.
ATOL = float(config.get("atol", 1e-6))

# Get the optional secondary parameter to sweep.
SWEEP_PARAM = config.get("sweep_parameter")
# Get the optional secondary sweep values.
SWEEP_VALUES = config.get("sweep_values")

# Convert configured sweep values to a numeric array when provided.
if SWEEP_VALUES is not None:
    # Store sweep values as floating-point numbers.
    SWEEP_VALUES = np.asarray(SWEEP_VALUES, dtype=float)


# ============================================================
# 4) VERIFY REQUIRED SBML IDS
# ============================================================

# Collect the species IDs available in the selected model.
species_ids = set(rr.species_ids)
# Collect the global parameter IDs available in the selected model.
parameter_ids = set(rr.parameter_ids)

# Start the required-species set with all fixed initial species.
required_species = set(FIXED_INITIAL_SPECIES)

# Require the ON/OFF input species when the input changes a species.
if INPUT_TYPE == "species_initial_concentration":
    # Add the input-control species to the validation set.
    required_species.add(INPUT_ID)

# Require the amplification-control species when amplification changes its initial value.
if AMP_TYPE == "species_initial_concentration":
    # Add the amplification-control species to the validation set.
    required_species.add(AMP_ID)

# Identify required species missing from the selected model.
missing_species = sorted(required_species.difference(species_ids))

# Report all missing required species in one error.
if missing_species:
    # List every missing species in the error message.
    raise ValueError(
        "The selected SBML model is missing required species:\n"
        + "\n".join(f"  - {species_id}" for species_id in missing_species)
    )

# Require the configured observable species.
if OBSERVABLE_ID not in species_ids:
    # Stop when the model cannot provide the requested output.
    raise ValueError(
        f"Observable {OBSERVABLE_ID!r} was not found in the selected SBML model."
    )

# Require the amplification parameter when amplification changes a parameter.
if AMP_TYPE == "parameter" and AMP_ID not in parameter_ids:
    # Report the missing amplification parameter and available alternatives.
    raise ValueError(
        f"Amplification parameter {AMP_ID!r} was not found.\n"
        f"Available global parameters: {sorted(parameter_ids)}"
    )

# Require the input parameter when the ON/OFF input changes a parameter.
if INPUT_TYPE == "parameter" and INPUT_ID not in parameter_ids:
    # Report the missing input parameter and available alternatives.
    raise ValueError(
        f"Input parameter {INPUT_ID!r} was not found.\n"
        f"Available global parameters: {sorted(parameter_ids)}"
    )

# Identify fixed parameters missing from the selected model.
missing_fixed_parameters = sorted(
    set(FIXED_PARAMETERS).difference(parameter_ids)
)

# Report all missing fixed parameters in one error.
if missing_fixed_parameters:
    # List every missing parameter in the error message.
    raise ValueError(
        "The selected SBML model is missing fixed parameters:\n"
        + "\n".join(
            f"  - {parameter_id}"
            for parameter_id in missing_fixed_parameters
        )
    )

# Collect parameters controlled directly by the ON/OFF or amplification values.
controlled_parameters = {
    parameter_id
    for control_type, parameter_id in (
        (INPUT_TYPE, INPUT_ID),
        (AMP_TYPE, AMP_ID),
    )
    if control_type == "parameter"
}

# Prevent fixed values from silently overriding a controlled parameter.
conflicting_fixed_parameters = sorted(
    controlled_parameters.intersection(FIXED_PARAMETERS)
)

# Report interface conflicts before any simulation begins.
if conflicting_fixed_parameters:
    # List every parameter assigned two different roles.
    raise ValueError(
        "Controlled parameters cannot also appear in fixed_parameters:\n"
        + "\n".join(
            f"  - {parameter_id}"
            for parameter_id in conflicting_fixed_parameters
        )
    )

# Require the optional secondary sweep parameter when configured.
if SWEEP_PARAM is not None and SWEEP_PARAM not in parameter_ids:
    # Report the missing sweep parameter and available alternatives.
    raise ValueError(
        f"Sweep parameter {SWEEP_PARAM!r} was not found.\n"
        f"Available global parameters: {sorted(parameter_ids)}"
    )


# ============================================================
# 5) SIMULATION HELPERS
# ============================================================

# Build the species and parameter dictionaries for one simulation condition.
def build_condition(input_value, amplification_value, params):
    """Build initial-value and parameter dictionaries for one simulation."""
    # Start with the species concentrations fixed across all conditions.
    initial_values = FIXED_INITIAL_SPECIES.copy()
    # Copy the parameter dictionary to keep the caller's values unchanged.
    parameters = params.copy()

    # Apply the requested ON or OFF value to a species when configured.
    if INPUT_TYPE == "species_initial_concentration":
        # Set the input-control species concentration.
        initial_values[INPUT_ID] = float(input_value)

    # Apply the requested ON or OFF value to a parameter when configured.
    elif INPUT_TYPE == "parameter":
        # Set the input-control parameter value.
        parameters[INPUT_ID] = float(input_value)

    # Apply amplification as an initial species concentration when configured.
    if AMP_TYPE == "species_initial_concentration":
        # Set the amplification-control species concentration.
        initial_values[AMP_ID] = float(amplification_value)

    # Apply amplification as a model parameter when configured.
    elif AMP_TYPE == "parameter":
        # Set the amplification-control parameter value.
        parameters[AMP_ID] = float(amplification_value)

    # Return independent dictionaries for the simulator.
    return initial_values, parameters


# Simulate one condition and return its observable and output rate.
def simulate_output_and_rate(input_value, amplification_value, params):
    """Simulate one condition and return the observable and d(output)/dt."""
    # Build the initial values and parameters for this condition.
    initial_values, parameters = build_condition(
        input_value=input_value,
        amplification_value=amplification_value,
        params=params,
    )

    # Run the selected SBML model over the configured time grid.
    return rr.simulate(
        initial_values=initial_values,
        parameters=parameters,
        fixed_parameters=FIXED_PARAMETERS,
        t_eval=t_model,
        observable_id=OBSERVABLE_ID,
        method="LSODA",
        rtol=RTOL,
        atol=ATOL,
        return_output_rate=True,
    )


# Compute normalized ON/OFF amplification for one parameter set.
def compute_amplification(params):
    """
    Compute normalized amplification across all configured amplification values.

                        v_ON(t,a) * v_OFF(t,a_ref)
        A(t,a) = -----------------------------------------
                        v_OFF(t,a) * v_ON(t,a_ref)
    """

    # Simulate the ON trajectory at the unamplified reference condition.
    _, v_on_ref = simulate_output_and_rate(
        INPUT_ON,
        AMP_REFERENCE,
        params,
    )

    # Simulate the OFF trajectory at the unamplified reference condition.
    _, v_off_ref = simulate_output_and_rate(
        INPUT_OFF,
        AMP_REFERENCE,
        params,
    )

    # Initialize the amplification matrix with undefined values.
    amplification = np.full(
        (len(AMP_VALUES), len(t_model)),
        np.nan,
        dtype=float,
    )

    # Evaluate every configured amplification-control value.
    for amp_index, amp_value in enumerate(AMP_VALUES):

        # Reuse the already-computed reference trajectories when possible.
        if np.isclose(amp_value, AMP_REFERENCE):
            # Reuse the reference ON rate.
            v_on = v_on_ref
            # Reuse the reference OFF rate.
            v_off = v_off_ref

        else:
            # Simulate the ON rate at this amplification-control value.
            _, v_on = simulate_output_and_rate(
                INPUT_ON,
                amp_value,
                params,
            )

            # Simulate the OFF rate at this amplification-control value.
            _, v_off = simulate_output_and_rate(
                INPUT_OFF,
                amp_value,
                params,
            )

        # Avoid artificial ratios when any relevant rate is effectively zero.
        valid = (
            (np.abs(v_on) >= RATE_FLOOR)
            & (np.abs(v_off) >= RATE_FLOOR)
            & (np.abs(v_on_ref) >= RATE_FLOOR)
            & (np.abs(v_off_ref) >= RATE_FLOOR)
        )

        # Calculate normalized amplification only where every rate is reliable.
        amplification[amp_index, valid] = (
            v_on[valid]
            * v_off_ref[valid]
        ) / (
            v_off[valid]
            * v_on_ref[valid]
        )

    # Return time courses for all amplification-control values.
    return amplification


# ============================================================
# 6) RUN OPTIONAL MODEL-PARAMETER SWEEP
# ============================================================

# Store one amplification matrix for each secondary sweep condition.
A_all = []
# Store the plot title associated with each amplification matrix.
row_labels = []

# Use the best-fit parameters directly when no secondary sweep is configured.
if SWEEP_PARAM is None:
    # Compute amplification from the unchanged best-fit parameters.
    A_all.append(compute_amplification(base_params.copy()))
    # Label the corresponding plot as the best-fit condition.
    row_labels.append("Best-fit parameters")

else:
    # Require values for a configured secondary sweep parameter.
    if SWEEP_VALUES is None:
        # Stop because the requested sweep has no values to evaluate.
        raise ValueError(
            "sweep_parameter is set but sweep_values is null."
        )

    # Compute a separate amplification matrix for each sweep value.
    for value in SWEEP_VALUES:
        # Change only the swept parameter while retaining all other best-fit values.
        params = base_params.copy()
        # Apply the current secondary parameter value.
        params[SWEEP_PARAM] = float(value)

        # Report progress for the current sweep value.
        print(f"Computing {SWEEP_PARAM} = {value:.3e}")

        # Store the amplification matrix for the current sweep value.
        A_all.append(compute_amplification(params))
        # Create the title for the current sweep panel.
        row_labels.append(f"{SWEEP_PARAM} = {value:.2e}")


# ============================================================
# 7) PLOT AMPLIFICATION TIME COURSES
# ============================================================

# Create one plot column for each parameter condition.
ncols = len(A_all)

# Create aligned amplification panels with shared axes.
fig, axes = plt.subplots(
    1,
    ncols,
    figsize=(5.2 * ncols, 4.8),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

# Normalize a single Axes object into an iterable array.
if ncols == 1:
    # Wrap the single plot axis in an array.
    axes = np.asarray([axes])

# Select the nearest simulated value for each requested plotted value.
amp_indices = [
    int(np.argmin(np.abs(AMP_VALUES - requested_value)))
    for requested_value in AMP_VALUES_TO_PLOT
]


# Format amplification-control values for the plot legend.
def format_control_value(value):
    # Include configured units when they are available.
    if AMP_UNITS:
        # Combine the numeric value with its units.
        return f"{value:g} {AMP_UNITS}"
    # Return a unitless value when no units are configured.
    return f"{value:g}"


# Populate one panel for each parameter condition.
for panel_index, (amplification, ax) in enumerate(zip(A_all, axes)):

    # Draw the unamplified reference as a dashed line at A = 1.
    ax.plot(
        t_minutes,
        np.ones_like(t_minutes),
        linestyle="--",
        linewidth=2.5,
        color="black",
        label=f"Reference: {format_control_value(AMP_REFERENCE)}",
    )

    # Plot each requested non-reference amplification value.
    for amp_index in amp_indices:
        # Get the simulated control value nearest the requested plot value.
        amp_value = AMP_VALUES[amp_index]

        # The reference is already represented by the dashed A = 1 line.
        if np.isclose(amp_value, AMP_REFERENCE):
            # Avoid duplicating the dashed reference line.
            continue

        # Draw the amplification time course for this control value.
        ax.plot(
            t_minutes,
            amplification[amp_index, :],
            linewidth=2.0,
            label=format_control_value(amp_value),
        )

    # Identify the parameter condition shown in this panel.
    ax.set_title(row_labels[panel_index])
    # Label time on the horizontal axis.
    ax.set_xlabel("Time (min)")
    # Limit the horizontal axis to the configured simulation duration.
    ax.set_xlim(0, T_END_MIN)
    # Add a light grid to make values easier to compare.
    ax.grid(True, alpha=0.3)

# Label amplification on the shared vertical axis.
axes[0].set_ylabel("Predicted amplification")

# Add the amplification-control legend to the final panel.
axes[-1].legend(
    title=AMP_LABEL,
    loc="upper right",
)

# Get the display units for the ON/OFF input.
input_units = str(INPUT_VALUES.get("units", ""))
# Get the display label for the ON/OFF input.
input_label = str(INPUT_VALUES.get("label", INPUT_ID))
# Format the ON input value for the figure title.
input_text = f"{INPUT_ON:g} {input_units}".strip()

# Summarize the selected model and ON condition above the panels.
fig.suptitle(
    f"Amplification prediction | {MODEL_NAME} | "
    f"{input_label} ON = {input_text}"
)

# Resolve the configured figure output path.
OUTPUT_FILE = resolve_path(config["output_file"])
# Create the output directory when it does not already exist.
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# Save a high-resolution copy of the amplification figure.
fig.savefig(
    OUTPUT_FILE,
    dpi=400,
    bbox_inches="tight",
)

# Report the saved figure location.
print(f"Saved figure -> {OUTPUT_FILE}")

# Display the completed figure interactively.
plt.show()
