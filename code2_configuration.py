# ============================================================
# CODE 2: SHARED MODEL DEFAULTS AND RUN CHOICES
# ============================================================
# Input: configs/model_interfaces.json and the chosen run settings.
# Output: one dictionary of model settings returned to Code 2.
# No files are written here.
#
# 1. Finds the selected model's shared defaults.
# 2. Applies this run's initial-concentration and fixed-parameter changes.
# 3. Reads which species to observe and which quantity to sweep.
# 4. Returns the combined settings, including any optional input mapping.
#
# A metric can use input_control to translate a named input into an SBML ID.
# The metric's requests decide the simulation conditions; this helper does not.
# Code 1 uses its own fitting configuration.
# ============================================================


#------------------------------------------------------------------------------
import json                                                                   #|
from pathlib import Path                                                      #|
#------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
def load_model_interface(config, model_name, base_dir):                                     #|
    """Load one model's defaults and apply explicit run-level settings."""                  #|
    # 1) Read the shared file, using the same project-relative path rule as Code 2.         #|
    if "model_interfaces" in config:                                                        #|
        raise ValueError(                                                                   #|
            "Move the inline model_interfaces block to a shared JSON file and "             #|
            "set model_interfaces_file. Put observable_id and analysis_control "            #|
            "in the run configuration, and use initial_species_overrides."                  #|
        )                                                                                   #|
    if "fixed_initial_species" in config or "input_control" in config:                      #|
        raise ValueError(                                                                   #|
            "Use initial_species_overrides in the run configuration; "                      #|
            "input_control belongs in the shared model interface."                          #|
        )                                                                                   #|
    if not config.get("model_interfaces_file"):                                             #|
        raise ValueError("The run configuration must specify model_interfaces_file.")       #|
    interface_file = Path(config["model_interfaces_file"])                                  #|
    if not interface_file.is_absolute():                                                    #|
        interface_file = Path(base_dir) / interface_file                                    #|
    with interface_file.open("r", encoding="utf-8") as handle:                              #|
        interfaces = json.load(handle)                                                      #|
    if not isinstance(interfaces, dict):                                                    #|
        raise ValueError("The shared model-interface file must contain a JSON object.")     #|
    if model_name not in interfaces:                                                        #|
        raise ValueError(                                                                   #|
            f"No shared interface exists for {model_name!r}. "                              #|
            f"Available models: {sorted(interfaces)}"                                       #|
        )                                                                                   #|
    shared = interfaces[model_name]                                                         #|
    if not isinstance(shared, dict):                                                        #|
        raise ValueError(f"The interface for {model_name!r} must be a JSON object.")        #|
    misplaced = {"observable_id", "observable_units", "analysis_control",                   #|
                 "fixed_initial_species", "fixed_parameters"}.intersection(shared)          #|
    if misplaced:                                                                           #|
        raise ValueError(                                                                   #|
            f"Unsupported fields in the shared interface: {sorted(misplaced)}. "            #|
            "Use default_initial_species, default_fixed_parameters, and "                   #|
            "default_observable_units for defaults; put analysis choices in the run JSON."  #|
        )                                                                                   #|
#--------------------------------------------------------------------------------------------


    #-----------------------------------------------------------------------------------
    # 2) Copy defaults, then replace only the keys explicitly overridden by this run.  #|
    initial_values = shared.get("default_initial_species", {})                         #|
    initial_overrides = config.get("initial_species_overrides", {})                    #|
    fixed_values = shared.get("default_fixed_parameters", {})                          #|
    fixed_overrides = config.get("fixed_parameters", {})                               #|
    for name, values in [                                                              #|
        ("default_initial_species", initial_values),                                   #|
        ("initial_species_overrides", initial_overrides),                              #|
        ("default_fixed_parameters", fixed_values),                                    #|
        ("fixed_parameters", fixed_overrides),                                         #|
    ]:                                                                                 #|
        if not isinstance(values, dict):                                               #|
            raise ValueError(f"{name} must be a JSON object.")                         #|
    initial_values = dict(initial_values)                                              #|
    initial_values.update(initial_overrides)                                           #|
    fixed_values = dict(fixed_values)                                                  #|
    fixed_values.update(fixed_overrides)                                               #|
    #-----------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------------------------
    # 3) The run chooses what to measure and vary; these are not model defaults.                     #|
    observable = config.get("observable_id")                                                         #|
    if not isinstance(observable, str) or not observable:                                            #|
        raise ValueError("The run configuration must specify a nonempty observable_id.")             #|
    control = config.get("analysis_control")                                                         #|
    if not isinstance(control, dict) or not control.get("id") or not control.get("type"):            #|
        raise ValueError("The run configuration must specify analysis_control with type and id.")    #|
    input_control = shared.get("input_control")                                                      #|
    if input_control is not None and (                                                               #|
        not isinstance(input_control, dict) or not input_control.get("id")                           #|
        or not input_control.get("type")                                                             #|
    ):                                                                                               #|
        raise ValueError("The shared input_control must have type and id, or be null.")              #|
                                                                                                     #|
    # Return the same plain interface dictionary used by the simulation code.                        #|
    return {                                                                                         #|
        "input_control": dict(input_control) if input_control is not None else None,                 #|
        "observable_id": observable,                                                                 #|
        "observable_units": config.get(                                                              #|
            "observable_units", shared.get("default_observable_units", "model concentration units")  #|
        ),                                                                                           #|
        "analysis_control": dict(control),                                                           #|
        "fixed_initial_species": initial_values,                                                     #|
        "fixed_parameters": fixed_values,                                                            #|
    }                                                                                                #|
    #-------------------------------------------------------------------------------------------------
