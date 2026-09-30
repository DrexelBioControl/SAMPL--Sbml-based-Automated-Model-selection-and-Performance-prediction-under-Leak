# ============================================================
# MODEL PROPERTIES: REQUESTS, CONDITIONS AND SWEEPS
# ============================================================
# Input: the metric's requests, fitted model, conditions and run settings.
# Output: one dictionary of property values and diagnostics per request.
# No files are written here; Code 2 saves the results.
#
# 1. Reads and checks the requested properties and model IDs.
# 2. Prepares the initial conditions and parameters for each sweep value.
# 3. Calls time_course_simulation.py for concentrations and rates.
# 4. Calls steady_state_simulation.py for steady states and output plateaus.
# 5. Collects the calculated values under the names chosen by the metric.
#
# Extra property calculations can be imported through property_definitions.
# collect_properties near the end runs these steps for each named request.
# ============================================================


#------------------------------------------------------------------------------------
import importlib.util                                                               #|
from copy import deepcopy                                                           #|
from pathlib import Path                                                            #|
import sys                                                                          #|
import numpy as np                                                                  #|
                                                                                    #|
from .time_course_simulation import simulate_time_points                            #|
from .steady_state_simulation import find_steady_state, read_steady_state_settings  #|
                                                                                    #|
                                                                                    #|
BUILTIN_PROPERTIES = ("concentrations", "rates", "steady_state", "output_plateau")  #|
LIMIT_MODES = {"steady_state": "full_system", "output_plateau": "output_plateau"}   #|
#------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
# 1) READS AND CHECKS REQUESTS                                                                                      #|
                                                                                                                    #|
def read_property_requests(metric, experiment, settings, definitions=None):                                         #|
    """Read named requests from the metric, without assuming any conditions."""                                     #|
    builder = getattr(metric, "build_requests", None)                                                               #|
    if not callable(builder):                                                                                       #|
        raise ValueError("The metric must define build_requests(experiment, settings).")                            #|
    requests = builder(experiment, settings)                                                                        #|
    if not isinstance(requests, dict) or not requests:                                                              #|
        raise ValueError("build_requests must return a nonempty dictionary of named requests.")                     #|
    definitions = {} if definitions is None else definitions                                                        #|
    if not isinstance(definitions, dict):                                                                           #|
        raise ValueError("property_definitions must be a JSON object.")                                             #|
    if set(definitions).intersection(BUILTIN_PROPERTIES):                                                           #|
        raise ValueError("Imported property definitions cannot replace built-in properties.")                       #|
                                                                                                                    #|
    checked = {}                                                                                                    #|
    for name, request in requests.items():                                                                          #|
        if not isinstance(name, str) or not name.isidentifier():                                                    #|
            raise ValueError("Request names must be simple identifiers, such as treated_rates.")                    #|
        if not isinstance(request, dict):                                                                           #|
            raise ValueError("Each request must be a dictionary.")                                                  #|
        unknown = set(request).difference({"property", "initial_values", "parameters", "control_values"})           #|
        if unknown:                                                                                                 #|
            raise ValueError("Unknown request fields: {}".format(sorted(unknown)))                                  #|
        property_name = request.get("property")                                                                     #|
        if not isinstance(property_name, str) or not property_name.isidentifier():                                  #|
            raise ValueError("Each request must name a property using a simple identifier.")                        #|
        if property_name not in BUILTIN_PROPERTIES and property_name not in definitions:                            #|
            raise ValueError("Unknown property {!r}; add its property_definitions entry.".format(property_name))    #|
                                                                                                                    #|
        # A request can use the main sweep or supply its own values.                                                #|
        values = np.asarray(request.get("control_values", experiment["control_values"]), dtype=float)               #|
        if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):                                 #|
            raise ValueError("Request control_values must be a nonempty finite one-dimensional list.")              #|
        result = {"property": property_name, "control_values": values.tolist()}                                     #|
        for field in ("initial_values", "parameters"):                                                              #|
            overrides = request.get(field, {})                                                                      #|
            if not isinstance(overrides, dict):                                                                     #|
                raise ValueError("Request {} must be a dictionary.".format(field))                                  #|
            result[field] = {}                                                                                      #|
            for target, value in overrides.items():                                                                 #|
                if not isinstance(target, str) or not target:                                                       #|
                    raise ValueError("Request overrides must use species or parameter IDs.")                        #|
                value = float(value)                                                                                #|
                if not np.isfinite(value) or (field == "initial_values" and value < 0):                             #|
                    raise ValueError("Request values must be finite; initial concentrations must be nonnegative.")  #|
                result[field][target] = value                                                                       #|
        checked[name] = result                                                                                      #|
    return checked                                                                                                  #|
#--------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------
def validate_request_targets(requests, experiment, species_ids, state_species_ids,                                         #|
                             parameter_ids, sweep_parameter=None):                                                         #|
    """Check the model IDs and avoid assigning two values to the same target."""                                           #|
    if experiment["observable_id"] not in state_species_ids:                                                               #|
        raise ValueError("The observable must be a dynamic species in the model.")                                         #|
    control = experiment["analysis_control"]                                                                               #|
    for name, request in requests.items():                                                                                 #|
        initial = request.get("initial_values", {})                                                                        #|
        parameters = request.get("parameters", {})                                                                         #|
        missing_species = set(initial).difference(species_ids)                                                             #|
        missing_parameters = set(parameters).difference(parameter_ids)                                                     #|
        if missing_species or missing_parameters:                                                                          #|
            raise ValueError("Request {!r} has unknown species {} or parameters {}.".format(                               #|
                name, sorted(missing_species), sorted(missing_parameters)))                                                #|
        overrides = initial if control["type"] == "species_initial_concentration" else parameters                          #|
        if control["id"] in overrides:                                                                                     #|
            raise ValueError("Request {!r} overrides the analysis_control; use its control_values instead.".format(name))  #|
        if set(parameters).intersection(experiment["fixed_parameters"]):                                                   #|
            raise ValueError("Request {!r} overrides fixed_parameters.".format(name))                                      #|
        if sweep_parameter is not None and sweep_parameter in parameters:                                                  #|
            raise ValueError("Request {!r} overrides the secondary sweep_parameter.".format(name))                         #|
#---------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------
# 2) LOADS ANY ADDITIONAL PROPERTY CALCULATIONS                                                                               #|
                                                                                                                              #|
def load_property_functions(names, definitions, base_dir):                                                                    #|
    """Import only requested external calculations; no automatic dependencies."""                                             #|
    functions = {}                                                                                                            #|
    for name in names:                                                                                                        #|
        if name in BUILTIN_PROPERTIES:                                                                                        #|
            continue                                                                                                          #|
        definition = definitions[name]                                                                                        #|
        if not isinstance(definition, dict) or not definition.get("definition_file"):                                         #|
            raise ValueError("property_definitions.{} needs definition_file.".format(name))                                   #|
        if not isinstance(definition.get("settings", {}), dict):                                                              #|
            raise ValueError("Imported property settings must be a JSON object.")                                             #|
        path = Path(definition["definition_file"])                                                                            #|
        if not path.is_absolute():                                                                                            #|
            path = Path(base_dir) / path                                                                                      #|
        spec = importlib.util.spec_from_file_location("model_property_" + name, str(path))                                    #|
        if spec is None or spec.loader is None:                                                                               #|
            raise ImportError("Cannot import property file: {}".format(path))                                                 #|
        module = importlib.util.module_from_spec(spec)                                                                        #|
        sys.modules[spec.name] = module                                                                                       #|
        spec.loader.exec_module(module)                                                                                       #|
        function = getattr(module, "calculate_property", None)                                                                #|
        if not callable(function):                                                                                            #|
            raise ValueError("{} must define calculate_property(simulator, experiment, parameters, settings).".format(path))  #|
        functions[name] = function                                                                                            #|
    return functions                                                                                                          #|
#------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# 3) PREPARES ONE CONDITION                                                      #|
                                                                                 #|
def build_condition(experiment, control_value, parameters):                      #|
    """Copy this request's settings and apply one value of the main control."""  #|
    initial_values = experiment["initial_values"].copy()                         #|
    parameters = parameters.copy()                                               #|
    target = experiment["analysis_control"]                                      #|
    if target["type"] == "species_initial_concentration":                        #|
        initial_values[target["id"]] = float(control_value)                      #|
    elif target["type"] == "parameter":                                          #|
        parameters[target["id"]] = float(control_value)                          #|
    else:                                                                        #|
        raise ValueError("Unknown control type: {}".format(target["type"]))      #|
    return initial_values, parameters                                            #|
#---------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
# 4) COLLECTS RESULTS ACROSS A SWEEP                                                             #|
                                                                                                 #|
def collect_time_properties(simulator, experiment, parameters):                                  #|
    """Calculate concentration and rate arrays for one requested experiment."""                  #|
    shape = (len(experiment["control_values"]), len(experiment["time_model"]))                   #|
    concentrations = np.empty(shape, dtype=float)                                                #|
    rates = np.empty(shape, dtype=float)                                                         #|
    for index, value in enumerate(experiment["control_values"]):                                 #|
        # Prepares one independent condition, then passes it to the time-course solver.          #|
        initial, condition_parameters = build_condition(experiment, value, parameters)           #|
        concentrations[index], rates[index] = simulate_time_points(                              #|
            simulator, initial, condition_parameters, experiment["observable_id"],               #|
            experiment["time_model"], fixed_parameters=experiment["fixed_parameters"],           #|
            rtol=experiment["rtol"], atol=experiment["atol"],                                    #|
        )                                                                                        #|
    return {                                                                                     #|
        "concentrations": {"values": concentrations, "units": experiment["output_units"]},       #|
        "rates": {"values": rates, "units": experiment["output_units"] + " / model time unit"},  #|
    }                                                                                            #|
#-------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------
def collect_limit_property(name, simulator, experiment, parameters):                                       #|
    """Calculate verified limits for one requested experiment."""                                          #|
    options = dict(experiment.get("limit_settings", {}))                                                   #|
    options["convergence_mode"] = LIMIT_MODES[name]                                                        #|
    read_steady_state_settings(options)                                                                    #|
    records = []                                                                                           #|
    for value in experiment["control_values"]:                                                             #|
        initial, condition_parameters = build_condition(experiment, value, parameters)                     #|
        record = find_steady_state(                                                                        #|
            simulator, initial, condition_parameters, experiment["observable_id"], options,                #|
            fixed_parameters=experiment["fixed_parameters"],                                               #|
            model_time_units_per_minute=experiment["model_time_units_per_minute"],                         #|
            rtol=experiment["rtol"], atol=experiment["atol"],                                              #|
        )                                                                                                  #|
        records.append(record)                                                                             #|
        print("  {}: {} = {:g}: {} at {:g} min".format(                                                    #|
            name, experiment["analysis_control"]["id"], value, record["status"], record["time_minutes"]))  #|
    result = {}                                                                                            #|
    for field in records[0]:                                                                               #|
        if field != "convergence_mode":                                                                    #|
            result["values" if field == "value" else field] = np.asarray([r[field] for r in records])      #|
    result["convergence_mode"] = options["convergence_mode"]                                               #|
    result["state_species_ids"] = np.asarray(simulator.state_species_ids)                                  #|
    result["units"] = experiment["output_units"]                                                           #|
    print("Verified {} values: {}/{}".format(name, np.count_nonzero(result["converged"]), len(records)))   #|
    return result                                                                                          #|
#-----------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
def validate_property_result(name, result, control_count):                                                          #|
    """Check that a property returns plain arrays that can be saved directly."""                                    #|
    if not isinstance(result, dict) or not isinstance(result.get("units"), str):                                    #|
        raise ValueError("Property {} must return a flat dictionary including units.".format(name))                 #|
    values = np.asarray(result.get("values"), dtype=float)                                                          #|
    if values.ndim == 0 or values.shape[0] != control_count:                                                        #|
        raise ValueError("Property {} values must have one row per control value.".format(name))                    #|
    if np.any(np.isinf(values)):                                                                                    #|
        raise ValueError("Property {} values contain infinity; use NaN for undefined values.".format(name))         #|
    for field, value in result.items():                                                                             #|
        if not isinstance(field, str) or not field.isidentifier():                                                  #|
            raise ValueError("Property result fields must be simple identifiers.")                                  #|
        if np.asarray(value).dtype.kind not in "biufUS":                                                            #|
            raise ValueError("Property {} field {} must be a numeric/string scalar or array.".format(name, field))  #|
    return result                                                                                                   #|
#--------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------------
# 5) CARRIES OUT ALL NAMED REQUESTS                                                                                                #|
                                                                                                                                   #|
def collect_properties(requests, simulator, experiment, parameters, custom_functions=None, definitions=None):                      #|
    """Carry out each named request and return its property arrays."""                                                             #|
    properties = {}                                                                                                                #|
    time_batches = {}                                                                                                              #|
    for name, request in requests.items():                                                                                         #|
        property_name = request["property"]                                                                                        #|
        # Copies the base experiment so requests cannot change each other's settings.                                              #|
        condition = deepcopy(experiment)                                                                                           #|
        condition["initial_values"] = dict(experiment["initial_values"], **request.get("initial_values", {}))                      #|
        condition["fixed_parameters"] = dict(experiment["fixed_parameters"])                                                       #|
        condition["control_values"] = np.asarray(request.get("control_values", experiment["control_values"]), dtype=float).copy()  #|
        requested_values = condition["control_values"].copy()                                                                      #|
        condition_parameters = dict(parameters, **request.get("parameters", {}))                                                   #|
        condition["parameters"] = dict(condition_parameters, **condition["fixed_parameters"])                                      #|
                                                                                                                                   #|
        if property_name in ("concentrations", "rates"):                                                                           #|
            # Concentrations and rates share a solve only for exactly equal conditions.                                            #|
            key = (tuple(sorted(condition["initial_values"].items())),                                                             #|
                   tuple(sorted(condition_parameters.items())), tuple(condition["control_values"]))                                #|
            if key not in time_batches:                                                                                            #|
                time_batches[key] = collect_time_properties(simulator, condition, condition_parameters)                            #|
            result = {field: value.copy() if isinstance(value, np.ndarray) else value                                              #|
                      for field, value in time_batches[key][property_name].items()}                                                #|
        elif property_name in LIMIT_MODES:                                                                                         #|
            result = collect_limit_property(property_name, simulator, condition, condition_parameters)                             #|
        else:                                                                                                                      #|
            if not custom_functions or property_name not in custom_functions:                                                      #|
                raise ValueError("No calculation function was loaded for property {}.".format(property_name))                      #|
            settings = deepcopy((definitions or {}).get(property_name, {}).get("settings", {}))                                    #|
            result = custom_functions[property_name](simulator, condition, condition_parameters, settings)                         #|
                                                                                                                                   #|
        validate_property_result(name, result, len(requested_values))                                                              #|
        # Records which control values belong to this particular request's rows.                                                   #|
        result["control_values"] = requested_values                                                                                #|
        properties[name] = result                                                                                                  #|
    return properties                                                                                                              #|
#-----------------------------------------------------------------------------------------------------------------------------------
