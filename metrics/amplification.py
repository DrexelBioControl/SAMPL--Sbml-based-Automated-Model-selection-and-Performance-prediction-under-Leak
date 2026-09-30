# ============================================================
# AMPLIFICATION PERFORMANCE METRIC
# ============================================================
# Compares ON/OFF output-rate discrimination with a reference.
#
# INPUT
#   metric.settings : input_values, reference_value, and rate_floor
#   experiment      : model input mapping and analysis-control grid
#   properties      : requested rates returned by properties/model_properties.py
#
# OUTPUT
#   Amplification curves, reference curve, and plot labels.
#   Code 2 saves the figure and numerical results.
#
# WHAT THE CODE DOES
#   1. Requests ON/OFF rates for the sweep and reference.
#   2. Calculates A = (v_ON * v_OFF_ref) / (v_OFF * v_ON_ref).
#   3. Returns NaN where any rate is below the chosen rate_floor.
#
# The input values and reference belong to this metric's setup.
# Rates and rate_floor use concentration per MODEL time unit.
# ============================================================


#------------------------------------------------------------------------------
import numpy as np                                                            #|
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------
# Define function to build the requests for the amplification metric                                  #|
def build_requests(experiment, settings):                                                             #|
    """Describe the four rate calculations used by amplification."""                                  #|
                                                                                                      #|
    # Reads the input values and the analysis-control reference.                                      #|
    inputs = settings.get("input_values", {})                                                         #|
    if not isinstance(inputs, dict) or not {"on_value", "off_value"}.issubset(inputs):                #|
        raise ValueError("Amplification metric.settings.input_values needs on_value and off_value.")  #|
    if "reference_value" not in settings:                                                             #|
        raise ValueError("Amplification metric.settings needs reference_value.")                      #|
    on_value = float(inputs["on_value"])                                                              #|
    off_value = float(inputs["off_value"])                                                            #|
    reference = float(settings["reference_value"])                                                    #|
    if not np.all(np.isfinite([on_value, off_value, reference])):                                     #|
        raise ValueError("Amplification input values and reference_value must be finite.")            #|
                                                                                                      #|
    # Finds whether the model's input changes a species or a parameter.                               #|
    target = experiment.get("input_control")                                                          #|
    on_changes, off_changes = {}, {}                                                                  #|
    if target is not None:                                                                            #|
        fields = {"species_initial_concentration": "initial_values", "parameter": "parameters"}       #|
        if target.get("type") not in fields or not target.get("id"):                                  #|
            raise ValueError("Amplification input_control needs a supported type and an SBML id.")    #|
        field = fields[target["type"]]                                                                #|
        on_changes[field] = {target["id"]: on_value}                                                  #|
        off_changes[field] = {target["id"]: off_value}                                                #|
                                                                                                      #|
    # The first two requests use the full sweep; the references use one value.                        #|
    return {                                                                                          #|
        "on_rates": dict(on_changes, property="rates"),                                               #|
        "off_rates": dict(off_changes, property="rates"),                                             #|
        "reference_on": dict(on_changes, property="rates", control_values=[reference]),               #|
        "reference_off": dict(off_changes, property="rates", control_values=[reference]),             #|
    }                                                                                                 #|
#------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# Define function to calculate the amplification curve for one time course            #|
def calculate_one_curve(on_rate, off_rate, on_reference, off_reference, rate_floor):  #|
    """Apply the original amplification definition to one time course."""             #|
                                                                                      #|
    # Start with undefined values, including times before any output is made.         #|
    values = np.full(on_rate.shape, np.nan, dtype=float)                              #|
                                                                                      #|
    # Preserve the original four-rate rule, including the >= boundary.                #|
    valid = (                                                                         #|
        (np.abs(on_rate) >= rate_floor)                                               #|
        & (np.abs(off_rate) >= rate_floor)                                            #|
        & (np.abs(on_reference) >= rate_floor)                                        #|
        & (np.abs(off_reference) >= rate_floor)                                       #|
    )                                                                                 #|
                                                                                      #|
    # Keep the original multiplication and division order.                            #|
    values[valid] = (on_rate[valid] * off_reference[valid]) / (                       #|
        off_rate[valid] * on_reference[valid]                                         #|
    )                                                                                 #|
                                                                                      #|
    return values                                                                     #|
#--------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# Define function to calculate amplification from the simulated ON/OFF and reference output rates  #|
def calculate_metric(properties, experiment, settings):                                            #|
    """Return amplification curves and their reference in a plain dictionary."""                   #|
    # 1) Reads the minimum permitted magnitude of each rate.                                       #|
    rate_floor = float(settings.get("rate_floor", 1e-8))                                           #|
    if not np.isfinite(rate_floor) or rate_floor <= 0:                                             #|
        raise ValueError("Amplification rate_floor must be a finite positive number.")             #|
#---------------------------------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # 2) Reads the requested arrays; a reference request contains one row.    #|
    on_rate = properties["on_rates"]["values"]                                #|
    off_rate = properties["off_rates"]["values"]                              #|
    on_reference = properties["reference_on"]["values"][0]                    #|
    off_reference = properties["reference_off"]["values"][0]                  #|
    #--------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    # 3) Calculates one curve per analysis-control value.                      #|
    metric_values = np.full(on_rate.shape, np.nan, dtype=float)                #|
    for control_index in range(len(experiment["control_values"])):             #|
        metric_values[control_index] = calculate_one_curve(                    #|
            on_rate[control_index], off_rate[control_index],                   #|
            on_reference, off_reference, rate_floor,                           #|
        )                                                                      #|
                                                                               #|
    # The reference equals one only where the same four-rate checks pass.      #|
    reference_values = calculate_one_curve(                                    #|
        on_reference, off_reference, on_reference, off_reference, rate_floor,  #|
    )                                                                          #|
    #---------------------------------------------------------------------------


    #--------------------------------------------------------------------------------------------------
    # 4) Supplies condition and reference labels for Code 2's general plot.                           #|
    reference = float(settings["reference_value"])                                                    #|
    reference_units = experiment.get("control_units", "")                                             #|
    reference_label = f"Reference: {reference:g} {reference_units}".rstrip()                          #|
    inputs = settings["input_values"]                                                                 #|
    if experiment.get("input_control") is None:                                                       #|
        condition_label = "Configured initial conditions"                                             #|
    else:                                                                                             #|
        input_label = str(inputs.get("label", experiment["input_control"]["id"]))                     #|
        input_units = str(inputs.get("units", ""))                                                    #|
        condition_label = f"{input_label} ON = {float(inputs['on_value']):g} {input_units}".rstrip()  #|
        condition_label += f"; OFF = {float(inputs['off_value']):g} {input_units}".rstrip()           #|
                                                                                                      #|
    return {                                                                                          #|
        "kind": "time_series",                                                                        #|
        "values": metric_values,                                                                      #|
        "label": "Predicted amplification",                                                           #|
        "units": "",                                                                                  #|
        "reference_values": reference_values,                                                         #|
        "reference_control_value": reference,                                                         #|
        "reference_label": reference_label,                                                           #|
        "condition_label": condition_label,                                                           #|
    }                                                                                                 #|
    #--------------------------------------------------------------------------------------------------
