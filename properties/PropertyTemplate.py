# ============================================================
# USER-DEFINED MODEL PROPERTY TEMPLATE
# ============================================================
# Example: the observable's initial concentration for each control.
#
# INPUT
#   simulator, experiment, parameters : supplied by properties/model_properties.py
#   settings : options from this property's property_definitions entry
#   properties/model_properties.py : build_condition prepares each condition
#
# OUTPUT
#   Initial concentrations and units returned to the chosen metric.
#   Code 2 saves the numerical results; this file writes no files.
#
# WHAT THE CODE DOES
#   1. Prepares each control value from the requested initial conditions.
#   2. Reads the initial observable amount and compartment volume.
#   3. Returns concentration = amount / volume.
#
# The run JSON registers initial_output in property_definitions
# with definition_file = properties/PropertyTemplate.py.
# The metric's build_requests then returns, for example:
#   {"initial": {"property": "initial_output"}}
# The metric reads properties["initial"]["values"].
# ============================================================


#------------------------------------------------------------------------------
import numpy as np                                                            #|
from properties.model_properties import build_condition                       #|
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------
def calculate_property(simulator, experiment, parameters, settings):                                       #|
    """Return the initial output for every requested control value."""                                     #|
    # Prepares each condition independently from its initial state.                                        #|
    values = []                                                                                            #|
    for control_value in experiment["control_values"]:                                                     #|
        initial, condition_parameters = build_condition(experiment, control_value, parameters)             #|
        model = simulator.prepare_instance(initial, condition_parameters, experiment["fixed_parameters"])  #|
        state = simulator.initial_amount_vector(model)                                                     #|
                                                                                                           #|
        # Converts the requested observable's amount into concentration.                                   #|
        observable = experiment["observable_id"]                                                           #|
        index = simulator.state_index[observable]                                                          #|
        volume = model.s[observable].compartment.size                                                      #|
        values.append(state[index] / volume)                                                               #|
                                                                                                           #|
    # Returns the quantity; the metric decides how to use it.                                              #|
    return {"values": np.asarray(values), "units": experiment["output_units"]}                             #|
#-----------------------------------------------------------------------------------------------------------
