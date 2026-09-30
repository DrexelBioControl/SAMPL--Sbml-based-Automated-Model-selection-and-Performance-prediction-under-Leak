# ============================================================
# USER-DEFINED PERFORMANCE METRIC TEMPLATE
# ============================================================
# Example: concentration at the final requested time.
#
# INPUT
#   experiment : species, initial conditions, control grid, and times
#   settings   : this metric's options from the run JSON
#   properties : calculated concentrations from properties/model_properties.py
#
# OUTPUT
#   One concentration per control value, with a label and units.
#   Code 2 saves the figure and numerical results.
#
# WHAT THE CODE DOES
#   1. Requests concentrations using the configured conditions.
#   2. Takes the final requested value from each time course.
#   3. Returns the values for plotting against the control grid.
#
# A new metric defines its setup and formula in these two functions.
# Each request has a name, property, and optional initial_values,
# parameters, or control_values overrides.
# The final sampled value is not automatically a steady state.
# ============================================================


#------------------------------------------------------------------------------
def build_requests(experiment, settings):                                     #|
    """Request concentrations for the configured experiment."""               #|
    return {"output": {"property": "concentrations"}}                         #|
#------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
def calculate_metric(properties, experiment, settings):                        #|
    """Example: read the final concentration for each control value."""        #|
    # Reads the calculated concentration trajectories.                         #|
    output = properties["output"]                                              #|
                                                                               #|
    # Takes the last time point; copy keeps the original arrays unchanged.     #|
    metric_values = output["values"][:, -1].copy()                             #|
                                                                               #|
    # Returns one value per control; time_series would return the full array.  #|
    return {                                                                   #|
        "kind": "per_control",                                                 #|
        "values": metric_values,                                               #|
        "label": "Final output concentration",                                 #|
        "units": output["units"],                                              #|
    }                                                                          #|
#-------------------------------------------------------------------------------
