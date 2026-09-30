# ============================================================
# STEADY-STATE OR OUTPUT-PLATEAU SENSITIVITY
# ============================================================
# Calculates S(theta) = (dx/dtheta) * (theta/x).
#
# INPUT
#   metric.settings : response_property, concentration_floor, absolute
#   experiment      : parameter grid and configured conditions
#   properties      : verified limiting outputs and convergence flags
#
# OUTPUT
#   Sensitivity values, dx/dtheta, labels, and units.
#   Code 2 saves the metric values and property diagnostics.
#   dx/dtheta is returned for inspection, but is not saved separately.
#
# WHAT THE CODE DOES
#   1. Requests steady_state or output_plateau for the parameter sweep.
#   2. Differentiates the limiting outputs with respect to theta.
#   3. Multiplies each derivative by that point's theta/x.
#
# A failed convergence check invalidates derivative stencils that use it.
# Full-system and output-only convergence remain distinct properties.
# Integration and convergence checks happen in the property helpers.
# ============================================================


#------------------------------------------------------------------------------
import numpy as np                                                            #|
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
def get_response_property(settings):                                                     #|
    """Tell Code 2 which limiting output this calculation needs."""                      #|
                                                                                         #|
    response_property = settings.get("response_property", "steady_state")                #|
    if response_property not in ("steady_state", "output_plateau"):                      #|
        raise ValueError(                                                                #|
            "Sensitivity response_property must be 'steady_state' or 'output_plateau'."  #|
        )                                                                                #|
    return response_property                                                             #|
#-----------------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def build_requests(experiment, settings):                                     #|
    """Describe the limiting output needed by sensitivity."""                 #|
    return {"response": {"property": get_response_property(settings)}}        #|
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
def calculate_metric(properties, experiment, settings):                                   #|
    """Return one normalized sensitivity per supplied parameter value."""                 #|
                                                                                          #|
    # 1) Read the settings specific to this metric.                                       #|
    concentration_floor = float(settings.get("concentration_floor", 1e-8))                #|
    if not np.isfinite(concentration_floor) or concentration_floor <= 0:                  #|
        raise ValueError("Sensitivity concentration_floor must be finite and positive.")  #|
                                                                                          #|
    absolute = settings.get("absolute", False)                                            #|
    if not isinstance(absolute, bool):                                                    #|
        raise ValueError("Sensitivity absolute must be true or false.")                   #|
#------------------------------------------------------------------------------------------


    #--------------------------------------------------------------------------------------------------
    # 2) Read the experimental parameter grid and the requested model property.                       #|
    # theta is supplied by the experiment; only the response requires solving.                        #|
    response_property = get_response_property(settings)                                               #|
    response = properties["response"]                                                                 #|
    theta = np.asarray(experiment["control_values"], dtype=float)                                     #|
    xss = np.asarray(response["values"], dtype=float)                                                 #|
    converged = np.asarray(response["converged"])                                                     #|
                                                                                                      #|
    mode = response.get("convergence_mode")                                                           #|
    expected_mode = "full_system" if response_property == "steady_state" else "output_plateau"        #|
    if mode != expected_mode:                                                                         #|
        raise ValueError(                                                                             #|
            f"Sensitivity {response_property} requires convergence_mode = '{expected_mode}'."         #|
        )                                                                                             #|
    control = experiment.get("analysis_control", {})                                                  #|
    if control.get("type") != "parameter":                                                            #|
        raise ValueError("Sensitivity requires analysis_control.type = 'parameter'.")                 #|
    if theta.ndim != 1 or theta.size < 3:                                                             #|
        raise ValueError("Sensitivity requires at least three parameter values in a 1D grid.")        #|
    if not np.all(np.isfinite(theta)) or np.any(theta <= 0):                                          #|
        raise ValueError("Sensitivity parameter values must be finite and positive.")                 #|
    if np.any(np.diff(theta) <= 0):                                                                   #|
        raise ValueError("Sensitivity parameter values must be strictly increasing.")                 #|
    if xss.shape != theta.shape or converged.shape != theta.shape:                                    #|
        raise ValueError("Steady-state values and convergence flags must match the parameter grid.")  #|
    if converged.dtype.kind != "b":                                                                   #|
        raise ValueError("Steady-state convergence flags must be boolean values.")                    #|
    #--------------------------------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # 3) Differentiate, then normalize at each parameter value.               #|
    # concentration_floor uses the same units as the selected species.        #|
    valid = converged & np.isfinite(xss) & (xss >= concentration_floor)       #|
    checked_xss = np.where(valid, xss, np.nan)                                #|
                                                                              #|
    # Every interior stencil uses i-1, i, i+1. Each endpoint uses the         #|
    # nearest three points. A failed point invalidates its whole stencil.     #|
    valid_stencil = valid.copy()                                              #|
    valid_stencil[1:-1] = valid[:-2] & valid[1:-1] & valid[2:]                #|
    valid_stencil[0] = np.all(valid[:3])                                      #|
    valid_stencil[-1] = np.all(valid[-3:])                                    #|
                                                                              #|
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):       #|
        dxss_dtheta = np.gradient(checked_xss, theta, edge_order=2)           #|
        metric_values = dxss_dtheta * theta / checked_xss                     #|
                                                                              #|
    dxss_dtheta[~valid_stencil | ~np.isfinite(dxss_dtheta)] = np.nan          #|
    metric_values[~valid_stencil | ~np.isfinite(metric_values)] = np.nan      #|
    if absolute:                                                              #|
        metric_values = np.abs(metric_values)                                 #|
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # 4) Return the values and the information needed for a general plot.     #|
    quantity = "steady-state" if mode == "full_system" else "output-plateau"  #|
    label = f"Normalized {quantity} sensitivity"                              #|
    if absolute:                                                              #|
        label = f"Absolute normalized {quantity} sensitivity"                 #|
    return {                                                                  #|
        "kind": "per_control",                                                #|
        "values": metric_values,                                              #|
        "label": label,                                                       #|
        "units": "",                                                          #|
        "dxss_dtheta": dxss_dtheta,                                           #|
    }                                                                         #|
    #--------------------------------------------------------------------------
