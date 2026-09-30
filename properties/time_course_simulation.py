# ============================================================
# TIME-COURSE CONCENTRATIONS AND RATES
# ============================================================
# Input: one prepared condition, the model, sampling times and tolerances.
# Output: observable concentrations and rates at those times.
# No files are written here.
#
# 1. read_time_values reads the requested sampling times from the run JSON.
# 2. simulate_time_points integrates one condition from time zero.
# 3. It returns concentrations and their ODE rates from the same trajectory.
#
# properties/model_properties.py prepares each condition and collects the sweep.
# Times passed to the solver use model time units; rates use those units too.
# ============================================================


#------------------------------------------------------------------------------
import numpy as np                                                            #|
#------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------
def read_time_values(config, required=False):                                                                            #|
    """Explicit times take precedence over a regular time_end_min/time_points grid."""                                   #|
    if "time_values_min" in config:                                                                                      #|
        times = np.asarray(config["time_values_min"], dtype=float)                                                       #|
    elif "time_end_min" in config or "time_points" in config:                                                            #|
        end = float(config.get("time_end_min", 0))                                                                       #|
        count = config.get("time_points", 0)                                                                             #|
        if not np.isfinite(end) or end <= 0 or isinstance(count, bool) or int(count) != count or count < 2:              #|
            raise ValueError("time_end_min must be positive and time_points must be an integer >= 2.")                   #|
        times = np.linspace(0.0, end, int(count))                                                                        #|
    elif required:                                                                                                       #|
        raise ValueError("Concentrations/rates require time_values_min or time_end_min and time_points.")                #|
    else:                                                                                                                #|
        return np.array([], dtype=float)                                                                                 #|
    if (times.ndim != 1 or times.size == 0 or not np.all(np.isfinite(times))                                             #|
            or np.any(times < 0) or np.any(np.diff(times) <= 0)):                                                        #|
        raise ValueError("Sampling times must be a nonempty, finite, strictly increasing list of nonnegative minutes.")  #|
    return times                                                                                                         #|
#-------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
def simulate_time_points(simulator, initial_values, parameters, observable_id, times,  #|
                         fixed_parameters=None, rtol=1e-6, atol=1e-6):                 #|
    """Return concentration and rate, always integrating from time zero."""            #|
    # Copies the supplied condition so this calculation keeps its inputs intact.       #|
    initial = initial_values.copy()                                                    #|
    parameters = parameters.copy()                                                     #|
    fixed_parameters = dict(fixed_parameters or {})                                    #|
    if times.size == 0:                                                                #|
        raise ValueError("A time property needs at least one sampling time.")          #|
    # A time-zero request reads the initial state and needs no integration.            #|
    if times.size == 1 and times[0] == 0:                                              #|
        model = simulator.prepare_instance(initial, parameters, fixed_parameters)      #|
        state = simulator.initial_amount_vector(model)                                 #|
        index = simulator.state_index[observable_id]                                   #|
        volume = model.s[observable_id].compartment.size                               #|
        return (np.asarray([state[index] / volume]),                                   #|
                np.asarray([simulator.rhs(model, 0.0, state)[index] / volume]))        #|
    # Starts at zero even when the requested measurements begin later.                 #|
    prepend_zero = times[0] > 0                                                        #|
    integration_times = np.r_[0.0, times] if prepend_zero else times                   #|
    output, rate = simulator.simulate(                                                 #|
        initial_values=initial, parameters=parameters,                                 #|
        fixed_parameters=fixed_parameters,                                             #|
        t_eval=integration_times, observable_id=observable_id,                         #|
        method="LSODA", rtol=rtol, atol=atol,                                          #|
        return_output_rate=True,                                                       #|
    )                                                                                  #|
    return (output[1:], rate[1:]) if prepend_zero else (output, rate)                  #|
#---------------------------------------------------------------------------------------
