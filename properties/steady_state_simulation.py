# ============================================================
# STEADY STATES AND OUTPUT PLATEAUS
# ============================================================
# Input: one prepared condition, the model, observable and convergence settings.
# Output: a verified limiting value (or NaN), final state and diagnostics.
# No files are written here; Code 2 saves the results.
#
# 1. Integrates the full model over successive time windows.
# 2. Checks concentration changes and rates against the chosen tolerances.
# 3. Requires several consecutive passing windows.
# 4. For an output plateau, extends the horizon to check that it persists.
#
# full_system: every dynamic species must settle.
# output_plateau: only the selected output must settle; all species are simulated.
# max_time_min includes the extension needed to confirm a plateau.
#
# properties/model_properties.py calls find_steady_state for each sweep value.
# Settings come from limit_settings in the run JSON.
# See README.md for the convergence rules and numerical interpretation.
# ============================================================


#------------------------------------------------------------------------------
import numpy as np                                                            #|
from scipy.integrate import solve_ivp                                         #|
#------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
def read_steady_state_settings(settings):                                                                #|
    """Read a small set of simulation settings with explicit units."""                                   #|
    mode = settings.get("convergence_mode", "full_system")                                               #|
    if mode not in ("full_system", "output_plateau"):                                                    #|
        raise ValueError("steady_state.convergence_mode must be full_system or output_plateau.")         #|
    extension_factor = float(settings.get("plateau_extension_factor", 2.0))                              #|
    if not np.isfinite(extension_factor) or extension_factor <= 1:                                       #|
        raise ValueError("steady_state.plateau_extension_factor must be finite and > 1.")                #|
                                                                                                         #|
    values = {                                                                                           #|
        "max_time_min": float(settings.get("max_time_min", 2000.0)),                                     #|
        "check_interval_min": float(settings.get("check_interval_min", 100.0)),                          #|
        "rate_tolerance": float(settings.get("rate_tolerance", 1e-8)),                                   #|
        "change_atol": float(settings.get("change_atol", 1e-6)),                                         #|
        "change_rtol": float(settings.get("change_rtol", 1e-6)),                                         #|
    }                                                                                                    #|
    for name, value in values.items():                                                                   #|
        if not np.isfinite(value) or value < 0:                                                          #|
            raise ValueError(f"steady_state.{name} must be finite and nonnegative.")                     #|
        if name != "change_rtol" and value == 0:                                                         #|
            raise ValueError(f"steady_state.{name} must be positive.")                                   #|
                                                                                                         #|
    for name, default, minimum in [("consecutive_checks", 3, 2), ("samples_per_check", 9, 3)]:           #|
        value = settings.get(name, default)                                                              #|
        if isinstance(value, bool) or not np.isfinite(value) or int(value) != value or value < minimum:  #|
            raise ValueError(f"steady_state.{name} must be an integer >= {minimum}.")                    #|
        values[name] = int(value)                                                                        #|
                                                                                                         #|
    minimum_time = values["check_interval_min"] * values["consecutive_checks"]                           #|
    if values["max_time_min"] < minimum_time:                                                            #|
        raise ValueError("max_time_min must allow all consecutive full check intervals.")                #|
    values["convergence_mode"] = mode                                                                    #|
    values["plateau_extension_factor"] = extension_factor                                                #|
    return values                                                                                        #|
#---------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
def find_steady_state(                                                                      #|
    simulator, initial_values, parameters, observable_id, settings,                         #|
    fixed_parameters=None, model_time_units_per_minute=60.0,                                #|
    rtol=1e-8, atol=1e-10,                                                                  #|
):                                                                                          #|
    """Return a steady-state or output-plateau estimate and its diagnostics.

    rate_tolerance: concentration per MODEL time unit.
    change_atol: concentration; change_rtol: dimensionless.
    For this function, solver atol is in concentration units and is
    converted to species amounts using each compartment volume.
    """
    options = read_steady_state_settings(settings)                                          #|
    mode = options["convergence_mode"]                                                      #|
    time_factor = float(model_time_units_per_minute)                                        #|
    if not np.isfinite(time_factor) or time_factor <= 0:                                    #|
        raise ValueError("model_time_units_per_minute must be finite and positive.")        #|
    if not np.isfinite(rtol) or rtol <= 0 or not np.isfinite(atol) or atol <= 0:            #|
        raise ValueError("Steady-state solver rtol and atol must be finite and positive.")  #|
#--------------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------
    # 1) Prepare one independent instance at the requested parameter value.                      #|
    model = simulator.prepare_instance(initial_values, parameters, fixed_parameters)             #|
    species_ids = simulator.state_species_ids                                                    #|
    observable_index = simulator.state_index[observable_id]                                      #|
    volumes = np.asarray([model.s[name].compartment.size for name in species_ids], dtype=float)  #|
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):                                 #|
        raise ValueError("Steady-state simulation requires finite positive compartment sizes.")  #|
                                                                                                 #|
    state = simulator.initial_amount_vector(model).copy()                                        #|
    time = 0.0                                                                                   #|
    max_time = options["max_time_min"] * time_factor                                             #|
    interval = options["check_interval_min"] * time_factor                                       #|
    if not np.isfinite(max_time) or not np.isfinite(interval):                                   #|
        raise ValueError("Steady-state times overflow after conversion to model time.")          #|
    passing_checks = 0                                                                           #|
    full_system_checks = 0                                                                       #|
    plateau_start_time = np.nan                                                                  #|
    plateau_target_time = np.nan                                                                 #|
    plateau_min = np.nan                                                                         #|
    plateau_max = np.nan                                                                         #|
    plateau_change = np.nan                                                                      #|
    output_rate = np.nan                                                                         #|
    output_change = np.nan                                                                       #|
    final_concentrations = state / volumes                                                       #|
    final_output = float(final_concentrations[observable_index])                                 #|
    max_rate = np.nan                                                                            #|
    max_change = np.nan                                                                          #|
    worst_species = ""                                                                           #|
    status = "max_time_reached"                                                                  #|
    converged = False                                                                            #|
    last_window_valid = False                                                                    #|
                                                                                                 #|
    if not np.all(np.isfinite(state)) or np.any(final_concentrations < 0):                       #|
        raise ValueError("Initial concentrations must be finite and nonnegative.")               #|
    #---------------------------------------------------------------------------------------------


    #----------------------------------------------------------------------------------
    # 2) Continue integration in fixed-duration windows.                              #|
    while time < max_time:                                                            #|
        last_window_valid = False                                                     #|
        end_time = min(time + interval, max_time)                                     #|
        sample_times = np.linspace(time, end_time, options["samples_per_check"])      #|
        try:                                                                          #|
            solution = solve_ivp(                                                     #|
                lambda t, y: simulator.rhs(model, t, y),                              #|
                (time, end_time),                                                     #|
                state,                                                                #|
                t_eval=sample_times,                                                  #|
                method="LSODA",                                                       #|
                rtol=rtol,                                                            #|
                atol=atol * volumes,                                                  #|
                # Sample the extended trajectory without taking oversized steps.      #|
                max_step=interval / (options["samples_per_check"] - 1)                #|
                if mode == "output_plateau" else np.inf,                              #|
            )                                                                         #|
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:  #|
            status = "integration_failed: " + str(exc)                                #|
            break                                                                     #|
                                                                                      #|
        if not solution.success:                                                      #|
            status = "integration_failed: " + solution.message                        #|
            break                                                                     #|
                                                                                      #|
        concentrations = solution.y / volumes[:, None]                                #|
        time = float(solution.t[-1])                                                  #|
        state = solution.y[:, -1].copy()                                              #|
        final_concentrations = concentrations[:, -1].copy()                           #|
        final_output = float(final_concentrations[observable_index])                  #|
                                                                                      #|
        if not np.all(np.isfinite(concentrations)):                                   #|
            status = "nonfinite_concentration"                                        #|
            break                                                                     #|
        if np.any(concentrations < -10.0 * atol):                                     #|
            status = "negative_concentration"                                         #|
            break                                                                     #|
    #----------------------------------------------------------------------------------


        #--------------------------------------------------------------------------
        # 3) Keep whole-system diagnostics in BOTH modes.                         #|
        rates = np.column_stack([                                                 #|
            simulator.rhs(model, float(t), solution.y[:, index]) / volumes        #|
            for index, t in enumerate(solution.t)                                 #|
        ])                                                                        #|
        if not np.all(np.isfinite(rates)):                                        #|
            status = "nonfinite_rate"                                             #|
            break                                                                 #|
                                                                                  #|
        last_window_valid = True                                                  #|
                                                                                  #|
        rates_by_species = np.max(np.abs(rates), axis=1)                          #|
        worst_index = int(np.argmax(rates_by_species))                            #|
        worst_species = species_ids[worst_index]                                  #|
        max_rate = float(rates_by_species[worst_index])                           #|
                                                                                  #|
        scale = options["change_atol"] + options["change_rtol"] * np.max(         #|
            np.abs(concentrations), axis=1,                                       #|
        )                                                                         #|
        changes_by_species = np.ptp(concentrations, axis=1) / scale               #|
        max_change = float(np.max(changes_by_species))                            #|
        output_rate = float(rates_by_species[observable_index])                   #|
        output_change = float(changes_by_species[observable_index])               #|
                                                                                  #|
        # The absolute rate check prevents a growing species from passing         #|
        # merely because its fractional change becomes small at large times.      #|
        full_window = end_time - sample_times[0] >= interval * (1.0 - 1e-12)      #|
        full_system_passed = (                                                    #|
            full_window                                                           #|
            and max_rate <= options["rate_tolerance"]                             #|
            and max_change <= 1.0                                                 #|
        )                                                                         #|
        full_system_checks = full_system_checks + 1 if full_system_passed else 0  #|
        #--------------------------------------------------------------------------


        #----------------------------------------------------------------------
        # 4) Full-system mode retains the original stopping rule.             #|
        if mode == "full_system":                                             #|
            if full_system_checks >= options["consecutive_checks"]:           #|
                converged = True                                              #|
                status = "converged"                                          #|
                break                                                         #|
            continue                                                          #|
        #----------------------------------------------------------------------


        #-----------------------------------------------------------------------------------
        # 5) Output mode first identifies an apparent plateau.                             #|
        output_passed = output_rate <= options["rate_tolerance"] and output_change <= 1.0  #|
        output_samples = concentrations[observable_index]                                  #|
        if np.isnan(plateau_start_time):                                                   #|
            passing_checks = passing_checks + 1 if full_window and output_passed else 0    #|
            if passing_checks >= options["consecutive_checks"]:                            #|
                plateau_start_time = time                                                  #|
                # Require a longer horizon AND several more complete windows.              #|
                target = max(                                                              #|
                    time * options["plateau_extension_factor"],                            #|
                    time + options["consecutive_checks"] * interval,                       #|
                )                                                                          #|
                # Round up to a full check window; never shorten to fit the time cap.      #|
                plateau_target_time = float(np.ceil(target / interval) * interval)         #|
                plateau_min = float(np.min(output_samples))                                #|
                plateau_max = float(np.max(output_samples))                                #|
                plateau_change = output_change                                             #|
            continue                                                                       #|
        #-----------------------------------------------------------------------------------


        #------------------------------------------------------------------------------------
        # 6) Check persistence throughout the extension, including excursions               #|
        # that return to the same endpoint and slow cumulative drift.                       #|
        plateau_min = min(plateau_min, float(np.min(output_samples)))                       #|
        plateau_max = max(plateau_max, float(np.max(output_samples)))                       #|
        plateau_scale = options["change_atol"] + options["change_rtol"] * max(              #|
            abs(plateau_min), abs(plateau_max),                                             #|
        )                                                                                   #|
        plateau_change = (plateau_max - plateau_min) / plateau_scale                        #|
        if not output_passed or plateau_change > 1.0:                                       #|
            # This candidate did not persist. A later plateau may still qualify.            #|
            passing_checks = 0                                                              #|
            plateau_start_time = np.nan                                                     #|
            plateau_target_time = np.nan                                                    #|
            plateau_min = np.nan                                                            #|
            plateau_max = np.nan                                                            #|
            continue                                                                        #|
                                                                                            #|
        if full_window and time >= plateau_target_time * (1.0 - 1e-12):                     #|
            converged = True                                                                #|
            status = "output_plateau_verified"                                              #|
            break                                                                           #|
                                                                                            #|
    if not converged and status == "max_time_reached" and np.isfinite(plateau_start_time):  #|
        status = "plateau_extension_incomplete"                                             #|
        #------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------------------
    # 7) Failed runs retain their endpoint, but never return it as a verified limit.                         #|
    return {                                                                                                 #|
        "value": final_output if converged else np.nan,                                                      #|
        "convergence_mode": mode,                                                                            #|
        "full_system_converged": last_window_valid and full_system_checks >= options["consecutive_checks"],  #|
        "output_max_abs_rate": output_rate,                                                                  #|
        "output_max_scaled_change": output_change,                                                           #|
        "plateau_start_time_minutes": plateau_start_time / time_factor,                                      #|
        "plateau_target_time_minutes": plateau_target_time / time_factor,                                    #|
        "plateau_max_scaled_change": plateau_change,                                                         #|
        "converged": converged,                                                                              #|
        "status": status,                                                                                    #|
        "time_minutes": time / time_factor,                                                                  #|
        "final_output": final_output,                                                                        #|
        "max_abs_rate": max_rate,                                                                            #|
        "max_scaled_change": max_change,                                                                     #|
        "worst_rate_species": worst_species,                                                                 #|
        "final_concentrations": final_concentrations,                                                        #|
    }                                                                                                        #|
    #---------------------------------------------------------------------------------------------------------
