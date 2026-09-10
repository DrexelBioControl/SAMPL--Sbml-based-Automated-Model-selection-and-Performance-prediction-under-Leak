# ============================================================
# SBMLtoODEpy -> SciPy solve_ivp ADAPTER
# ============================================================
#
# SPIRIT OF THIS FILE
# -------------------
# SBMLtoODEpy reads an SBML file and writes a normal Python class containing
# the reactions and the resulting ODE right-hand side. Its generated private
# method has the form:
#
#     model._SolveReactions(y, t)
#
# because SBMLtoODEpy originally calls scipy.integrate.odeint.
# Our established workflow instead uses solve_ivp, which expects:
#
#     rhs(t, y)
#
# This adapter bridges those two conventions without rewriting the biochemical
# equations by hand. SBML remains the source of truth for the mechanism.
# ============================================================

# Postpone evaluation of type hints so this file remains friendly to Python 3.7.
from __future__ import annotations

# Load Python modules dynamically from their generated file paths.
import importlib.util
# Read the source code of generated functions when determining species order.
import inspect
# Search generated source text with regular expressions.
import re
# Register dynamically imported modules with the running Python interpreter.
import sys
# Represent file and folder paths in an operating-system-independent way.
from pathlib import Path

# NumPy supplies arrays and numerical operations for state vectors.
import numpy as np
# solve_ivp numerically integrates the model's differential equations.
from scipy.integrate import solve_ivp

# Try to import the package that converts SBML into Python equations.
try:
    import sbmltoodepy
# Replace a vague import failure with instructions specific to this workflow.
except ImportError as exc:
    raise RuntimeError(
        "SBMLtoODEpy is not installed in this Python environment.\n"
        "This package is old and officially targets Python 3.5-3.7.\n"
        "Create the supplied Python 3.7 conda environment before running."
    ) from exc


class SBMLtoODEpySolveIVPModel:
    # This class generates, imports, and simulates one SBML model with solve_ivp.
    #
    # It supports the reaction-network features used by the current ctRSD
    # models: reaction-controlled species, global parameters, local reaction
    # parameters, compartments, and function definitions.
    #
    # Models containing rate rules, events, or algebraic rules need separate
    # validation before this adapter can safely use them.

    def __init__(self, sbml_file: Path, generated_models_dir: Path):
        # Convert the supplied SBML filename into a complete absolute path.
        self.sbml_file = Path(sbml_file).resolve()
        # Use the SBML filename without .xml as the model's short name.
        self.model_name = self.sbml_file.stem
        # Convert the generated-model folder into a complete absolute path.
        self.generated_models_dir = Path(generated_models_dir).resolve()
        # Create the generated-model folder if it does not already exist.
        self.generated_models_dir.mkdir(parents=True, exist_ok=True)

        # Replace characters that Python class and module names cannot contain.
        safe_stem = re.sub(r"\W|^(?=\d)", "_", self.model_name)
        # Build the path where SBMLtoODEpy will write the generated Python model.
        self.generated_file = self.generated_models_dir / f"{safe_stem}_generated.py"
        # Build a unique and valid name for the generated Python class.
        self.class_name = f"Generated_{safe_stem}"

        # Generate the Python model only when its cached version is outdated.
        self._generate_when_needed()
        # Correct SBMLtoODEpy's extra leading matrix delimiter before importing.
        self._repair_generated_stoichiometric_matrix()
        # Import the generated class and keep it for later model construction.
        self._model_class = self._import_generated_class()

        # Make a temporary model instance so its contents can be inspected.
        probe = self._model_class()
        # Save the ordered identifiers of all species found in the model.
        self.species_ids = list(probe.s.keys())
        # Save the ordered identifiers of all global model parameters.
        self.parameter_ids = list(probe.p.keys())
        # Save the ordered identifiers of all model reactions.
        self.reaction_ids = list(probe.r.keys())
        # Recover the exact species order used in the generated ODE state vector.
        self.state_species_ids = self._find_state_species_order(probe)

        # Confirm that every SBML species appears exactly once in the state vector.
        if set(self.state_species_ids) != set(self.species_ids):
            # List species present in SBML but missing from the generated state.
            missing = sorted(set(self.species_ids) - set(self.state_species_ids))
            # List unexpected state entries that are not ordinary SBML species.
            extra = sorted(set(self.state_species_ids) - set(self.species_ids))
            raise NotImplementedError(
                f"{self.sbml_file.name}: generated state vector is not a simple "
                "species-only vector. This usually means the SBML contains rate "
                "rules or another feature not yet supported by this adapter. "
                f"Missing species={missing}; extra entries={extra}."
            )

        # Create a quick species-ID-to-row-number lookup for solution arrays.
        self.state_index = {
            # Store the numerical position associated with each species ID.
            species_id: index
            # Visit every species ID together with its zero-based position.
            for index, species_id in enumerate(self.state_species_ids)
        }

    def _generate_when_needed(self) -> None:
        # Regenerate the Python equations if no generated file exists or if the
        # SBML file was edited more recently than the generated Python file.
        # The cache needs rebuilding if it is absent or older than the SBML.
        needs_generation = (
            not self.generated_file.exists()
            or self.generated_file.stat().st_mtime < self.sbml_file.stat().st_mtime
        )

        # Stop here when the existing generated file is already current.
        if not needs_generation:
            return

        # Ask SBMLtoODEpy to translate the SBML into an importable Python class.
        sbmltoodepy.ParseAndCreateModel(
            str(self.sbml_file),
            outputFilePath=str(self.generated_file),
            className=self.class_name,
        )


    def _repair_generated_stoichiometric_matrix(self) -> None:
        # ============================================================
        # WHY THIS REPAIR IS NEEDED
        # ============================================================
        # Some SBMLtoODEpy versions begin every generated matrix row with
        # an empty entry. For example, a seven-reaction row is written as:
        #
        #     [,1,0,0,-1,1,0,0]
        #
        # The empty entry is not an omitted zero. It is an extra delimiter.
        # Replacing it with zero would create eight columns for seven reaction
        # rates and would make every solve fail inside the LSODA callback.
        # Removing the delimiter gives the intended seven-column row:
        #
        #     [1,0,0,-1,1,0,0]
        #
        # Only the generated Python cache is repaired. The SBML mechanism
        # remains the source of truth and is never changed here.
        # ============================================================
        # Read the complete generated Python file as text.
        text = self.generated_file.read_text(encoding="utf-8")
        # Collect every original or corrected line in its original order.
        repaired_lines = []
        # Remember whether a correction was actually made.
        changed = False

        # Examine the generated file one line at a time while retaining newlines.
        for line in text.splitlines(keepends=True):
            # Only the line defining the stoichiometric matrix needs this repair.
            if "stoichiometricMatrix = np.array" in line:
                # Remove the empty first entry without adding another matrix column.
                repaired = line.replace("[,", "[")
                # Record whether this line differed from the generated version.
                changed = changed or repaired != line
                # Keep the corrected version of this line.
                line = repaired
            # Add the current line to the reconstructed generated file.
            repaired_lines.append(line)

        # Rewrite the generated file only if a bad delimiter was found.
        if changed:
            self.generated_file.write_text(
                "".join(repaired_lines),
                encoding="utf-8",
            )
            # Tell the user that the automatic compatibility repair was applied.
            print(
                f"Removed extra leading delimiters from generated stoichiometric "
                f"matrix: {self.generated_file.name}"
            )

    def _import_generated_class(self):
        # Import the generated Python file directly from its path. This avoids
        # requiring the generated_models directory to be a Python package.
        # Make a unique internal module name so different SBML files cannot clash.
        module_key = f"_generated_{self.generated_file.stem}_{abs(hash(self.generated_file))}"
        # Create an import description that points at the generated Python file.
        spec = importlib.util.spec_from_file_location(module_key, self.generated_file)

        # Stop with a clear message if Python cannot create a usable import loader.
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not import generated model: {self.generated_file}")

        # Create an empty module object using the import description.
        module = importlib.util.module_from_spec(spec)
        # Register the module so Python treats it like a normal imported module.
        sys.modules[module_key] = module
        # Execute the generated file and populate the module object.
        spec.loader.exec_module(module)

        # Try to retrieve the expected generated model class from the module.
        try:
            return getattr(module, self.class_name)
        # Explain which class is absent instead of returning a vague attribute error.
        except AttributeError as exc:
            raise ImportError(
                f"Generated class {self.class_name!r} was not found in "
                f"{self.generated_file.name}."
            ) from exc

    @staticmethod
    def _find_state_species_order(model_instance) -> list[str]:
        # Read the exact species order used by generated _SolveReactions().
        #
        # SBMLtoODEpy writes an assignment such as:
        #
        #     self.s['A'].amount, self.s['B'].amount = y
        #
        # Reading that assignment is safer than assuming dictionary order is
        # the same as the ODE state-vector order.
        # Read the generated _SolveReactions method as ordinary source text.
        source = inspect.getsource(model_instance._SolveReactions)
        # Keep the assignment's left side, which lists species in state order.
        left_side = source.split("= y", 1)[0]
        # Extract every species ID appearing before its .amount attribute.
        ids = re.findall(r"self\.s\[['\"]([^'\"]+)['\"]\]\.amount", left_side)

        # A missing list means the generated layout is not one this adapter knows.
        if not ids:
            raise RuntimeError(
                "Could not determine the SBMLtoODEpy state order from "
                "_SolveReactions()."
            )

        # Return species IDs in exactly the order used by the generated equations.
        return ids

    def new_instance(self):
        # Return a fresh model so one experimental condition cannot change the
        # state used by another experimental condition.
        # Calling the stored class constructs a completely new model object.
        return self._model_class()

    @staticmethod
    def _set_parameter(model, parameter_id: str, value: float) -> None:
        # Reject misspelled or unavailable parameters before attempting a fit.
        if parameter_id not in model.p:
            raise KeyError(f"Parameter {parameter_id!r} is absent from generated model.")

        # SBML kinetic parameters are normally marked constant. Constant means
        # constant during one simulation, not forbidden from being estimated.
        # SBMLtoODEpy's public setter blocks changes to constant parameters, so
        # fitting intentionally updates the stored value directly before solve_ivp.
        # Convert the fitted number to float and store it in the generated model.
        model.p[parameter_id]._value = float(value)

    @staticmethod
    def _set_species_concentration(model, species_id: str, value: float) -> None:
        # Reject misspelled or unavailable species before starting a simulation.
        if species_id not in model.s:
            raise KeyError(f"Species {species_id!r} is absent from generated model.")

        # Retrieve the generated Species object that must be initialized.
        species = model.s[species_id]
        # Convert spreadsheet or NumPy scalar input into an ordinary Python float.
        concentration = float(value)
        # Store the requested initial concentration inside the Species object.
        species._concentration = concentration
        # Convert concentration to amount because the ODE state stores amounts.
        species._amount = concentration * species.compartment.size

    def prepare_instance(self, initial_values, parameters, fixed_parameters=None):
        # Create and configure one independent model for one data condition.
        # Begin with a clean generated model whose state has not been simulated.
        model = self.new_instance()

        # Apply every condition-specific initial species concentration.
        for species_id, value in initial_values.items():
            self._set_species_concentration(model, species_id, value)

        # Apply parameters currently being varied by the optimizer.
        for parameter_id, value in parameters.items():
            self._set_parameter(model, parameter_id, value)

        # Apply model-specific values that remain fixed during optimization.
        for parameter_id, value in (fixed_parameters or {}).items():
            self._set_parameter(model, parameter_id, value)

        # Reset simulation time so every experimental condition starts at zero.
        model.time = 0.0
        # Evaluate any generated assignment rules once after initialization.
        model.AssignmentRules()
        # Give the fully configured independent model back to the caller.
        return model

    def initial_amount_vector(self, model) -> np.ndarray:
        # Collect species amounts in the exact order expected by _SolveReactions.
        return np.asarray(
            [model.s[species_id].amount for species_id in self.state_species_ids],
            dtype=float,
        )

    def rhs(self, model, t: float, y: np.ndarray) -> np.ndarray:
        # Convert SBMLtoODEpy's odeint calling order, (y, t), into the calling
        # order used by solve_ivp, which is (t, y).
        # Evaluate the generated biochemical equations and expose useful errors.
        try:
            derivative = np.asarray(model._SolveReactions(y, t), dtype=float)
        # A matrix multiplication error usually means generated dimensions disagree.
        except ValueError as exc:
            # Give a useful error before SciPy reduces the original matrix
            # problem to "Call-back ... failed" and "capi_return is NULL".
            # Replace SciPy's cryptic callback message with a corrective action.
            if "matmul" in str(exc):
                raise RuntimeError(
                    f"{self.sbml_file.name}: generated stoichiometric matrix "
                    "does not match the reaction-velocity vector. Delete the "
                    "generated_models folder so the corrected adapter can "
                    "regenerate and repair the model."
                ) from exc
            # Preserve unrelated ValueErrors exactly as they were originally raised.
            raise

        # The derivative must contain one value for every state-vector value.
        if derivative.shape != y.shape:
            raise RuntimeError(
                f"Generated derivative has shape {derivative.shape}, expected {y.shape}."
            )

        # Return dy/dt to solve_ivp so it can advance the numerical solution.
        return derivative

    def simulate(
        self,
        initial_values,
        parameters,
        t_eval,
        observable_id,
        fixed_parameters=None,
        method="LSODA",
        rtol=1e-6,
        atol=1e-6,
        return_output_rate=False,
    ):
        # Simulate at the exact requested times and return the observable.
        #
        # The generated ODE state stores species amounts. The experimental data
        # use concentration-like values, so the calculated amount is divided by
        # its compartment size before this method returns it.
        # Convert the requested output times into a one-dimensional NumPy array.
        t_eval = np.asarray(t_eval, dtype=float)

        # Require a usable time series containing at least a start and end time.
        if t_eval.ndim != 1 or len(t_eval) < 2:
            raise ValueError("t_eval must contain at least two one-dimensional times.")

        # Require every requested time to be later than the preceding time.
        if np.any(np.diff(t_eval) <= 0):
            raise ValueError("t_eval must be strictly increasing.")

        # Confirm that the requested experimental output is a simulated species.
        if observable_id not in self.state_index:
            raise KeyError(f"Observable {observable_id!r} is not a dynamic species.")

        # Build an independent model configured for this experimental condition.
        model = self.prepare_instance(
            initial_values=initial_values,
            parameters=parameters,
            fixed_parameters=fixed_parameters,
        )
        # Create the initial ODE state vector in generated-model species order.
        y0 = self.initial_amount_vector(model)

        # Numerically integrate all biochemical ODEs over the requested time range.
        solution = solve_ivp(
            # Adapt solve_ivp's (time, state) call to the generated model equations.
            fun=lambda t, y: self.rhs(model, t, y),
            # Integrate from the first requested time through the final requested time.
            t_span=(float(t_eval[0]), float(t_eval[-1])),
            # Supply the starting amount of every dynamic species.
            y0=y0,
            # Ask the solver to report values at the experimental measurement times.
            t_eval=t_eval,
            # Use the selected numerical integration algorithm, normally LSODA.
            method=method,
            # Set the allowed relative numerical integration error.
            rtol=float(rtol),
            # Set the allowed absolute numerical integration error.
            atol=float(atol),
        )

        # Do not allow an incomplete numerical trajectory to enter model fitting.
        if not solution.success:
            raise RuntimeError(solution.message)

        # Find the solution-array row belonging to the requested observable.
        observable_index = self.state_index[observable_id]
        # Read the physical compartment volume used by that observable species.
        compartment_size = model.s[observable_id].compartment.size
        # Convert the simulated observable amount back into concentration units.
        output = solution.y[observable_index] / compartment_size

        # Return only the observable trajectory unless its rate was also requested.
        if not return_output_rate:
            return output

        # Allocate an output array with the same shape as the observable trajectory.
        output_rate = np.empty_like(output)
        # Re-evaluate dy/dt at every reported solution time.
        for index, time in enumerate(solution.t):
            # Calculate all species derivatives at this saved state and time.
            derivative = self.rhs(model, float(time), solution.y[:, index])
            # Select the observable derivative and convert it to concentration/time.
            output_rate[index] = derivative[observable_index] / compartment_size

        # Return both the observable concentration and its instantaneous rate.
        return output, output_rate
