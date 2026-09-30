"""Scientific checks for the built-in steady-state collector (no SBML files)."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "properties" / "steady_state_simulation.py"
SPEC = importlib.util.spec_from_file_location("steady_state_under_test", MODULE_PATH)
steady_state = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(steady_state)


class TinySimulator:
    """Minimal adapter double: equations are supplied in concentration units."""

    def __init__(self, ids, concentration_rhs, volumes=None):
        self.state_species_ids = list(ids)
        self.state_index = {name: index for index, name in enumerate(ids)}
        self.volumes = np.asarray(volumes or [1.0] * len(ids), dtype=float)
        self.concentration_rhs = concentration_rhs
        self.prepared_parameters = []

    def prepare_instance(self, initial_values, parameters, fixed_parameters=None):
        prepared = dict(parameters)
        prepared.update(fixed_parameters or {})
        self.prepared_parameters.append(prepared.copy())
        species = {
            name: SimpleNamespace(compartment=SimpleNamespace(size=self.volumes[index]))
            for index, name in enumerate(self.state_species_ids)
        }
        initial = np.asarray([initial_values.get(name, 0.0)
                              for name in self.state_species_ids])
        return SimpleNamespace(s=species, parameters=prepared,
                               initial_amounts=initial * self.volumes)

    def initial_amount_vector(self, model):
        return model.initial_amounts.copy()

    def rhs(self, model, time, amounts):
        concentration_rates = np.asarray(
            self.concentration_rhs(time, amounts / self.volumes, model.parameters),
            dtype=float,
        )
        return concentration_rates * self.volumes


def birth_death(time, concentrations, parameters):
    return [parameters["production"] - parameters["decay"] * concentrations[0]]


class TestSteadyStateSimulation(unittest.TestCase):

    def settings(self, **changes):
        settings = {
            "max_time_min": 40.0,
            "check_interval_min": 2.0,
            "rate_tolerance": 1e-8,
            "change_atol": 1e-7,
            "change_rtol": 1e-7,
            "consecutive_checks": 2,
            "samples_per_check": 5,
        }
        settings.update(changes)
        return settings

    def run_model(self, simulator, initial=None, parameters=None, settings=None,
                  **kwargs):
        return steady_state.find_steady_state(
            simulator, initial or {}, parameters or {}, "x",
            self.settings() if settings is None else settings,
            model_time_units_per_minute=1.0, rtol=1e-10, atol=1e-12,
            **kwargs
        )

    def test_birth_death_matches_analytic_equilibrium(self):
        simulator = TinySimulator(["x"], birth_death)
        result = self.run_model(simulator, {"x": 0.0},
                                {"production": 6.0, "decay": 2.0})
        self.assertTrue(result["converged"])
        self.assertEqual(result["status"], "converged")
        self.assertAlmostEqual(result["value"], 3.0, places=8)
        self.assertLessEqual(result["max_abs_rate"], 1e-8)
        self.assertLessEqual(result["max_scaled_change"], 1.0)
        self.assertEqual(len(simulator.prepared_parameters), 1)

    def test_hidden_accumulation_prevents_false_equilibrium(self):
        simulator = TinySimulator(["x", "hidden"],
                                  lambda t, c, p: [0.0, 1.0])
        result = self.run_model(simulator, {"x": 5.0}, settings=self.settings(
            max_time_min=6.0))
        self.assertFalse(result["converged"])
        self.assertTrue(np.isnan(result["value"]))
        self.assertEqual(result["final_output"], 5.0)
        self.assertEqual(result["worst_rate_species"], "hidden")
        self.assertAlmostEqual(result["max_abs_rate"], 1.0)
        self.assertAlmostEqual(result["final_concentrations"][1], 6.0)

    def test_compartment_volume_preserves_concentrations_and_rate_threshold(self):
        results = []
        for volume in [1.0, 1e-3]:
            simulator = TinySimulator(["x"], birth_death, volumes=[volume])
            results.append(self.run_model(simulator, parameters={
                "production": 4.0, "decay": 2.0}))
        for result in results:
            self.assertTrue(result["converged"])
            self.assertAlmostEqual(result["value"], 2.0, places=8)
        self.assertEqual(results[0]["time_minutes"], results[1]["time_minutes"])
        self.assertAlmostEqual(results[0]["max_abs_rate"],
                               results[1]["max_abs_rate"], places=9)

    def test_fixed_parameters_are_applied_by_adapter(self):
        simulator = TinySimulator(["x"], birth_death)
        result = self.run_model(simulator, parameters={"production": 1, "decay": 3},
                                fixed_parameters={"production": 6})
        self.assertTrue(result["converged"])
        self.assertAlmostEqual(result["value"], 2.0, places=8)
        self.assertEqual(simulator.prepared_parameters,
                         [{"production": 6, "decay": 3}])

    def test_timeout_keeps_finite_endpoint_but_not_a_steady_state(self):
        simulator = TinySimulator(["x"], birth_death)
        result = self.run_model(simulator, parameters={"production": 1, "decay": 0.01},
                                settings=self.settings(max_time_min=4.0))
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"], "max_time_reached")
        self.assertTrue(np.isnan(result["value"]))
        self.assertTrue(np.isfinite(result["final_output"]))
        self.assertAlmostEqual(result["final_output"],
                               100 * (1 - np.exp(-0.04)), places=7)
        self.assertEqual(result["time_minutes"], 4.0)

    def test_negative_concentration_is_not_an_equilibrium(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [-1.0])
        result = self.run_model(simulator, {"x": 0.5})
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"], "negative_concentration")
        self.assertTrue(np.isnan(result["value"]))
        self.assertLess(result["final_output"], 0.0)

    def test_negative_initial_concentration_is_rejected(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.0])
        with self.assertRaisesRegex(ValueError, "Initial concentrations"):
            self.run_model(simulator, {"x": -1.0})

    def test_consecutive_checks_require_multiple_full_windows(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.0])
        result = self.run_model(simulator, {"x": 7.0}, settings=self.settings(
            consecutive_checks=3, check_interval_min=2.0))
        self.assertTrue(result["converged"])
        self.assertEqual(result["time_minutes"], 6.0)
        self.assertEqual(result["value"], 7.0)

    def test_absolute_rate_prevents_large_accumulation_passing_relative_check(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.01])
        result = self.run_model(simulator, {"x": 1e8}, settings=self.settings(
            max_time_min=4.0, change_atol=1.0, change_rtol=0.01))
        self.assertLess(result["max_scaled_change"], 1.0)
        self.assertGreater(result["max_abs_rate"], 1e-8)
        self.assertFalse(result["converged"])

    def test_invalid_settings_raise_clear_value_errors(self):
        cases = [
            {"max_time_min": 0},
            {"max_time_min": 1, "check_interval_min": 2},
            {"rate_tolerance": -1},
            {"rate_tolerance": float("nan")},
            {"change_atol": 0},
            {"change_rtol": -1},
            {"samples_per_check": 2},
            {"samples_per_check": True},
            {"consecutive_checks": 1},
            {"consecutive_checks": 2.5},
        ]
        for changes in cases:
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    steady_state.read_steady_state_settings(self.settings(**changes))

    def test_invalid_compartment_volume_is_rejected(self):
        for volume in [0.0, -1.0, float("nan")]:
            with self.subTest(volume=volume):
                simulator = TinySimulator(["x"], lambda t, c, p: [0.0], [volume])
                with self.assertRaisesRegex(ValueError, "compartment sizes"):
                    self.run_model(simulator)


if __name__ == "__main__":
    unittest.main()
