"""Checks for output plateaus, horizon extension, and full-system safety."""

from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from test_steady_state_simulation import TinySimulator, steady_state
from test_sensitivity_metric import METRIC, make_inputs


class OutputPlateauTests(unittest.TestCase):

    def settings(self, **changes):
        settings = {
            "convergence_mode": "output_plateau",
            "plateau_extension_factor": 2.0,
            "max_time_min": 40.0,
            "check_interval_min": 2.0,
            "rate_tolerance": 1e-8,
            "change_atol": 1e-7,
            "change_rtol": 1e-7,
            "consecutive_checks": 2,
            "samples_per_check": 9,
        }
        settings.update(changes)
        return settings

    def run_model(self, simulator, initial=None, settings=None):
        return steady_state.find_steady_state(
            simulator, {"x": 5.0} if initial is None else initial, {}, "x",
            self.settings() if settings is None else settings,
            model_time_units_per_minute=1.0, rtol=1e-10, atol=1e-12,
        )

    def test_exact_plateau_requires_extension_beyond_initial_checks(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.0])
        result = self.run_model(simulator)
        self.assertTrue(result["converged"])
        self.assertTrue(result["full_system_converged"])
        self.assertEqual(result["status"], "output_plateau_verified")
        self.assertEqual(result["plateau_start_time_minutes"], 4.0)
        self.assertEqual(result["plateau_target_time_minutes"], 8.0)
        self.assertEqual(result["time_minutes"], 8.0)
        self.assertEqual(result["value"], 5.0)
        self.assertEqual(result["plateau_max_scaled_change"], 0.0)

    def test_extension_requires_at_least_another_set_of_check_windows(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.0])
        result = self.run_model(simulator, settings=self.settings(
            plateau_extension_factor=1.1))
        # factor * candidate_time = 4.4, but two extra windows require 8.
        self.assertEqual(result["plateau_target_time_minutes"], 8.0)
        self.assertEqual(result["time_minutes"], 8.0)

    def test_cap_is_not_extended_silently_to_finish_a_candidate(self):
        for maximum, expected in [(6.0, False), (8.0, True)]:
            with self.subTest(maximum=maximum):
                simulator = TinySimulator(["x"], lambda t, c, p: [0.0])
                result = self.run_model(simulator, settings=self.settings(
                    max_time_min=maximum))
                self.assertEqual(result["converged"], expected)
                self.assertEqual(result["time_minutes"], maximum)
                self.assertEqual(result["plateau_target_time_minutes"], 8.0)
                if not expected:
                    self.assertEqual(result["status"], "plateau_extension_incomplete")
                    self.assertTrue(np.isnan(result["value"]))
                    self.assertEqual(result["final_output"], 5.0)

    def test_hidden_accumulation_passes_only_output_mode(self):
        for mode, expected in [("full_system", False), ("output_plateau", True)]:
            with self.subTest(mode=mode):
                simulator = TinySimulator(["x", "hidden"],
                                          lambda t, c, p: [0.0, 1.0])
                result = self.run_model(simulator, settings=self.settings(
                    convergence_mode=mode, max_time_min=8.0))
                self.assertEqual(result["converged"], expected)
                self.assertFalse(result["full_system_converged"])
                self.assertEqual(result["worst_rate_species"], "hidden")
                self.assertAlmostEqual(result["max_abs_rate"], 1.0)
                self.assertGreater(result["max_scaled_change"], 1.0)
                self.assertEqual(result["output_max_abs_rate"], 0.0)
                self.assertEqual(result["output_max_scaled_change"], 0.0)

    def test_delayed_departure_invalidates_initial_plateau(self):
        simulator = TinySimulator(
            ["x"], lambda t, c, p: [0.0 if t < 5.0 else 1.0])
        result = self.run_model(simulator, settings=self.settings(max_time_min=8.0))
        self.assertFalse(result["converged"])
        self.assertTrue(np.isnan(result["value"]))
        self.assertGreater(result["final_output"], 7.0)

    def test_new_plateau_after_departure_gets_its_own_extension(self):
        simulator = TinySimulator(
            ["x"], lambda t, c, p: [1.0 if 5.0 <= t < 9.0 else 0.0])
        result = self.run_model(simulator)
        self.assertTrue(result["converged"])
        self.assertEqual(result["status"], "output_plateau_verified")
        self.assertGreater(result["plateau_start_time_minutes"], 9.0)
        self.assertEqual(result["plateau_target_time_minutes"],
                         2.0 * result["plateau_start_time_minutes"])
        self.assertGreaterEqual(result["time_minutes"],
                                result["plateau_target_time_minutes"])
        self.assertAlmostEqual(result["value"], 9.0, places=6)

    def test_excursion_returning_to_same_endpoint_is_not_a_plateau(self):
        def transient(time, concentrations, parameters):
            rate = (np.sin(2.0 * np.pi * (time - 4.4) / 3.2)
                    if 4.4 < time < 7.6 else 0.0)
            return [rate]

        result = self.run_model(TinySimulator(["x"], transient),
                                settings=self.settings(max_time_min=8.0))
        self.assertAlmostEqual(result["final_output"], 5.0, places=6)
        self.assertFalse(result["converged"])
        self.assertTrue(np.isnan(result["value"]))

    def test_plateau_horizon_uses_minutes_and_rates_use_model_time(self):
        simulator = TinySimulator(["x", "hidden"],
                                  lambda t, c, p: [0.0, 1.0], volumes=[0.001, 100.0])
        result = steady_state.find_steady_state(
            simulator, {"x": 5.0}, {}, "x", self.settings(),
            model_time_units_per_minute=60.0, rtol=1e-10, atol=1e-12,
        )
        self.assertTrue(result["converged"])
        self.assertEqual(result["time_minutes"], 8.0)
        self.assertEqual(result["plateau_target_time_minutes"], 8.0)
        self.assertAlmostEqual(result["value"], 5.0)
        self.assertAlmostEqual(result["max_abs_rate"], 1.0)
        self.assertAlmostEqual(result["final_concentrations"][1], 480.0)

    def test_cumulative_drift_rejects_locally_flat_windows(self):
        simulator = TinySimulator(["x"], lambda t, c, p: [0.3])
        result = self.run_model(simulator, settings=self.settings(
            max_time_min=8.0, rate_tolerance=1.0,
            change_atol=1.0, change_rtol=0.0))
        # Each 2-minute window drifts only 0.6; the extension drifts 1.2.
        self.assertLess(result["output_max_scaled_change"], 1.0)
        self.assertLess(result["output_max_abs_rate"], 1.0)
        self.assertFalse(result["converged"])
        self.assertTrue(np.isnan(result["value"]))

    def test_negative_hidden_species_blocks_output_plateau(self):
        simulator = TinySimulator(["x", "hidden"], lambda t, c, p: [0.0, -1.0])
        result = self.run_model(simulator, initial={"x": 5.0, "hidden": 0.5})
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"], "negative_concentration")
        self.assertTrue(np.isnan(result["value"]))

    def test_invalid_late_state_cannot_keep_stale_full_system_flag(self):
        simulator = TinySimulator(
            ["x", "hidden"], lambda t, c, p: [0.0, 0.0 if t < 5.0 else -1.0])
        result = self.run_model(simulator, initial={"x": 5.0, "hidden": 0.5})
        self.assertEqual(result["status"], "negative_concentration")
        self.assertFalse(result["converged"])
        self.assertFalse(result["full_system_converged"])
        self.assertTrue(np.isnan(result["value"]))

    def test_nonfinite_hidden_species_blocks_output_plateau(self):
        simulator = TinySimulator(["x", "hidden"], lambda t, c, p: [0.0, 0.0])
        solution = SimpleNamespace(
            success=True, t=np.linspace(0.0, 2.0, 9),
            y=np.array([[5.0] * 9, [0.0] * 8 + [np.nan]]),
        )
        with patch.object(steady_state, "solve_ivp", return_value=solution):
            result = self.run_model(simulator)
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"], "nonfinite_concentration")
        self.assertTrue(np.isnan(result["value"]))

    def test_invalid_mode_and_extension_factor_raise_value_error(self):
        cases = [{"convergence_mode": mode} for mode in [None, "output", "Full_system"]]
        cases += [{"plateau_extension_factor": value}
                  for value in [True, 1.0, 0.0, -2.0, np.nan, np.inf]]
        for changes in cases:
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    steady_state.read_steady_state_settings(self.settings(**changes))

    def test_omitted_mode_keeps_full_system_behavior(self):
        settings = self.settings(convergence_mode="full_system")
        explicit = self.run_model(TinySimulator(["x"], lambda t, c, p: [0.0]),
                                  settings=settings)
        settings.pop("convergence_mode")
        implicit = self.run_model(TinySimulator(["x"], lambda t, c, p: [0.0]),
                                  settings=settings)
        for key in ["status", "converged", "time_minutes", "value",
                    "max_abs_rate", "max_scaled_change", "full_system_converged"]:
            self.assertEqual(implicit[key], explicit[key], key)
        self.assertEqual(implicit["status"], "converged")
        self.assertEqual(implicit["time_minutes"], 4.0)

    def test_metric_labels_distinguish_plateau_from_full_system_sensitivity(self):
        full_inputs = make_inputs([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
        full = METRIC.calculate_metric(*full_inputs, {})
        plateau_inputs = make_inputs(
            [1.0, 2.0, 3.0], [2.0, 4.0, 6.0], response_property="output_plateau"
        )
        plateau = METRIC.calculate_metric(
            *plateau_inputs, {"response_property": "output_plateau"}
        )
        self.assertIn("plateau", plateau["label"].lower())
        self.assertNotEqual(plateau["label"], full["label"])
        np.testing.assert_allclose(plateau["values"], full["values"])
        np.testing.assert_allclose(plateau["values"], 1.0)


if __name__ == "__main__":
    unittest.main()
