"""Numerical and validity checks for the external sensitivity definition."""

import copy
import importlib.util
from pathlib import Path
import unittest

import numpy as np


METRIC_FILE = Path(__file__).resolve().parents[1] / "metrics" / "sensitivity.py"
SPEC = importlib.util.spec_from_file_location("sensitivity_metric", str(METRIC_FILE))
METRIC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(METRIC)


def make_inputs(theta, xss, converged=None, response_property="steady_state"):
    theta = np.asarray(theta, dtype=float)
    properties = {
        "response": {
            "values": np.asarray(xss, dtype=float),
            "converged": (
                np.ones(theta.shape, dtype=bool)
                if converged is None else np.asarray(converged)
            ),
            "convergence_mode": (
                "full_system" if response_property == "steady_state" else "output_plateau"
            ),
            "units": "nM",
        },
    }
    experiment = {
        "control_values": theta,
        "analysis_control": {"type": "parameter", "id": "theta"},
        "observable_id": "x",
    }
    return properties, experiment


class SensitivityMetricTests(unittest.TestCase):
    def test_linear_response_has_unit_sensitivity_on_nonuniform_grid(self):
        theta = np.array([0.2, 0.5, 1.1, 2.0, 4.0])
        result = METRIC.calculate_metric(*make_inputs(theta, 3 * theta), {})
        np.testing.assert_allclose(result["values"], 1.0, rtol=1e-13)
        np.testing.assert_allclose(result["dxss_dtheta"], 3.0, rtol=1e-13)
        self.assertEqual(result["kind"], "per_control")
        self.assertEqual(result["units"], "")

    def test_quadratic_response_has_sensitivity_two_including_endpoints(self):
        theta = np.array([0.3, 0.6, 1.2, 1.8, 3.2])
        result = METRIC.calculate_metric(*make_inputs(theta, theta ** 2), {})
        np.testing.assert_allclose(result["values"], 2.0, rtol=1e-13)
        np.testing.assert_allclose(result["dxss_dtheta"], 2 * theta, rtol=1e-13)

    def test_constant_response_has_zero_sensitivity(self):
        theta = np.array([0.3, 0.6, 1.2, 1.8, 3.2])
        result = METRIC.calculate_metric(*make_inputs(theta, np.full(5, 4.0)), {})
        np.testing.assert_allclose(result["values"], 0.0, atol=1e-13)

    def test_inverse_response_preserves_negative_sign_and_refines(self):
        errors = []
        for count in (21, 41):
            theta = np.linspace(1.0, 2.0, count)
            result = METRIC.calculate_metric(*make_inputs(theta, 5.0 / theta), {})
            self.assertTrue(np.all(result["values"] < 0))
            errors.append(np.max(np.abs(result["values"] + 1.0)))
        self.assertLess(errors[1], errors[0] / 3.0)
        self.assertLess(errors[1], 0.002)

    def test_absolute_option_does_not_change_derivative(self):
        theta = np.linspace(1, 2, 11)
        properties, experiment = make_inputs(theta, 1 / theta)
        signed = METRIC.calculate_metric(properties, experiment, {})
        absolute = METRIC.calculate_metric(properties, experiment, {"absolute": True})
        np.testing.assert_allclose(absolute["values"], -signed["values"])
        np.testing.assert_allclose(absolute["dxss_dtheta"], signed["dxss_dtheta"])
        self.assertIn("Absolute", absolute["label"])

    def test_failed_point_invalidates_interior_and_endpoint_stencils(self):
        theta = np.arange(1.0, 8.0)
        # Point 2 belongs to both the first endpoint stencil and points 1-3.
        # Point 4 similarly belongs to points 3-5 and the last endpoint.
        for bad_index, invalid_indices in ((2, [0, 1, 2, 3]), (4, [3, 4, 5, 6])):
            with self.subTest(bad_index=bad_index):
                converged = np.ones(theta.size, dtype=bool)
                converged[bad_index] = False
                result = METRIC.calculate_metric(*make_inputs(theta, theta, converged), {})
                expected_invalid = np.zeros(theta.size, dtype=bool)
                expected_invalid[invalid_indices] = True
                np.testing.assert_array_equal(np.isnan(result["values"]), expected_invalid)
                np.testing.assert_array_equal(np.isnan(result["dxss_dtheta"]), expected_invalid)
                np.testing.assert_allclose(result["values"][~expected_invalid], 1.0)

    def test_invalid_concentration_does_not_get_bridged(self):
        theta = np.arange(1.0, 8.0)
        for bad_value in (0.0, -1.0, 1e-10, np.nan, np.inf):
            with self.subTest(bad_value=bad_value):
                xss = theta.copy()
                xss[3] = bad_value
                result = METRIC.calculate_metric(*make_inputs(theta, xss), {})
                np.testing.assert_array_equal(
                    np.isnan(result["values"]),
                    [False, False, True, True, True, False, False],
                )

    def test_floor_boundary_is_valid(self):
        theta = np.array([1.0, 2.0, 3.0])
        result = METRIC.calculate_metric(
            *make_inputs(theta, theta * 1e-8), {"concentration_floor": 1e-8}
        )
        np.testing.assert_allclose(result["values"], 1.0)

    def test_rejects_invalid_grid(self):
        grids = ([1, 2], [1, 1, 2], [3, 2, 1], [0, 1, 2], [-1, 1, 2],
                 [1, 2, np.inf], [1, 2, np.nan], [[1, 2, 3]])
        for theta in grids:
            with self.subTest(theta=theta):
                with self.assertRaises(ValueError):
                    METRIC.calculate_metric(*make_inputs(theta, np.ones(np.shape(theta))), {})

    def test_rejects_mismatched_arrays_and_nonboolean_convergence(self):
        for xss, flags in (([1, 2], [True] * 3),
                           ([1, 2, 3], [True] * 2),
                           ([1, 2, 3], [1, 1, 0])):
            with self.subTest(xss=xss, flags=flags):
                with self.assertRaises(ValueError):
                    METRIC.calculate_metric(*make_inputs([1, 2, 3], xss, flags), {})

    def test_rejects_nonparameter_control_and_invalid_settings(self):
        properties, experiment = make_inputs([1, 2, 3], [1, 2, 3])
        for floor in (0.0, -1.0, np.inf, np.nan):
            with self.subTest(floor=floor):
                with self.assertRaises(ValueError):
                    METRIC.calculate_metric(properties, experiment, {"concentration_floor": floor})
        for absolute in ("false", 1, None):
            with self.subTest(absolute=absolute):
                with self.assertRaises(ValueError):
                    METRIC.calculate_metric(properties, experiment, {"absolute": absolute})
        experiment["analysis_control"]["type"] = "species"
        with self.assertRaises(ValueError):
            METRIC.calculate_metric(properties, experiment, {})

    def test_requirements_select_explicit_response_property(self):
        self.assertEqual(METRIC.build_requests({}, {}), {"response": {"property": "steady_state"}})
        for name in ("steady_state", "output_plateau"):
            with self.subTest(name=name):
                settings = {"response_property": name}
                self.assertEqual(METRIC.build_requests({}, settings), {"response": {"property": name}})
                result = METRIC.calculate_metric(
                    *make_inputs([1, 2, 3], [2, 4, 6], response_property=name), settings
                )
                np.testing.assert_allclose(result["values"], 1.0)
        for name in ("rates", "concentrations", "unknown", None):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, "response_property"):
                    METRIC.build_requests({}, {"response_property": name})

    def test_requested_property_requires_matching_convergence_criterion(self):
        for name, wrong_mode in (("steady_state", "output_plateau"),
                                 ("output_plateau", "full_system"),
                                 ("steady_state", None)):
            with self.subTest(name=name, wrong_mode=wrong_mode):
                properties, experiment = make_inputs(
                    [1, 2, 3], [2, 4, 6], response_property=name
                )
                properties["response"]["convergence_mode"] = wrong_mode
                with self.assertRaisesRegex(ValueError, "requires convergence_mode"):
                    METRIC.calculate_metric(
                        properties, experiment, {"response_property": name}
                    )

    def test_inputs_are_unchanged(self):
        properties, experiment = make_inputs([1, 2, 3, 4], [2, 0, 6, 8])
        settings = {"absolute": True, "concentration_floor": 1e-6}
        saved_properties = copy.deepcopy(properties)
        saved_experiment = copy.deepcopy(experiment)
        saved_settings = settings.copy()
        METRIC.calculate_metric(properties, experiment, settings)
        for key in ("values", "converged"):
            np.testing.assert_array_equal(
                properties["response"][key], saved_properties["response"][key]
            )
        np.testing.assert_array_equal(
            experiment["control_values"], saved_experiment["control_values"]
        )
        self.assertEqual(experiment["analysis_control"], saved_experiment["analysis_control"])
        self.assertEqual(properties["response"]["convergence_mode"], "full_system")
        self.assertEqual(settings, saved_settings)


if __name__ == "__main__":
    unittest.main()
