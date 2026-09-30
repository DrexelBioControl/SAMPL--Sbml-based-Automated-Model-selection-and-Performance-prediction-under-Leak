"""Numerical checks for metric formulae and external experiment setup."""

import copy
import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_metric(filename):
    path = ROOT / "metrics" / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AMPLIFICATION = load_metric("amplification.py")
TEMPLATE = load_metric("MetricTemplate.py")
AMP_SETTINGS = {
    "input_values": {"on_value": 2.5, "off_value": 0.0, "label": "Input", "units": "nM"},
    "reference_value": 0.0,
}


def amplification(properties, experiment, settings):
    return AMPLIFICATION.calculate_metric(properties, experiment, dict(AMP_SETTINGS, **settings))


def sample_inputs():
    # The second control doubles discrimination relative to the reference.
    properties = {
        "on_rates": {"values": np.array([[4., 4., 4.], [8., 8., 8.]]), "units": "nM/time"},
        "off_rates": {"values": np.array([[2., 2., 2.], [2., 2., 2.]]), "units": "nM/time"},
        "reference_on": {"values": np.array([[4., 4., 4.]]), "units": "nM/time"},
        "reference_off": {"values": np.array([[2., 2., 2.]]), "units": "nM/time"},
        "output": {"values": np.array([[0., 4., 8.], [0., 8., 16.]]), "units": "nM"},
    }
    experiment = {
        "time_minutes": np.array([0.0, 1.0, 2.0]),
        "control_values": np.array([0.0, 5.0]),
        "control_units": "nM",
        "output_units": "nM",
        "input_control": {"type": "species_initial_concentration", "id": "IN"},
    }
    return properties, experiment


class PerformanceMetricTests(unittest.TestCase):
    def test_known_amplification_and_reference(self):
        result = amplification(*sample_inputs(), {"rate_floor": 1.0})
        np.testing.assert_array_equal(result["values"], [[1, 1, 1], [2, 2, 2]])
        np.testing.assert_array_equal(result["reference_values"], [1, 1, 1])
        self.assertEqual(result["kind"], "time_series")
        self.assertEqual(result["units"], "")
        self.assertEqual(result["reference_control_value"], 0.0)
        self.assertEqual(result["reference_label"], "Reference: 0 nM")
        self.assertEqual(result["condition_label"], "Input ON = 2.5 nM; OFF = 0 nM")

    def test_each_of_the_four_rates_can_make_a_point_undefined(self):
        properties, experiment = sample_inputs()
        experiment["control_values"] = np.array([5.0])
        experiment["time_minutes"] = np.arange(5, dtype=float)
        properties["on_rates"]["values"] = np.array([[0.5, 2, 2, 2, 2]])
        properties["off_rates"]["values"] = np.array([[2, 0.5, 2, 2, 2]])
        properties["reference_on"]["values"] = np.array([[2, 2, 0.5, 2, 2]])
        properties["reference_off"]["values"] = np.array([[2, 2, 2, 0.5, 2]])
        result = amplification(properties, experiment, {"rate_floor": 1.0})
        np.testing.assert_array_equal(
            np.isnan(result["values"]), [[True, True, True, True, False]]
        )
        self.assertEqual(result["values"][0, -1], 1.0)
        np.testing.assert_array_equal(
            np.isnan(result["reference_values"]), [False, False, True, True, False]
        )

    def test_floor_is_inclusive_and_rates_keep_their_sign(self):
        properties, experiment = sample_inputs()
        properties["on_rates"]["values"][1] = [1.0, -1.0, 0.999]
        result = amplification(properties, experiment, {"rate_floor": 1.0})
        np.testing.assert_allclose(result["values"][1, :2], [0.25, -0.25])
        self.assertTrue(np.isnan(result["values"][1, 2]))

    def test_zero_and_nan_are_undefined(self):
        properties, experiment = sample_inputs()
        properties["on_rates"]["values"][1] = [0.0, np.nan, 8.0]
        result = amplification(properties, experiment, {})
        np.testing.assert_array_equal(np.isnan(result["values"][1]), [True, True, False])
        self.assertEqual(result["values"][1, 2], 2.0)

    def test_reference_does_not_need_to_be_in_control_grid(self):
        properties, experiment = sample_inputs()
        experiment["control_values"] = np.array([5.0])
        for name in ("on_rates", "off_rates"):
            properties[name]["values"] = properties[name]["values"][1:]
        result = amplification(properties, experiment, {})
        np.testing.assert_array_equal(result["values"], [[2, 2, 2]])
        np.testing.assert_array_equal(result["reference_values"], [1, 1, 1])

    def test_invalid_floor_is_rejected(self):
        for rate_floor in [0.0, -1.0, np.nan, np.inf]:
            with self.subTest(rate_floor=rate_floor):
                with self.assertRaises(ValueError):
                    amplification({}, {}, {"rate_floor": rate_floor})

    def test_metric_does_not_change_property_or_experiment_arrays(self):
        properties, experiment = sample_inputs()
        original = copy.deepcopy(properties)
        original_experiment = copy.deepcopy(experiment)
        amplification(properties, experiment, {})
        for name in properties:
            np.testing.assert_array_equal(properties[name]["values"], original[name]["values"])
        for key in ("control_values", "time_minutes"):
            np.testing.assert_array_equal(experiment[key], original_experiment[key])

    def test_amplification_requests_only_its_four_rate_conditions(self):
        _, experiment = sample_inputs()
        requests = AMPLIFICATION.build_requests(experiment, AMP_SETTINGS)
        self.assertEqual(set(requests), {"on_rates", "off_rates", "reference_on", "reference_off"})
        self.assertTrue(all(request["property"] == "rates" for request in requests.values()))
        self.assertEqual(requests["on_rates"]["initial_values"], {"IN": 2.5})
        self.assertEqual(requests["off_rates"]["initial_values"], {"IN": 0.0})
        self.assertNotIn("control_values", requests["on_rates"])
        for name in ("reference_on", "reference_off"):
            self.assertEqual(requests[name]["control_values"], [0.0])

    def test_parameter_input_is_an_explicit_request_override(self):
        _, experiment = sample_inputs()
        experiment["input_control"] = {"type": "parameter", "id": "stimulus"}
        requests = AMPLIFICATION.build_requests(experiment, AMP_SETTINGS)
        self.assertEqual(requests["on_rates"]["parameters"], {"stimulus": 2.5})
        self.assertEqual(requests["off_rates"]["parameters"], {"stimulus": 0.0})

    def test_autonomous_model_needs_no_input_override(self):
        properties, experiment = sample_inputs()
        experiment["input_control"] = None
        requests = AMPLIFICATION.build_requests(experiment, AMP_SETTINGS)
        for request in requests.values():
            self.assertNotIn("initial_values", request)
            self.assertNotIn("parameters", request)
        result = amplification(properties, experiment, {})
        self.assertEqual(result["condition_label"], "Configured initial conditions")

    def test_setup_requires_explicit_finite_input_and_reference_values(self):
        _, experiment = sample_inputs()
        invalid_settings = [
            {}, {"input_values": {"on_value": 1}, "reference_value": 0},
            {"input_values": AMP_SETTINGS["input_values"]},
            dict(AMP_SETTINGS, reference_value=np.nan),
            dict(AMP_SETTINGS, input_values={"on_value": np.inf, "off_value": 0}),
        ]
        for settings in invalid_settings:
            with self.subTest(settings=settings):
                with self.assertRaises(ValueError):
                    AMPLIFICATION.build_requests(experiment, settings)

    def test_setup_leaves_settings_and_experiment_unchanged(self):
        _, experiment = sample_inputs()
        settings = copy.deepcopy(AMP_SETTINGS)
        saved = copy.deepcopy(settings)
        original = copy.deepcopy(experiment)
        AMPLIFICATION.build_requests(experiment, settings)
        self.assertEqual(settings, saved)
        self.assertEqual(experiment["input_control"], original["input_control"])
        np.testing.assert_array_equal(experiment["control_values"], original["control_values"])

    def test_template_needs_only_one_ordinary_condition(self):
        self.assertEqual(TEMPLATE.build_requests({}, {}), {"output": {"property": "concentrations"}})
        properties, experiment = sample_inputs()
        result = TEMPLATE.calculate_metric({"output": properties["output"]}, experiment, {})
        np.testing.assert_array_equal(result["values"], [8, 16])
        self.assertNotIn("reference_values", result)
        self.assertEqual(result["kind"], "per_control")
        self.assertEqual(result["units"], "nM")
        result["values"][0] = 99
        self.assertEqual(properties["output"]["values"][0, -1], 8)


if __name__ == "__main__":
    unittest.main()

