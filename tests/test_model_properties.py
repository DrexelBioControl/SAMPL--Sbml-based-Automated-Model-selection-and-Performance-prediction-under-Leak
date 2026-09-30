"""Check generic property requests and independent numerical experiments."""

from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from properties import model_properties as properties
from properties import time_course_simulation as time_course


class LinearSimulator:
    """Small adapter stand-in: dX/dt = k + U, with a compartment volume of 2."""

    state_species_ids = ["X"]
    state_index = {"X": 0}

    def __init__(self):
        self.calls = []
        self.prepared = []

    def prepare_instance(self, initial_values, parameters, fixed_parameters):
        self.prepared.append((initial_values.copy(), parameters.copy(), fixed_parameters.copy()))
        return SimpleNamespace(initial=initial_values["X"],
                               rate=parameters["k"] + initial_values.get("U", 0.0),
                               s={"X": SimpleNamespace(compartment=SimpleNamespace(size=2.0))})

    def initial_amount_vector(self, model):
        return np.asarray([2.0 * model.initial])

    def rhs(self, model, time, state):
        return np.asarray([2.0 * model.rate])

    def simulate(self, **kwargs):
        self.calls.append(deepcopy(kwargs))
        initial = kwargs["initial_values"]
        rate = kwargs["parameters"]["k"] + initial.get("U", 0.0)
        times = kwargs["t_eval"]
        # Like the real adapter, the first sample defines the integration start.
        return initial["X"] + rate * (times - times[0]), np.full(times.shape, rate)


def example_experiment(times=(0.0, 5.0, 10.0)):
    # Deliberately contains no ON/OFF values or input-control mapping.
    return {
        "initial_values": {"X": 3.0},
        "analysis_control": {"type": "parameter", "id": "k"},
        "fixed_parameters": {}, "observable_id": "X",
        "control_values": np.asarray([1.0, 2.0]), "control_reference": 1.0,
        "time_minutes": np.asarray(times, dtype=float), "time_model": np.asarray(times, dtype=float),
        "model_time_units_per_minute": 1.0, "output_units": "nM",
        "rtol": 1e-8, "atol": 1e-10, "limit_settings": {},
    }


def requests_for(experiment, raw, definitions=None):
    metric = SimpleNamespace(build_requests=lambda experiment, settings: raw)
    return properties.read_property_requests(metric, experiment, {}, definitions)


class TestPropertyRequests(unittest.TestCase):

    def test_metric_selects_property_and_condition_using_experiment_and_settings(self):
        experiment = example_experiment()
        metric = SimpleNamespace(build_requests=lambda experiment, settings: {
            "response": {"property": settings["response_property"],
                         "initial_values": {"U": experiment["control_reference"]}}
        })
        for name in ("steady_state", "output_plateau"):
            result = properties.read_property_requests(metric, experiment, {"response_property": name})
            self.assertEqual(set(result), {"response"})
            self.assertEqual(result["response"]["property"], name)
            self.assertEqual(result["response"]["initial_values"], {"U": 1.0})
            self.assertEqual(result["response"]["parameters"], {})
            np.testing.assert_array_equal(result["response"]["control_values"], [1.0, 2.0])

    def test_normalized_requests_copy_overrides_and_default_or_explicit_grids(self):
        experiment = example_experiment()
        raw = {
            "sweep": {"property": "rates", "initial_values": {"U": 2.0}, "parameters": {"q": 4.0}},
            "reference": {"property": "rates", "control_values": [0.0]},
        }
        result = requests_for(experiment, raw)
        result["sweep"]["initial_values"]["U"] = 99.0
        result["sweep"]["parameters"]["q"] = 99.0
        result["sweep"]["control_values"][0] = 99.0
        result["reference"]["control_values"][0] = 99.0
        self.assertEqual(raw["sweep"]["initial_values"], {"U": 2.0})
        self.assertEqual(raw["sweep"]["parameters"], {"q": 4.0})
        self.assertEqual(raw["reference"]["control_values"], [0.0])
        np.testing.assert_array_equal(experiment["control_values"], [1.0, 2.0])

    def test_missing_noncallable_and_invalid_requests_fail(self):
        experiment = example_experiment()
        for metric in (SimpleNamespace(), SimpleNamespace(build_requests="rates")):
            with self.subTest(metric=vars(metric)), self.assertRaises(ValueError):
                properties.read_property_requests(metric, experiment, {})
        invalid = [
            None, [], {}, "rates", {"bad name": {"property": "rates"}},
            {"response": "rates"}, {"response": {}},
            {"response": {"property": "unknown_quantity"}},
            {"response": {"property": "bad name"}},
            {"response": {"property": "rates", "input_on": 1.0}},
            {"response": {"property": "rates", "control_values": []}},
            {"response": {"property": "rates", "control_values": [[1.0, 2.0]]}},
            {"response": {"property": "rates", "control_values": [np.nan]}},
            {"response": {"property": "rates", "control_values": [np.inf]}},
            {"response": {"property": "rates", "initial_values": []}},
            {"response": {"property": "rates", "initial_values": {"U": -1.0}}},
            {"response": {"property": "rates", "initial_values": {"U": np.nan}}},
            {"response": {"property": "rates", "parameters": []}},
            {"response": {"property": "rates", "parameters": {"q": np.inf}}},
        ]
        for raw in invalid:
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                requests_for(experiment, raw)

    def test_generic_grids_allow_zero_singletons_and_negative_parameter_values(self):
        for grid in ([0.0], [-1.0, 0.0, 1.0]):
            result = requests_for(example_experiment(), {
                "response": {"property": "rates", "control_values": grid, "parameters": {"q": -2.0}}
            })
            np.testing.assert_array_equal(result["response"]["control_values"], grid)

    def test_imports_cannot_shadow_builtin_properties(self):
        for name in properties.BUILTIN_PROPERTIES:
            with self.subTest(name=name), self.assertRaises(ValueError):
                requests_for(example_experiment(), {"response": {"property": "rates"}},
                             {name: {"definition_file": "unused.py"}})

    def test_explicit_times_and_regular_grid_validation(self):
        actual = time_course.read_time_values({"time_values_min": [5, 10], "time_end_min": 100, "time_points": 3})
        np.testing.assert_array_equal(actual, [5, 10])
        np.testing.assert_array_equal(time_course.read_time_values({"time_end_min": 10, "time_points": 3}), [0, 5, 10])
        self.assertEqual(time_course.read_time_values({}).size, 0)
        with self.assertRaises(ValueError):
            time_course.read_time_values({}, required=True)
        for times in ([], [-1, 0], [0, 0], [10, 5], [0, np.nan], [0, np.inf], [[0, 1]]):
            with self.subTest(times=times), self.assertRaises(ValueError):
                time_course.read_time_values({"time_values_min": times})


class TestRequestTargets(unittest.TestCase):

    def validate(self, raw, experiment=None, sweep_parameter=None):
        experiment = example_experiment() if experiment is None else experiment
        requests = requests_for(experiment, raw)
        properties.validate_request_targets(requests, experiment, ["X", "U", "S"], ["X"],
                                            ["k", "q", "p"], sweep_parameter=sweep_parameter)

    def test_valid_species_and_parameter_overrides(self):
        self.validate({"response": {"property": "rates", "initial_values": {"U": 2.0}, "parameters": {"q": 3.0}}})

    def test_unknown_species_and_parameters_fail(self):
        for field in ("initial_values", "parameters"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate({"response": {"property": "rates", field: {"missing": 1.0}}})

    def test_request_cannot_override_analysis_control(self):
        with self.assertRaises(ValueError):
            self.validate({"response": {"property": "rates", "parameters": {"k": 3.0}}})
        experiment = example_experiment()
        experiment["analysis_control"] = {"type": "species_initial_concentration", "id": "S"}
        with self.assertRaises(ValueError):
            self.validate({"response": {"property": "rates", "initial_values": {"S": 3.0}}}, experiment)

    def test_request_cannot_override_fixed_or_secondary_sweep_parameter(self):
        raw = {"response": {"property": "rates", "parameters": {"q": 3.0}}}
        experiment = example_experiment()
        experiment["fixed_parameters"] = {"q": 4.0}
        with self.assertRaises(ValueError):
            self.validate(raw, experiment)
        with self.assertRaises(ValueError):
            self.validate(raw, sweep_parameter="q")

    def test_observable_must_be_dynamic(self):
        experiment = example_experiment()
        experiment["observable_id"] = "U"
        with self.assertRaises(ValueError):
            self.validate({"response": {"property": "concentrations"}}, experiment)


class TestTimeProperties(unittest.TestCase):

    def test_three_arbitrary_conditions_and_reference_grid_have_separate_results(self):
        simulator = LinearSimulator()
        experiment = example_experiment()
        requests = requests_for(experiment, {
            "low": {"property": "concentrations", "initial_values": {"U": 0.0}},
            "medium": {"property": "concentrations", "initial_values": {"U": 2.0}},
            "high": {"property": "concentrations", "initial_values": {"U": 4.0}},
            "reference": {"property": "rates", "initial_values": {"U": 2.0}, "control_values": [0.0]},
        })
        result = properties.collect_properties(requests, simulator, experiment, {"k": 9.0})
        self.assertEqual(set(result), set(requests))
        self.assertEqual(len(simulator.calls), 7)
        for name, offset in (("low", 0.0), ("medium", 2.0), ("high", 4.0)):
            expected = 3.0 + (np.asarray([1.0, 2.0]) + offset)[:, None] * experiment["time_model"]
            np.testing.assert_array_equal(result[name]["values"], expected)
            np.testing.assert_array_equal(result[name]["control_values"], [1.0, 2.0])
        np.testing.assert_array_equal(result["reference"]["values"], [[2.0, 2.0, 2.0]])
        np.testing.assert_array_equal(result["reference"]["control_values"], [0.0])
        self.assertEqual(experiment["initial_values"], {"X": 3.0})

    def test_rates_and_concentrations_share_exactly_matching_trajectory_batch(self):
        simulator = LinearSimulator()
        experiment = example_experiment()
        requests = requests_for(experiment, {
            "output": {"property": "concentrations", "initial_values": {"U": 1.0}},
            "speed": {"property": "rates", "initial_values": {"U": 1.0}},
        })
        result = properties.collect_properties(requests, simulator, experiment, {"k": 9.0})
        self.assertEqual(len(simulator.calls), 2)
        np.testing.assert_array_equal(result["output"]["values"], [[3, 13, 23], [3, 18, 33]])
        np.testing.assert_array_equal(result["speed"]["values"], [[2, 2, 2], [3, 3, 3]])
        self.assertEqual(result["speed"]["units"], "nM / model time unit")
        self.assertEqual(result["output"]["units"], "nM")

    def test_distinct_tiny_grids_and_overrides_are_not_merged(self):
        simulator = LinearSimulator()
        experiment = example_experiment()
        requests = requests_for(experiment, {
            "first": {"property": "rates", "control_values": [1e-10, 2e-10]},
            "second": {"property": "rates", "control_values": [2e-10, 4e-10]},
            "third": {"property": "rates", "control_values": [1e-10, 2e-10], "initial_values": {"U": 1e-10}},
        })
        result = properties.collect_properties(requests, simulator, experiment, {"k": 9.0})
        self.assertEqual(len(simulator.calls), 6)
        for name, expected in (("first", [1e-10, 2e-10]), ("second", [2e-10, 4e-10]), ("third", [2e-10, 3e-10])):
            np.testing.assert_allclose(result[name]["values"][:, 0], expected, rtol=1e-15, atol=0)

    def test_samples_after_zero_integrate_from_the_initial_condition(self):
        for times in ([5.0, 10.0], [10.0]):
            simulator = LinearSimulator()
            experiment = example_experiment(times)
            experiment["initial_values"]["U"] = 3.0
            output, rate = time_course.simulate_time_points(
                simulator, experiment["initial_values"], {"k": 2.0},
                experiment["observable_id"], experiment["time_model"],
                fixed_parameters=experiment["fixed_parameters"],
                rtol=experiment["rtol"], atol=experiment["atol"],
            )
            np.testing.assert_array_equal(simulator.calls[0]["t_eval"], [0.0] + times)
            np.testing.assert_array_equal(output, 3.0 + 5.0 * np.asarray(times))
            np.testing.assert_array_equal(rate, np.full(len(times), 5.0))

    def test_zero_only_reads_initial_concentration_and_rhs_without_integration(self):
        simulator = LinearSimulator()
        experiment = example_experiment([0.0])
        experiment["initial_values"]["U"] = 3.0
        output, rate = time_course.simulate_time_points(
            simulator, experiment["initial_values"], {"k": 2.0},
            experiment["observable_id"], experiment["time_model"],
            fixed_parameters=experiment["fixed_parameters"],
            rtol=experiment["rtol"], atol=experiment["atol"],
        )
        np.testing.assert_array_equal(output, [3.0])
        np.testing.assert_array_equal(rate, [5.0])
        self.assertEqual(len(simulator.calls), 0)
        self.assertEqual(len(simulator.prepared), 1)

    def test_condition_building_copies_then_applies_analysis_control(self):
        experiment = example_experiment()
        parameters = {"k": 9.0, "q": 7.0}
        initial, changed = properties.build_condition(experiment, 2.0, parameters)
        self.assertEqual(initial, {"X": 3.0})
        self.assertEqual(changed, {"k": 2.0, "q": 7.0})
        initial["X"] = -99.0
        changed["q"] = -99.0
        self.assertEqual(experiment["initial_values"], {"X": 3.0})
        self.assertEqual(parameters, {"k": 9.0, "q": 7.0})
        experiment["analysis_control"] = {"type": "species_initial_concentration", "id": "S"}
        initial, changed = properties.build_condition(experiment, 0.0, parameters)
        self.assertEqual(initial, {"X": 3.0, "S": 0.0})
        self.assertEqual(changed, parameters)

    def test_request_parameter_overrides_do_not_leak_into_other_requests(self):
        simulator = LinearSimulator()
        experiment = example_experiment()
        experiment["analysis_control"] = {"type": "species_initial_concentration", "id": "S"}
        requests = requests_for(experiment, {
            "altered": {"property": "rates", "parameters": {"k": 4.0}},
            "original": {"property": "rates"},
        })
        parameters = {"k": 9.0}
        result = properties.collect_properties(requests, simulator, experiment, parameters)
        np.testing.assert_array_equal(result["altered"]["values"], np.full((2, 3), 4.0))
        np.testing.assert_array_equal(result["original"]["values"], np.full((2, 3), 9.0))
        self.assertEqual(parameters, {"k": 9.0})


class TestLimitProperties(unittest.TestCase):

    def setUp(self):
        self.records = []

    def fake_limit(self, simulator, initial, parameters, observable, options, **kwargs):
        self.records.append((initial.copy(), parameters.copy(), options.copy()))
        mode = options["convergence_mode"]
        converged = mode == "output_plateau"
        # Mutates local inputs to expose sharing between independent runs.
        initial["X"] = -99.0
        parameters["k"] = -99.0
        return {"value": 8.0 if converged else np.nan, "converged": converged,
                "status": "converged" if converged else "max_time_reached", "time_minutes": 400.0,
                "final_output": 8.0, "final_concentrations": np.asarray([8.0]), "convergence_mode": mode}

    def test_both_limit_properties_accept_species_controls_and_independent_conditions(self):
        experiment = example_experiment()
        experiment["analysis_control"] = {"type": "species_initial_concentration", "id": "S"}
        experiment["control_values"] = np.asarray([0.0, 2.0])
        requests = requests_for(experiment, {
            "whole_system": {"property": "steady_state", "initial_values": {"U": 1.0}},
            "output_only": {"property": "output_plateau", "initial_values": {"U": 3.0}},
        })
        parameters = {"k": 9.0}
        with patch.object(properties, "find_steady_state", side_effect=self.fake_limit), redirect_stdout(StringIO()):
            result = properties.collect_properties(requests, LinearSimulator(), experiment, parameters)
        self.assertEqual([record[0]["S"] for record in self.records], [0, 2, 0, 2])
        self.assertEqual([record[0]["U"] for record in self.records], [1, 1, 3, 3])
        self.assertTrue(all(record[0]["X"] == 3 and record[1]["k"] == 9 for record in self.records))
        self.assertEqual(result["whole_system"]["convergence_mode"], "full_system")
        self.assertEqual(result["output_only"]["convergence_mode"], "output_plateau")
        self.assertTrue(np.all(np.isnan(result["whole_system"]["values"])))
        np.testing.assert_array_equal(result["output_only"]["values"], [8, 8])
        np.testing.assert_array_equal(result["output_only"]["control_values"], [0, 2])
        self.assertEqual(result["whole_system"]["converged"].dtype.kind, "b")
        self.assertEqual(parameters, {"k": 9.0})
        self.assertEqual(experiment["initial_values"], {"X": 3.0})

    def test_mixed_time_and_limit_requests_start_independent_conditions(self):
        simulator = LinearSimulator()
        experiment = example_experiment()
        requests = requests_for(experiment, {
            "limit": {"property": "output_plateau", "control_values": [0.0]},
            "speed": {"property": "rates"}, "output": {"property": "concentrations"},
        })
        with patch.object(properties, "find_steady_state", side_effect=self.fake_limit), redirect_stdout(StringIO()):
            result = properties.collect_properties(requests, simulator, experiment, {"k": 9.0})
        self.assertEqual(set(result), {"limit", "speed", "output"})
        self.assertEqual(len(simulator.calls), 2)
        self.assertEqual(len(self.records), 1)
        self.assertTrue(all(call["initial_values"]["X"] == 3 for call in simulator.calls))
        self.assertEqual(self.records[0][1]["k"], 0.0)

    def test_limit_criterion_cannot_be_overridden_by_old_setting(self):
        experiment = example_experiment()
        experiment["limit_settings"] = {"convergence_mode": "output_plateau"}
        with patch.object(properties, "find_steady_state", side_effect=self.fake_limit), redirect_stdout(StringIO()):
            result = properties.collect_limit_property("steady_state", LinearSimulator(), experiment, {"k": 9.0})
        self.assertEqual(result["convergence_mode"], "full_system")
        self.assertTrue(all(record[2]["convergence_mode"] == "full_system" for record in self.records))


class TestImportedProperties(unittest.TestCase):

    def test_relative_and_absolute_imports_receive_resolved_conditions_and_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "my_property.py"
            path.write_text(
                "def calculate_property(simulator, experiment, parameters, settings):\n"
                "    return {'values': [parameters['q'] * value + experiment['initial_values']['U'] + settings['offset']\n"
                "                       for value in experiment['control_values']],\n"
                "            'units': experiment['output_units'], 'note': 'custom result'}\n", encoding="utf-8")
            for pointer in (path.name, str(path)):
                definitions = {"my_quantity": {"definition_file": pointer, "settings": {"offset": 4.0}}}
                experiment = example_experiment()
                requests = requests_for(experiment, {
                    "response": {"property": "my_quantity", "initial_values": {"U": 2.0},
                                 "parameters": {"q": 3.0}, "control_values": [0.0, 5.0, 10.0]}
                }, definitions)
                functions = properties.load_property_functions(["my_quantity"], definitions, directory)
                parameters = {"k": 9.0, "q": 1.0}
                result = properties.collect_properties(requests, LinearSimulator(), experiment, parameters, functions, definitions)
                np.testing.assert_array_equal(result["response"]["values"], [6, 21, 36])
                np.testing.assert_array_equal(result["response"]["control_values"], [0, 5, 10])
                self.assertEqual(result["response"]["note"], "custom result")
                self.assertEqual(parameters, {"k": 9.0, "q": 1.0})
                self.assertEqual(experiment["initial_values"], {"X": 3.0})

    def test_custom_mutations_do_not_change_later_requests_or_shared_inputs(self):
        experiment = example_experiment()
        parameters = {"k": 9.0}
        definitions = {"mutating_quantity": {"definition_file": "unused.py", "settings": {"offset": 4.0}}}
        requests = requests_for(experiment, {
            "custom": {"property": "mutating_quantity", "control_values": [0.0]},
            "ordinary": {"property": "rates"},
        }, definitions)

        def mutate(simulator, resolved_experiment, resolved_parameters, settings):
            resolved_experiment["initial_values"]["X"] = -99.0
            resolved_experiment["control_values"][0] = -99.0
            resolved_parameters["k"] = -99.0
            settings["offset"] = -99.0
            return {"values": [7.0], "units": "nM"}

        result = properties.collect_properties(requests, LinearSimulator(), experiment, parameters,
                                               {"mutating_quantity": mutate}, definitions)
        np.testing.assert_array_equal(result["ordinary"]["values"][:, 0], [1.0, 2.0])
        np.testing.assert_array_equal(result["custom"]["control_values"], [0.0])
        np.testing.assert_array_equal(requests["custom"]["control_values"], [0.0])
        np.testing.assert_array_equal(experiment["control_values"], [1.0, 2.0])
        self.assertEqual(experiment["initial_values"], {"X": 3.0})
        self.assertEqual(parameters, {"k": 9.0})
        self.assertEqual(definitions["mutating_quantity"]["settings"], {"offset": 4.0})

    def test_import_without_calculation_function_fails_early(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty_property.py"
            path.write_text("description = 'No calculation function'\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "calculate_property"):
                properties.load_property_functions(["custom"], {"custom": {"definition_file": path.name}}, directory)

    def test_flat_results_are_savable_without_pickle(self):
        valid = {"values": np.asarray([1.0, np.nan]), "units": "nM", "valid": [True, False], "status": ["ok", "failed"]}
        checked = properties.validate_property_result("custom", valid, 2)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "property.npz"
            np.savez(str(path), **checked)
            with np.load(str(path), allow_pickle=False) as saved:
                np.testing.assert_array_equal(saved["values"], valid["values"])
                np.testing.assert_array_equal(saved["valid"], [True, False])
                self.assertEqual(saved["units"].item(), "nM")

    def test_unsavable_or_mismatched_results_are_rejected(self):
        invalid = [
            {"values": [1, 2]}, {"values": [1], "units": "nM"},
            {"values": [1, np.inf], "units": "nM"},
            {"values": [1, 2], "units": "nM", "metadata": {"nested": True}},
            {"values": [1, 2], "units": "nM", "extra": np.asarray([object()], dtype=object)},
            {"values": [1, 2], "units": "nM", "bad name": 1},
        ]
        for result in invalid:
            with self.subTest(result=result), self.assertRaises(ValueError):
                properties.validate_property_result("custom", result, 2)


if __name__ == "__main__":
    unittest.main()
