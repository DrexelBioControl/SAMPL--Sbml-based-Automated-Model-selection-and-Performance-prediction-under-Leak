"""Checks shared generated equations without importing the SBML packages."""

import ast
import hashlib
import importlib.util
import os
from pathlib import Path
import py_compile
import sys
import tempfile
from types import SimpleNamespace
import unittest


ADAPTER_FILE = Path(__file__).resolve().parents[1] / "sbmltoodepy_solveivp_adapter.py"


def read_cache_methods(converter):
    # Reads just the two file-handling methods; no model simulation is needed.
    tree = ast.parse(ADAPTER_FILE.read_text(encoding="utf-8"))
    adapter = next(node for node in tree.body if isinstance(node, ast.ClassDef))
    methods = [node for node in adapter.body if isinstance(node, ast.FunctionDef)
               and node.name in ("_generate_when_needed", "_import_generated_class")]
    module = ast.Module(body=methods)
    module.type_ignores = []
    namespace = {"hashlib": hashlib, "sbmltoodepy": converter,
                 "importlib": importlib, "sys": sys}
    exec(compile(module, str(ADAPTER_FILE), "exec"), namespace)
    return namespace["_generate_when_needed"], namespace["_import_generated_class"]


class GeneratedModelCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.generated = self.folder / "generated_models" / "model_generated.py"
        self.generated.parent.mkdir()
        self.calls = []

        def convert(source, outputFilePath, className):
            self.calls.append(source)
            value = Path(source).read_text(encoding="utf-8")
            Path(outputFilePath).write_text(
                f"class {className}:\n    value = {value!r}\n", encoding="utf-8")

        self.generate, self.load = read_cache_methods(
            SimpleNamespace(ParseAndCreateModel=convert))
        self.current = self.folder / "current" / "model.xml"
        self.current.parent.mkdir()
        self.current.write_text("A", encoding="utf-8")
        self.model = SimpleNamespace(sbml_file=self.current,
                                     generated_file=self.generated,
                                     class_name="Generated_model")

    def test_unchanged_contents_reuse_equations_even_if_timestamp_changes(self):
        self.generate(self.model)
        original = self.generated.read_bytes()
        original_mtime = self.generated.stat().st_mtime_ns
        newer = self.generated.stat().st_mtime + 100
        os.utime(str(self.current), (newer, newer))
        self.generate(self.model)
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(self.generated.read_bytes(), original)
        self.assertEqual(self.generated.stat().st_mtime_ns, original_mtime)

    def test_old_different_sbml_with_same_filename_replaces_cached_equations(self):
        self.generate(self.model)
        self.assertEqual(self.load(self.model).value, "A")
        archived = self.folder / "archived" / "model.xml"
        archived.parent.mkdir()
        archived.write_text("B", encoding="utf-8")
        older = self.generated.stat().st_mtime - 100
        os.utime(str(archived), (older, older))
        self.model.sbml_file = archived
        self.generate(self.model)
        self.assertEqual(self.load(self.model).value, "B")
        self.model.sbml_file = self.current
        self.generate(self.model)
        self.assertEqual(self.load(self.model).value, "A")
        self.assertEqual(len(self.calls), 3)

    def test_regeneration_removes_compiled_copies_with_same_timestamp_and_size(self):
        self.generate(self.model)
        original_time = self.generated.stat().st_mtime
        original_size = self.generated.stat().st_size
        compiled = Path(importlib.util.cache_from_source(str(self.generated)))
        py_compile.compile(str(self.generated), cfile=str(compiled), doraise=True)
        other_compiled = compiled.parent / "model_generated.cpython-37.opt-1.pyc"
        other_compiled.write_bytes(b"old compiled equations")
        legacy_compiled = self.generated.with_suffix(".pyc")
        legacy_compiled.write_bytes(b"old compiled equations")
        unrelated = compiled.parent / "other_model.cpython-37.pyc"
        unrelated.write_bytes(b"keep unrelated model")
        self.current.write_text("B", encoding="utf-8")
        self.generate(self.model)
        os.utime(str(self.generated), (original_time, original_time))
        self.assertEqual(self.generated.stat().st_size, original_size)
        self.assertFalse(compiled.exists())
        self.assertFalse(other_compiled.exists())
        self.assertFalse(legacy_compiled.exists())
        self.assertTrue(unrelated.exists())
        self.assertEqual(self.load(self.model).value, "B")

    def test_old_generated_file_without_hash_is_rebuilt_once(self):
        self.generated.write_text("# old generated model\n", encoding="utf-8")
        self.generate(self.model)
        expected_hash = hashlib.sha256(self.current.read_bytes()).hexdigest()
        self.assertTrue(self.generated.read_text(encoding="utf-8").startswith(
            "# SBML SHA256: " + expected_hash + "\n"))
        self.generate(self.model)
        self.assertEqual(len(self.calls), 1)


if __name__ == "__main__":
    unittest.main()
