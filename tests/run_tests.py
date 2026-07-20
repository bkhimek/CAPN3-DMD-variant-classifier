#!/usr/bin/env python3
"""Standalone pytest-style test runner.

Milestone 1 was built in a sandbox with no PyPI access, so pytest itself
could not be installed to verify the test suite. This script discovers
tests/unit/test_*.py modules, calls every test_* function in them, and
reports pass/fail — enough to self-verify the suite without the pytest
package. If pytest IS available in your environment (it will be, once you
`pip install -r requirements.txt` on your own machine), just run `pytest`
instead; the test files need no changes either way.

Usage: python3 tests/run_tests.py   (run from the repository root)
"""

import importlib
import sys
import traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
UNIT_TEST_DIR = REPO_ROOT / "tests" / "unit"


def discover_test_modules():
    sys.path.insert(0, str(SRC_DIR))
    sys.path.insert(0, str(UNIT_TEST_DIR))
    modules = []
    for path in sorted(UNIT_TEST_DIR.glob("test_*.py")):
        module_name = path.stem
        modules.append(importlib.import_module(module_name))
    return modules


def run() -> int:
    modules = discover_test_modules()
    passed = 0
    failed = []

    for module in modules:
        test_functions = [
            (name, getattr(module, name))
            for name in dir(module)
            if name.startswith("test_") and callable(getattr(module, name))
        ]
        for name, func in test_functions:
            full_name = f"{module.__name__}.{name}"
            try:
                func()
            except Exception as exc:  # noqa: BLE001 - a failing test can raise anything
                failed.append((full_name, exc, traceback.format_exc()))
                print(f"FAIL  {full_name}: {exc}")
            else:
                passed += 1
                print(f"PASS  {full_name}")

    total = passed + len(failed)
    print()
    print(f"{passed}/{total} passed")
    if failed:
        print()
        print("Failures:")
        for full_name, exc, tb in failed:
            print(f"--- {full_name} ---")
            print(tb)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(run())
