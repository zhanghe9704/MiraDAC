"""Every script in python/examples runs to completion (plan T8.2)."""
import pathlib
import subprocess
import sys

import pytest

EXAMPLES = sorted((pathlib.Path(__file__).parents[1] / "examples").glob("*.py"))


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.stem)
def test_example_runs(path, tmp_path):
    # In a subprocess, since each example calls init() and clear(); in
    # tmp_path, since examples.py writes da_output.txt.
    r = subprocess.run([sys.executable, str(path)], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "FAIL" not in r.stdout


def test_all_examples_present():
    names = {p.stem for p in EXAMPLES}
    assert names == {"examples", "example_interop", "example_complex_da",
                     "example_1_symbolic", "example_complex_symbolic"}
