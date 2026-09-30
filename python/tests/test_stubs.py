"""Stubs, typing and docstrings (plan T8.1, T8.3)."""
import os
import pathlib
import subprocess
import sys

import pytest

import miradac._core as core

ROOT = pathlib.Path(__file__).parents[2]
STUB = ROOT / "python" / "miradac" / "_core.pyi"
PATTERNS = ROOT / "python" / "stubgen_patterns.txt"


@pytest.mark.skipif(not core.HAS_SYMBOLIC, reason="the committed stub describes the symbolic build")
def test_stub_is_current(tmp_path):
    """The committed stub equals a fresh nanobind.stubgen run (see PATTERNS)."""
    pytest.importorskip("nanobind.stubgen")
    out = tmp_path / "_core.pyi"
    subprocess.run([sys.executable, "-m", "nanobind.stubgen", "-q", "-m", "miradac._core", "-P",
                    "-p", str(PATTERNS), "-o", str(out)],
                   check=True, capture_output=True, cwd=ROOT)
    assert out.read_text() == STUB.read_text(), "regenerate the stub, see " + PATTERNS.name


def test_py_typed():
    assert (STUB.parent / "py.typed").exists()


def test_mypy_strict():
    """mypy --strict passes on the examples and on the package itself."""
    pytest.importorskip("mypy")
    # MYPYPATH makes mypy check the package as source, not as an installed
    # (and therefore silenced) one.
    env = dict(os.environ, MYPYPATH=str(ROOT / "python"))
    for target in (["-p", "miradac"], ["python/examples"]):
        r = subprocess.run([sys.executable, "-m", "mypy", "--strict", *target],
                           cwd=ROOT, env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


def _has_doc(f):
    sigs = getattr(f, "__nb_signature__", None)
    if sigs is None:
        return bool(getattr(f, "__doc__", None))
    return any(doc for _, doc, *_ in sigs)


def test_every_binding_has_a_docstring():
    # The *List classes' methods come from nanobind's bind_vector.
    missing = []
    for name, obj in vars(core).items():
        if name.startswith("__"):
            continue
        if isinstance(obj, type):
            if not obj.__doc__:
                missing.append(name)
            if name.endswith("List") or issubclass(obj, BaseException):
                continue
            for key, member in vars(obj).items():
                if isinstance(member, property):
                    ok = bool(member.__doc__)
                elif isinstance(member, staticmethod):
                    ok = _has_doc(member.__func__)
                elif type(member).__name__ in ("nb_func", "nb_method"):
                    ok = _has_doc(member)
                else:
                    continue
                if not ok:
                    missing.append(f"{name}.{key}")
        elif type(obj).__name__ == "nb_func" and not _has_doc(obj):
            missing.append(name)
    assert missing == []
