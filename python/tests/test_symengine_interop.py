import glob
import os
import subprocess
import sys
import time

import pytest

import miradac as da

pytestmark = pytest.mark.skipif(not da.HAS_SYMBOLIC, reason="built without symbolic support")
symengine = pytest.importorskip("symengine")


def test_status_shared():
    st = da.symengine_interop_status()
    assert st["mode"] == "shared", st["reason"]
    assert st["reason"] is None
    assert st["symengine_version"] == "0.14.1"
    assert "libsymengine.so.0.14" in st["wrapper_needed"]
    assert st["loaded_libsymengine"] == [st["expected_libsymengine"]]
    assert st["layout_selftest"] is True
    assert da.symengine_interop_status() == st


def test_round_trip_keeps_address():
    x = symengine.Symbol("x")
    e = da.Expr.from_symengine(x)
    assert da._rcp_address(e) == da._rcp_address(x)
    back = e.to_symengine()
    assert isinstance(back, symengine.Symbol)
    assert da._rcp_address(back) == da._rcp_address(x)
    assert back == x


def test_conversions_equal():
    a, b = da.symbols("a b")
    e = da.Expr("sin(a)*b**2 + 1/3")
    s = e.to_symengine()
    assert s == symengine.sympify("sin(a)*b**2 + 1/3")
    assert da.Expr.from_symengine(s * 2) == 2 * e
    with pytest.raises(TypeError):
        da.Expr.from_symengine(1.0)
    with pytest.raises(TypeError):
        da._rcp_address(1.0)


def test_refcounts_balanced():
    a, b, c = da.symbols("a b c")
    e = ((a + b + c) ** 20).expand()
    addr, count = da._rcp_address(e), da._rcp_use_count(e)
    for _ in range(1000):
        e2 = da.Expr.from_symengine(e.to_symengine())
        assert da._rcp_address(e2) == addr
    del e2
    assert da._rcp_address(e) == addr
    assert da._rcp_use_count(e) == count


def _best(f, n=5, reps=20):
    best = float("inf")
    for _ in range(n):
        t = time.perf_counter()
        for _ in range(reps):
            f()
        best = min(best, time.perf_counter() - t)
    return best / reps


def test_shared_much_faster_than_string():
    a, b, c = da.symbols("a b c")
    e = ((a + b + c) ** 20).expand()
    shared = _best(lambda: da.Expr.from_symengine(e.to_symengine()))
    # What string mode does: str() one way, parse the other.
    string = _best(lambda: da.Expr(str(symengine.sympify(str(e)))), reps=2)
    assert string / shared >= 100, (shared, string)


@pytest.mark.skipif(not glob.glob(os.path.expanduser("~/.local/lib/libsymengine.so.0.15*")),
                    reason="needs another libsymengine to swap in")
def test_import_rejects_other_libsymengine():
    other = sorted(glob.glob(os.path.expanduser("~/.local/lib/libsymengine.so.0.15*")))[-1]
    r = subprocess.run([sys.executable, "-c", "import miradac"], capture_output=True, text=True,
                       env={**os.environ, "LD_PRELOAD": other})
    assert r.returncode != 0
    assert "ImportError: miradac: loaded libsymengine" in r.stderr
