import math
import time
from pathlib import Path

import numpy as np
import pytest

import miradac as da

DATA = Path(__file__).resolve().parents[2] / "test"
EPS = 1e-14


def ref(name):
    return str(DATA / name)


@pytest.fixture
def x(env):
    return 1.0 + da.var(0) + 2.0 * da.var(1) + 5.0 * da.var(2)


def read(name):
    d = da.NDA()
    da.read_da_from_file(ref(name), d)
    return d


# --- T2.1 opaque lists -------------------------------------------------------

def test_ndalist_behaves_like_a_list(x):
    lst = da.NDAList([x, da.exp(x)])
    lst.append(da.sin(x))
    assert len(lst) == 3
    assert da.compare_da_vectors(lst[1], da.exp(x))
    lst[0] += 1.0                              # getitem copy, iadd, setitem
    assert lst[0].con == pytest.approx(x.con + 1.0)


@pytest.mark.parametrize("as_list", [da.NDAList, list])
def test_composition_accepts_ndalist_and_plain_list(x, as_list):
    lx = as_list([x, da.exp(x), da.sinh(x)])
    lu = as_list([da.sin(x), da.cos(x), da.tan(x)])
    ly = da.NDAList([da.NDA() for _ in range(3)])
    da.compose(lx, lu, ly)
    for i in range(3):
        assert da.compare_da_with_file(ref(f"da_composition_{i}.txt"), ly[i], EPS)


def test_plain_list_is_rejected_as_output(x):
    m = da.NDAList([x] * 3)
    with pytest.raises(TypeError):
        da.compose(m, m, [da.NDA() for _ in range(3)])


# --- T2.2 algorithms (mirrors ../tpsa/python-wrapper/tests/tests.py) ---------

def test_der_and_int(x):
    y = da.log(x)
    assert da.compare_da_vectors(read("da_der.txt"), da.der(y, 1), EPS)
    assert da.compare_da_vectors(read("da_int.txt"), da.int_(y, 1), EPS)
    out = da.NDA()
    da.da_der(y, 1, out)
    assert da.compare_da_vectors(read("da_der.txt"), out, EPS)
    da.da_int(y, 1, out)
    assert da.compare_da_vectors(read("da_int.txt"), out, EPS)
    with pytest.raises(IndexError):
        da.der(y, 3)


def test_substitute_number(x):
    z = da.NDA()
    da.da_substitute_const(da.exp(x), 0, 1.0, z)
    assert da.compare_da_with_file(ref("substitute_number.txt"), z, EPS)
    da.substitute(da.exp(x), 0, 1, z)
    assert da.compare_da_with_file(ref("substitute_number.txt"), z, EPS)


def test_substitute_da_vector(x):
    z = da.NDA()
    da.substitute(da.exp(x), 0, x, z)
    assert da.compare_da_with_file(ref("substitute_da_vector.txt"), z, EPS)


def test_substitute_multiple_da_vectors(x):
    z = da.NDA()
    da.substitute(da.exp(x), [0, 1], [da.sin(x), da.cos(x)], z)
    assert da.compare_da_with_file(ref("substitute_multiple_da_vectors.txt"), z, EPS)
    with pytest.raises(ValueError):
        da.substitute(da.exp(x), [0, 0], [da.sin(x), da.cos(x)], z)
    with pytest.raises(IndexError):
        da.substitute(da.exp(x), [0, 3], [da.sin(x), da.cos(x)], z)


def test_bunch_substitution(x):
    lx = da.NDAList([x, da.exp(x), da.sinh(x)])
    ly = da.NDAList([da.NDA() for _ in range(3)])
    da.substitute(lx, [0, 1], [da.sin(x), da.cos(x)], ly)
    for i in range(3):
        assert da.compare_da_with_file(ref(f"bunch_substitution_{i}.txt"), ly[i], EPS)


def test_composition_with_numbers(x):
    m = da.NDAList([x, da.exp(x), da.sinh(x)])
    pt = [0.1, -0.2, 0.3]
    y = da.compose(m, pt)
    assert y[0] == pytest.approx(1.0 + 0.1 - 0.4 + 1.5, rel=1e-15)
    assert y[1] == pytest.approx(sum(1.2**n / math.factorial(n) for n in range(5)) * math.e)
    zc = da.compose(m, [0.1j, 0.0, 0.0])
    assert zc[0] == pytest.approx(1.0 + 0.1j)
    with pytest.raises(ValueError):
        da.compose(m, [0.1, 0.2])


def test_composition_rejects_aliased_output(x):
    m = da.NDAList([x, da.exp(x), da.sinh(x)])
    with pytest.raises(ValueError):
        da.compose(m, m, m)


def test_inv_map():
    da.init(5, 2, 1000)
    try:
        x, y = da.var(0), da.var(1)
        m = da.NDAList([2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y,
                        -0.4 * x + 1.5 * y + 0.2 * y * y])
        inv = da.NDAList([da.NDA(), da.NDA()])
        da.inv_map(m, 2, inv)
        out = da.NDAList([da.NDA(), da.NDA()])
        da.compose(m, inv, out)
        for i in range(2):
            assert (out[i] - da.var(i)).norm() < 1e-12
        with pytest.raises(ValueError):
            da.inv_map(da.NDAList([m[0] + 1.0, m[1]]), 2, inv)
        with pytest.raises(ValueError):
            da.inv_map(m, 3, inv)
        del x, y, m, inv, out
    finally:
        da.clear()


def test_devide_by_element_and_compare(x):
    y = da.exp(x)
    q = da.devide_by_element(y, y)
    assert q.con == 1.0 and q.element([1, 0, 0]) == 1.0
    assert da.compare_da_vectors(y, y.copy())
    assert not da.compare_da_vectors(y, y + 1e-6)


def test_aliases(env):
    assert da.der is da.da_der and da.int_ is da.da_int
    assert da.substitute is da.da_substitute and da.compose is da.da_composition


# --- T2.3 batch evaluation ---------------------------------------------------

def test_evaluate_map_and_call(x):
    m = da.NDAList([x, da.exp(x), da.sinh(x)])
    pts = np.random.default_rng(1).uniform(-1, 1, (1000, 3))
    out = da.evaluate_map(m, pts)
    assert out.shape == (1000, 3)
    loop = np.array([da.compose(m, p) for p in pts.tolist()])
    np.testing.assert_allclose(out, loop, rtol=1e-15, atol=0)
    assert m[1](pts[7].tolist()) == pytest.approx(out[7, 1], rel=1e-15)
    with pytest.raises(ValueError):
        da.evaluate_map(m, np.zeros((4, 2)))
    with pytest.raises(ValueError):
        x([0.1, 0.2])


def test_evaluate_map_is_20x_faster_than_a_python_loop(x):
    m = da.NDAList([x, da.exp(x), da.sinh(x)])
    pts = np.random.default_rng(2).uniform(-1, 1, (100_000, 3))

    def best(f, n):
        t = []
        for _ in range(n):
            t0 = time.perf_counter()
            f()
            t.append(time.perf_counter() - t0)
        return min(t)

    rows = pts.tolist()
    t_batch = best(lambda: da.evaluate_map(m, pts), 5)
    t_loop = best(lambda: [da.compose(m, p) for p in rows], 2)
    assert t_loop / t_batch >= 20, f"only {t_loop / t_batch:.1f}x"


# --- T2.4 exponent table -----------------------------------------------------

def test_exponents(env):
    e = da.exponents()
    assert e.shape == (da.full_length(), da.nvars()) and e.dtype == np.int32
    v = da.NDA.from_coeffs(np.arange(1.0, da.full_length() + 1))    # dense
    for i in np.random.default_rng(3).integers(0, da.full_length(), 10):
        assert list(e[i]) == v.index_element(int(i))[0]


def test_list_from_another_env_raises_env_error():
    da.init(4, 3, 400)
    stale = da.NDAList([da.var(0), da.var(1), da.var(2)])
    da.clear()
    da.init(4, 3, 400)
    try:
        m = da.NDAList([da.var(0)])
        out = da.NDAList([da.NDA()])
        with pytest.raises(da.EnvError):
            da.compose(m, stale, out)
        with pytest.raises(da.EnvError):
            da.evaluate_map(stale, np.zeros((1, 3)))
        del m, out
    finally:
        da.clear()
