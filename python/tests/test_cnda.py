import cmath
import copy
import math
from pathlib import Path

import numpy as np
import pytest

import miradac as da

DATA = Path(__file__).resolve().parents[2] / "test"
EPS = 1e-14


def ref(name):
    return str(DATA / name)


def dense(v):
    out = np.zeros(da.full_length())
    c = v.coeffs()
    out[: len(c)] = c
    return out


def assert_nda(v, expected, tol=1e-13):
    np.testing.assert_allclose(dense(v), dense(expected), rtol=tol, atol=tol)


def assert_cnda(c, re, im, tol=1e-13):
    assert isinstance(c, da.CNDA)
    assert_nda(c.real, re, tol)
    assert_nda(c.imag, im, tol)


@pytest.fixture
def xy(env):
    """The x1, x2 of test_numeric.cc's CD test, plus a second pair."""
    t = da.var(0) + 2.0 * da.var(1) + 3.0 * da.var(2)
    s = 0.5 * da.var(0) + 4.0 * da.var(1) + 2.7 * da.var(2)
    return da.cos(t), da.sin(t), da.cos(s), da.sin(s)


@pytest.fixture
def y(xy):
    x1, x2, x3, x4 = xy
    return da.CNDA(x1, x2), da.CNDA(x3, x4)


# --- T3.1 class ---------------------------------------------------------------

def test_constructors(xy):
    x1, x2, _, _ = xy
    zero = da.CNDA()
    assert zero.real.iszero() and zero.imag.iszero()
    c = da.CNDA(x1, x2)
    assert_cnda(c, x1, x2, 0)
    assert_cnda(da.CNDA(x1), x1, da.NDA(), 0)
    assert_cnda(da.CNDA(1.5 - 2j), da.NDA(1.5), da.NDA(-2.0), 0)
    assert_cnda(da.CNDA(3.0), da.NDA(3.0), da.NDA(), 0)


def test_real_imag_are_copies_and_assignable(y, xy):
    c, _ = y
    x1, x2, x3, _ = xy
    r = c.real
    r += 1.0
    assert_nda(c.real, x1, 0)                  # getter returned a copy
    c.imag = x3
    assert_cnda(c, x1, x3, 0)
    c.real = x2
    assert_cnda(c, x2, x3, 0)


def test_copy_is_independent(y):
    c, _ = y
    for d in (c.copy(), copy.copy(c), copy.deepcopy(c)):
        d += 1.0
        assert d.real.con == pytest.approx(c.real.con + 1.0)


def test_repr_and_str(y):
    c, _ = y
    assert repr(c) == f"CNDA(order=4, nvars=3, nonzero=({c.real.n_element}, {c.imag.n_element}))"
    s = str(c)
    lines = s.splitlines()
    assert "V [" in lines[0] and "Base" in lines[0]
    union = np.count_nonzero((dense(c.real) != 0) | (dense(c.imag) != 0))
    assert len([l for l in lines[2:] if l.strip()]) == union
    assert not hasattr(c, "conj")              # no conj in da.h


def test_parts_must_share_the_env(env):
    old = da.var(0)
    da.clear()
    da.init(4, 3, 400)
    with pytest.raises(da.EnvError):
        da.CNDA(da.var(0), old)
    c = da.CNDA(da.var(0))
    with pytest.raises(da.EnvError):
        c.imag = old


# --- T3.2 arithmetic ------------------------------------------------------------

def test_cnda_cnda_operators(y):
    a, b = y
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    assert_cnda(a + b, ar + br, ai + bi)
    assert_cnda(a - b, ar - br, ai - bi)
    assert_cnda(a * b, ar * br - ai * bi, ar * bi + ai * br)
    n = br * br + bi * bi
    assert_cnda(a / b, (ar * br + ai * bi) / n, (ai * br - ar * bi) / n)


def test_cnda_scalar_operators(y):
    a, _ = y
    ar, ai = a.real, a.imag
    z = 0.5 - 2.0j
    assert_cnda(a + z, ar + z.real, ai + z.imag)
    assert_cnda(z + a, ar + z.real, ai + z.imag)
    assert_cnda(a - z, ar - z.real, ai - z.imag)
    assert_cnda(z - a, z.real - ar, z.imag - ai)
    assert_cnda(a * z, ar * z.real - ai * z.imag, ar * z.imag + ai * z.real)
    assert_cnda(z * a, ar * z.real - ai * z.imag, ar * z.imag + ai * z.real)
    w = 1 / z
    assert_cnda(a / z, ar * w.real - ai * w.imag, ar * w.imag + ai * w.real)
    n = ar * ar + ai * ai
    q = z * da.CNDA(ar / n, -1.0 * ai / n)
    assert_cnda(z / a, q.real, q.imag)

    assert_cnda(a + 2.0, ar + 2.0, ai)
    assert_cnda(2.0 + a, ar + 2.0, ai)
    assert_cnda(a - 2.0, ar - 2.0, ai)
    assert_cnda(2.0 - a, 2.0 - ar, -1.0 * ai)
    assert_cnda(a * 2.0, ar * 2.0, ai * 2.0)
    assert_cnda(2 * a, ar * 2.0, ai * 2.0)
    assert_cnda(a / 2.0, ar / 2.0, ai / 2.0)
    assert_cnda(2.0 / a, 2.0 * ar / n, -2.0 * ai / n)
    assert_cnda(-a, -ar, -ai, 0)
    assert_cnda(+a, ar, ai, 0)
    with pytest.raises(ValueError):
        a / 0.0


def test_nda_complex_operators(xy):
    x, _, _, _ = xy
    z = 0.5 - 2.0j
    assert_cnda(x + z, x + z.real, da.NDA(z.imag))
    assert_cnda(z + x, x + z.real, da.NDA(z.imag))
    assert_cnda(x - z, x - z.real, da.NDA(z.imag))       # as da.h defines it
    assert_cnda(z - x, z.real - x, da.NDA(z.imag))
    assert_cnda(x * z, x * z.real, x * z.imag)
    assert_cnda(z * x, x * z.real, x * z.imag)
    w = 1 / z
    assert_cnda(x / z, x * w.real, x * w.imag)
    assert_cnda(z / x, z.real / x, z.imag / x)
    assert type(x + 2) is da.NDA and type(x * 2.0) is da.NDA  # float overload first


@pytest.mark.parametrize("op", ["__add__", "__radd__", "__sub__", "__rsub__", "__mul__", "__rmul__",
                                "__truediv__", "__rtruediv__", "__iadd__", "__isub__", "__imul__",
                                "__itruediv__"])
def test_float_overload_before_complex(op):
    # An int fails every overload before conversion; the first one that
    # converts it must be the float one (A.5 rule 1).
    doc = getattr(da.CNDA, op).__doc__
    assert doc.index("arg: float") < doc.index("arg: complex")


def test_cnda_nda_operators(y, xy):
    a, _ = y
    x = xy[0]
    for op in ("__add__", "__sub__", "__mul__", "__truediv__"):
        rop = op.replace("__", "__r", 1)
        xc = da.CNDA(x)                                  # x + 0i
        for got, want in ((getattr(a, op)(x), getattr(a, op)(xc)),
                          (getattr(a, rop)(x), getattr(xc, op)(a))):
            assert isinstance(got, da.CNDA)
            d = got - want                               # lengths may differ
            assert d.real.norm() < 1e-12 and d.imag.norm() < 1e-12


@pytest.mark.parametrize("rhs", ["cnda", 0.5 - 2.0j, 2.0])
def test_inplace_operators(y, rhs):
    a, b = y
    other = b if rhs == "cnda" else rhs
    for op, iop in (("__add__", "__iadd__"), ("__sub__", "__isub__"),
                    ("__mul__", "__imul__"), ("__truediv__", "__itruediv__")):
        c = a.copy()
        before = id(c)
        expected = getattr(a, op)(other)
        c = getattr(c, iop)(other)
        assert id(c) == before
        assert_cnda(c, expected.real, expected.imag, 0)


def test_pow_operator(y):
    a, _ = y
    assert_cnda(a ** 3, (a * a * a).real, (a * a * a).imag)
    s = a ** 0.5
    assert_cnda(s * s, a.real, a.imag, 1e-12)


def test_env_mismatch_raises(y):
    a, _ = y
    da.clear()
    da.init(4, 3, 400)
    with pytest.raises(da.EnvError):
        a + 1.0
    with pytest.raises(da.EnvError):
        da.CNDA(da.var(0)) + a


# --- T3.3 functions -------------------------------------------------------------

def test_functions_keep_nda_first(xy):
    x = xy[0]
    assert type(da.exp(x)) is da.NDA
    assert type(da.abs(x)) is float
    assert type(da.pow(x, 2)) is da.NDA


def test_complex_functions(y):
    a, _ = y
    ar, ai = a.real, a.imag
    e = da.exp(a)
    assert_cnda(e, da.exp(ar) * da.cos(ai), da.exp(ar) * da.sin(ai))
    s = da.sqrt(a)
    assert_cnda(s * s, ar, ai)
    lg = da.log(a)
    assert_cnda(da.exp(lg), ar, ai)
    assert_cnda(da.pow(a, 3), (a * a * a).real, (a * a * a).imag)
    assert_cnda(da.pow(a, 0.5), s.real, s.imag)
    assert da.abs(a) == max(ar.norm(), ai.norm())
    b = 0.3 * a + (0.1 + 0.2j)                 # away from asin's branch point at 1
    t = da.asin(b) + da.acos(b)
    assert_cnda(t, da.NDA(math.pi / 2), da.NDA(), 1e-12)
    with pytest.raises(TypeError):
        da.sin(a)                              # no complex sin in C++


@pytest.mark.parametrize("name,con", [("atan", 0.3), ("asinh", 0.3), ("atanh", 0.3),
                                      ("acosh", 1.7), ("asin", 0.3), ("acos", 0.3)])
def test_complex_functions_on_real_input(env, name, con):
    x = con + da.var(0) + 0.5 * da.var(1)
    f = getattr(da, name)
    r = f(da.CNDA(x))
    assert_cnda(r, f(x), da.NDA(), 1e-12)
    assert r.real.con == pytest.approx(getattr(cmath, name)(con).real)


# --- T3.4 algorithms ------------------------------------------------------------

@pytest.mark.parametrize("i,op", list(enumerate(["__add__", "__sub__", "__mul__", "__truediv__"])))
def test_cd_calculation_files(y, i, op):
    y1, y2 = y
    r = getattr(y1, op)(y2)
    assert da.compare_cd_with_file(ref(f"cd_calculation_{i}.txt"), r, EPS)


def test_read_and_compare_cd(y):
    y1, y2 = y
    r = da.CNDA()
    assert da.read_cd_from_file(ref("cd_calculation_2.txt"), r)
    assert da.compare_cd_vectors(r, y1 * y2, EPS)
    assert not da.compare_cd_vectors(r, y1 + y2, EPS)


@pytest.fixture
def maps(xy, y):
    x1, x2, _, _ = xy
    y1, y2 = y
    mmap = da.NDAList([x1, x2])
    cnmap = da.CNDAList([y1, y2, y1 * y2])
    cmmap = da.CNDAList([da.CNDA(x1, da.exp(x1)), da.CNDA(x2, da.exp(x2))])
    return mmap, cnmap, cmmap


def outputs(n=2):
    return da.CNDAList([da.CNDA() for _ in range(n)])


@pytest.mark.parametrize("as_list", [False, True])
def test_cd_composition_files(maps, xy, as_list):
    mmap, cnmap, cmmap = maps
    wrap = list if as_list else (lambda l: l)
    out = outputs()
    da.cd_composition(wrap(mmap), wrap(cnmap), out)
    for i in range(2):
        assert da.compare_cd_with_file(ref(f"da_composition_cd_{i}.txt"), out[i], EPS)
    da.cd_composition(wrap(cmmap), wrap(cnmap), out)
    for i in range(2):
        assert da.compare_cd_with_file(ref(f"cd_composition_cd_{i}.txt"), out[i], EPS)
    x1, x2 = mmap
    mmap.append(x1 + 0.33 * x2)
    da.cd_composition(wrap(cmmap), wrap(mmap), out)
    for i in range(2):
        assert da.compare_cd_with_file(ref(f"cd_composition_da_{i}.txt"), out[i], EPS)


def test_cd_composition_checks(maps):
    mmap, cnmap, cmmap = maps
    with pytest.raises(TypeError):
        da.cd_composition(mmap, cnmap, [da.CNDA(), da.CNDA()])   # output list must be a CNDAList
    with pytest.raises(ValueError):
        da.cd_composition(mmap, da.CNDAList(cnmap[:2]), outputs())
    with pytest.raises(ValueError):
        da.cd_composition(mmap, da.CNDAList(), outputs())
    with pytest.raises(ValueError):
        da.cd_composition(mmap, cnmap, outputs(3))
    with pytest.raises(ValueError):
        da.cd_composition(da.NDAList(list(mmap) + [mmap[0]]), cnmap, cnmap)   # output is also v
    out = da.CNDAList(cmmap)
    da.cd_composition(cmmap, cnmap, out)       # the CNDA-map forms allow it
    ref_out = outputs()
    da.cd_composition(cmmap, cnmap, ref_out)
    for a, b in zip(out, ref_out):
        assert da.compare_cd_vectors(a, b, 0)


def test_da_composition_at_complex_points(maps):
    mmap, _, _ = maps
    pt = [0.1 + 0.2j, -0.3 + 0.05j, 0.2 - 0.1j]
    vals = da.da_composition(mmap, pt)
    out = outputs()
    da.cd_composition(mmap, da.CNDAList([da.CNDA(z) for z in pt]), out)
    for v, c in zip(vals, out):
        assert v == pytest.approx(complex(c.real.con, c.imag.con), rel=1e-14)


def test_cndalist_behaves_like_a_list(y):
    y1, y2 = y
    lst = da.CNDAList([y1])
    lst.append(y2)
    assert len(lst) == 2
    lst[0] += 1.0
    assert lst[0].real.con == pytest.approx(y1.real.con + 1.0)
