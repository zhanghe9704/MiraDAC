import copy
import math
import re
from pathlib import Path

import numpy as np
import pytest

import miradac as da

DATA = Path(__file__).resolve().parents[2] / "test"


def coeff(v, *exps):
    return v.element(list(exps))


def norm_slot(text):
    return re.sub(r"V \[\d+\]", "V [#]", text)


# Output of the C++ `da::base[0].print()` at order 4 with 3 variables.
BASE0_PRINT = (
    " I          V [0]              Base  [ 2 / 35 ]\n"
    "------------------------------------------------\n"
    " 1   1.000000000000000e+00     1 0 0     1\n"
    "\n"
)


# --- T1.3 construction and inspection ---------------------------------------

def test_str_matches_cpp_print(env):
    assert norm_slot(str(da.var(0))) == norm_slot(BASE0_PRINT)


def test_construction(env):
    z = da.NDA()
    assert z.con == 0.0 and z.iszero()
    c = da.NDA(2.5)
    assert c.con == 2.5 and c.length == 1 and c.n_element == 1
    assert (c.nvars, c.order) == (3, 4)
    assert repr(da.var(0) + 1.0) == "NDA(order=4, nvars=3, nonzero=2)"


def test_from_coeffs_roundtrip(env):
    v = da.exp(1.0 + da.var(0) + 2 * da.var(1))
    a = v.coeffs()
    assert a.dtype == np.float64 and a.shape == (v.length,)
    w = da.NDA.from_coeffs(a)
    assert (w - v).iszero(0.0)
    assert np.array_equal(w.coeffs(), a)
    a[0] = 99.0                        # a copy, not a view
    assert v.con != 99.0
    with pytest.raises(ValueError):
        da.NDA.from_coeffs(np.zeros(36))


def test_copy_is_independent(env):
    v = da.var(0) + 1.0
    for w in (v.copy(), copy.copy(v), copy.deepcopy(v)):
        assert w is not v and (w - v).iszero(0.0)
        w += 1.0
        assert v.con == 1.0


def test_con_setter_resets_to_constant(env):
    v = da.var(0) + 1.0
    v.con = 3.0
    assert v.con == 3.0 and v.n_element == 1


def test_elements(env):
    v = 1.0 + da.var(0) + 2 * da.var(1)
    assert v.element(2) == 2.0             # positional index
    assert coeff(v, 0, 1, 0) == 2.0        # exponent vector
    v.set_element([1, 0, 2], 7.0)
    assert coeff(v, 1, 0, 2) == 7.0
    exps, value = v.index_element(2)
    assert exps == [0, 1, 0] and value == 2.0
    assert v.n_element == 4


def test_norms_zero_clean_reset(env):
    v = 3.0 - 4.0 * da.var(0)
    assert v.norm() == 4.0
    assert v.weighted_norm(1.0) == pytest.approx(4.0)
    v.set_element([0, 1, 0], 1e-20)
    v.clean()
    assert coeff(v, 0, 1, 0) == 0.0
    v.clean(3.5)
    assert v.n_element == 1
    assert not v.iszero() and v.iszero(5.0)
    v.reset()
    assert v.iszero(0.0)


def test_base_proxy(env):
    assert len(da.base) == 3
    assert (da.base[2] - da.var(2)).iszero(0.0)


# --- T1.4 arithmetic --------------------------------------------------------

def test_binary_da_da(env):
    x = 1.0 + da.var(0) + 2 * da.var(1)    # 1 + u
    y = 3.0 + da.var(2)
    s = x + y
    assert (s.con, coeff(s, 1, 0, 0), coeff(s, 0, 1, 0), coeff(s, 0, 0, 1)) == (4, 1, 2, 1)
    d = x - y
    assert (d.con, coeff(d, 0, 0, 1)) == (-2, -1)
    p = x * y                              # 3 + 3a + 6b + c + ac + 2bc
    assert (p.con, coeff(p, 1, 0, 0), coeff(p, 0, 1, 0), coeff(p, 0, 0, 1)) == (3, 3, 6, 1)
    assert (coeff(p, 1, 0, 1), coeff(p, 0, 1, 1)) == (1, 2)
    q = y / x                              # (3 + c)(1 - u + u^2 - ...)
    assert q.con == 3 and coeff(q, 1, 0, 0) == -3 and coeff(q, 0, 0, 1) == 1
    assert coeff(q, 2, 0, 0) == pytest.approx(3) and coeff(q, 1, 0, 1) == pytest.approx(-1)


def test_binary_with_float(env):
    x = 1.0 + da.var(0)
    assert (x + 2.0).con == 3 and (2.0 + x).con == 3
    assert (x - 2.0).con == -1 and (2.0 - x).con == 1 and coeff(2.0 - x, 1, 0, 0) == -1
    assert coeff(x * 3.0, 1, 0, 0) == 3 and coeff(3 * x, 1, 0, 0) == 3
    assert coeff(x / 4.0, 1, 0, 0) == 0.25
    r = 2.0 / x                            # 2(1 - a + a^2 - a^3 + a^4)
    assert [coeff(r, k, 0, 0) for k in range(5)] == [2, -2, 2, -2, 2]


def test_unary_and_pow(env):
    x = 1.0 + da.var(0) + 2 * da.var(1)
    n = -x
    assert (n.con, coeff(n, 0, 1, 0)) == (-1, -2)
    p = +x
    assert p is not x and (p - x).iszero(0.0)
    sq = x ** 2                            # 1 + 2u + u^2
    assert (coeff(sq, 1, 0, 0), coeff(sq, 2, 0, 0), coeff(sq, 1, 1, 0), coeff(sq, 0, 2, 0)) == (2, 1, 4, 4)
    r = x ** 0.5                           # 1 + u/2 - u^2/8
    assert coeff(r, 1, 0, 0) == pytest.approx(0.5)
    assert coeff(r, 2, 0, 0) == pytest.approx(-0.125)
    assert coeff(r, 1, 1, 0) == pytest.approx(-0.5)


def test_inplace_keeps_identity(env):
    x = 1.0 + da.var(0)
    y = da.var(1)
    i = id(x)
    x += y; x += 1.0
    assert id(x) == i and x.con == 2 and coeff(x, 0, 1, 0) == 1
    x -= y; x -= 1.0
    assert id(x) == i and x.con == 1 and coeff(x, 0, 1, 0) == 0
    x *= y
    assert id(x) == i and coeff(x, 1, 1, 0) == 1 and coeff(x, 0, 1, 0) == 1
    x *= 2.0
    assert id(x) == i and coeff(x, 1, 1, 0) == 2
    x /= 2.0
    assert id(x) == i and coeff(x, 1, 1, 0) == 1
    z = 1.0 + da.var(0)
    j = id(z)
    z /= z.copy()
    assert id(z) == j and z.con == pytest.approx(1) and z.iszero(1e-14) is False
    assert (z - 1.0).iszero(1e-14)


# --- T1.5 math functions ----------------------------------------------------

@pytest.fixture
def xfile(env):
    return 1.0 + da.var(0) + 2.0 * da.var(1) + 5.0 * da.var(2)


@pytest.mark.parametrize("name, fn", [
    ("exp_da.txt", da.exp),
    ("log_da.txt", da.log),
    ("sqrt_da.txt", da.sqrt),
    ("pow0p3_da.txt", lambda x: da.pow(x, 0.3)),
])
def test_functions_match_reference_files(xfile, name, fn):
    assert da.compare_da_with_file(str(DATA / name), fn(xfile), 1e-14)


def test_pow3_matches_reference_file(xfile):
    # compare_da_with_file cannot read pow3_da.txt: read_da_from_file reports
    # success only when the file lists terms up to full_length() - 1, and x**3
    # stops at index 19 of 35 (the C++ test ignores that return value too).
    y = da.pow(xfile, 3)
    lines = (DATA / "pow3_da.txt").read_text().splitlines()[2:]
    terms = [line.split() for line in lines if line.strip()]
    assert len(terms) == y.n_element == 20
    for t in terms:
        assert y.element([int(e) for e in t[2:-1]]) == pytest.approx(float(t[1]), rel=1e-14)


def test_function_constants(env):
    x = 0.5 + da.var(0)
    for name in ("sqrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan",
                 "sinh", "cosh", "tanh", "asinh", "atanh", "erf"):
        assert getattr(da, name)(x).con == pytest.approx(getattr(math, name)(0.5)), name
    assert da.acosh(1.5 + da.var(0)).con == pytest.approx(math.acosh(1.5))
    assert da.abs(-2.0 + da.var(0)) == pytest.approx(2.0)
    assert coeff(da.sin(x), 1, 0, 0) == pytest.approx(math.cos(0.5))
