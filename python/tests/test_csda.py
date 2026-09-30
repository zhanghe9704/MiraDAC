"""CSDA bindings (plan T6.1-T6.3): the port of test/test_symbolic_cd.cc, the
checks of examples/example_complex_symbolic.cc, and the binding behavior."""

import copy
import math

import pytest

import miradac as da

pytestmark = pytest.mark.skipif(not da.HAS_SYMBOLIC, reason="built without symbolic support")


@pytest.fixture
def env32():
    """Order 3, 2 variables, as test_symbolic_cd.cc."""
    da.init(3, 2, 800)
    yield
    da.clear()


def make_csda(c0r, b1r, b2r, c0i, b1i, b2i):
    re = da.SDA(c0r) + b1r * da.svar(0) + b2r * da.svar(1)
    im = da.SDA(c0i) + b1i * da.svar(0) + b2i * da.svar(1)
    return da.CSDA(re, im)


def make_cnda(c0r, b1r, b2r, c0i, b1i, b2i):
    re = da.NDA(c0r) + b1r * da.var(0) + b2r * da.var(1)
    im = da.NDA(c0i) + b1i * da.var(0) + b2i * da.var(1)
    return da.CNDA(re, im)


def same(z, n, eps):
    return da.compare_cd_vectors(z, n, eps)


def coeffs(v):
    return [str(e) for e in v.coeffs()]


# --- port of test_symbolic_cd.cc -------------------------------------------------

def test_structural_check(env32):
    (b,) = da.symbols("b")
    z = make_csda(1.0, b, 0.0, 0.5, 0.2, 0.0)
    assert float(z.real.con) == pytest.approx(1.0)
    assert float(z.imag.con) == pytest.approx(0.5)
    diff = z - z
    assert float(diff.real.con) == 0.0 and float(diff.imag.con) == 0.0


def test_arithmetic_cross_check(env32):
    b1, b2, c1, c2, e1, e2, f1, f2 = syms = da.symbols("b1 b2 c1 c2 e1 e2 f1 f2")
    z1 = make_csda(0.8, b1, b2, 0.3, c1, c2)
    z2 = make_csda(0.5, e1, e2, 0.6, f1, f2)
    v = [0.3, -0.1, 0.2, 0.4, -0.15, 0.05, 0.1, -0.2]
    vals = dict(zip(syms, v))
    n1 = make_cnda(0.8, v[0], v[1], 0.3, v[2], v[3])
    n2 = make_cnda(0.5, v[4], v[5], 0.6, v[6], v[7])
    c = 0.5 - 0.3j
    for rs, rn in [(z1 + z2, n1 + n2), (z1 - z2, n1 - n2), (z1 * z2, n1 * n2),
                   (z1 / z2, n1 / n2), (z1 + c, n1 + c), (z1 * 2.5, n1 * 2.5),
                   (1.0 / z2, 1.0 / n2), ((-1.0) * z1, (-1.0) * n1)]:
        assert same(da.evaluate(rs, vals), rn, 1e-10)


@pytest.mark.parametrize("fn, cons", [
    ("exp", (0.5, 0.3)), ("sqrt", (0.5, 0.3)), ("log", (0.5, 0.3)),
    ("asin", (0.3, 0.2)), ("acos", (0.3, 0.2)), ("atan", (0.3, 0.2)),
    ("asinh", (0.5, 0.3)), ("atanh", (0.5, 0.3)), ("acosh", (1.5, 0.2))])
def test_function_cross_check(env32, fn, cons):
    b1, c1 = da.symbols("b1 c1")
    z = make_csda(cons[0], b1, 0.0, cons[1], c1, 0.0)
    vals = {b1: 0.08, c1: 0.05}
    n = make_cnda(cons[0], 0.08, 0.0, cons[1], 0.05, 0.0)
    f = getattr(da, fn)
    assert same(da.evaluate(f(z), vals), f(n), 1e-9)


def test_exp_sqrt_log_two_symbols(env32):
    b1, b2, c1, c2 = da.symbols("b1 b2 c1 c2")
    z = make_csda(0.5, b1, b2, 0.3, c1, c2)
    vals = {b1: 0.10, b2: 0.05, c1: 0.08, c2: -0.04}
    n = make_cnda(0.5, 0.10, 0.05, 0.3, 0.08, -0.04)
    for f in (da.exp, da.sqrt, da.log):
        assert same(da.evaluate(f(z), vals), f(n), 1e-9)


def test_pow_cross_check(env32):
    b1, c1 = da.symbols("b1 c1")
    z = make_csda(0.5, b1, 0.0, 0.3, c1, 0.0)
    vals = {b1: 0.08, c1: 0.05}
    n = make_cnda(0.5, 0.08, 0.0, 0.3, 0.05, 0.0)
    assert same(da.evaluate(da.pow(z, 3), vals), da.pow(n, 3), 1e-9)
    assert same(da.evaluate(da.pow(z, 2.5), vals), da.pow(n, 2.5), 1e-9)
    assert same(da.evaluate(z ** 3, vals), n ** 3, 1e-9)


def test_abs_is_symbolic_magnitude(env32):
    b1, c1 = da.symbols("b1 c1")
    z = make_csda(0.6, b1, 0.0, 0.4, c1, 0.0)
    n = make_cnda(0.6, 0.1, 0.0, 0.4, 0.08, 0.0)
    r = da.abs(z)
    assert isinstance(r, da.SDA)
    mag = da.sqrt(n.real * n.real + n.imag * n.imag)
    assert da.compare_da_vectors(da.evaluate(r, {b1: 0.1, c1: 0.08}), mag, 1e-10)


def test_promote_evaluate_round_trip(env32):
    orig = da.CNDA(1.5 + 0.3 * da.var(0) - 0.1 * da.var(1),
                   0.7 + 0.2 * da.var(0) - 0.05 * da.var(1))
    p = da.promote(orig)
    assert isinstance(p, da.CSDA)
    assert same(da.evaluate(p, {}), orig, 1e-14)


def test_cd_composition_cross_check(env32):
    b1, c1 = da.symbols("b1 c1")
    z0 = make_csda(1.2, b1, 0.0, 0.4, c1, 0.0)
    z1 = make_csda(0.5, 0.0, 0.0, 0.2, 0.0, 0.0)
    v0_re, v0_im = 0.7 + 0.1 * da.var(0), 0.05 * da.var(0)
    v1_re, v1_im = 0.3 + 0.05 * da.var(1), 0.02 * da.var(1)
    cvec = da.CSDAList([da.CSDA(da.promote(v0_re), da.promote(v0_im)),
                        da.CSDA(da.promote(v1_re), da.promote(v1_im))])
    out = da.CSDAList([da.CSDA(), da.CSDA()])
    da.cd_composition([z0, z1], cvec, out)

    vals = {b1: 0.3, c1: 0.15}
    n0 = da.CNDA(1.2 + 0.3 * da.var(0), 0.4 + 0.15 * da.var(0))
    n1 = da.CNDA(da.NDA(0.5), da.NDA(0.2))
    nout = da.CNDAList([da.CNDA(), da.CNDA()])
    da.cd_composition([n0, n1], [da.CNDA(v0_re, v0_im), da.CNDA(v1_re, v1_im)], nout)
    for s, n in zip(out, nout):
        assert same(da.evaluate(s, vals), n, 1e-9)


def test_cd_composition_sda_into_complex_map(env32):
    b1, c1 = da.symbols("b1 c1")
    ivec = da.SDAList([1.1 + b1 * da.svar(0), 0.5 + c1 * da.svar(1)])
    cv0re, cv1re = 0.5 + 0.1 * da.var(0), da.NDA(0.3)
    cv0im, cv1im = da.NDA(0.2), 0.05 * da.var(1)
    cmap = [da.CSDA(da.promote(cv0re), da.promote(cv0im)),
            da.CSDA(da.promote(cv1re), da.promote(cv1im))]
    out = da.CSDAList([da.CSDA(), da.CSDA()])
    da.cd_composition(ivec, cmap, out)

    nivec = [1.1 + 0.3 * da.var(0), 0.5 + 0.2 * da.var(1)]
    nout = da.CNDAList([da.CNDA(), da.CNDA()])
    da.cd_composition(nivec, [da.CNDA(cv0re, cv0im), da.CNDA(cv1re, cv1im)], nout)
    vals = {b1: 0.3, c1: 0.2}
    for s, n in zip(out, nout):
        assert same(da.evaluate(s, vals), n, 1e-9)


# --- checks of example_complex_symbolic.cc --------------------------------------

def test_example_complex_symbolic(env32):
    a, b = da.symbols("a b")
    re = da.SDA(1.5) + a * da.svar(0) + da.SDA(0.1) * da.svar(1)
    im = da.SDA(0.5) + b * da.svar(0) + da.SDA(-0.05) * da.svar(1)
    z = da.CSDA(re, im)

    z2 = z * z
    assert float(z2.real.con) == pytest.approx(2.0)
    assert float(z2.imag.con) == pytest.approx(1.5)

    ez = da.exp(z)
    assert float(ez.real.con) == pytest.approx(math.exp(1.5) * math.cos(0.5), rel=1e-14)
    assert float(ez.imag.con) == pytest.approx(math.exp(1.5) * math.sin(0.5), rel=1e-14)

    vals = {a: 0.3, b: 0.2}
    n = da.CNDA(1.5 + 0.3 * da.var(0) + 0.1 * da.var(1),
                0.5 + 0.2 * da.var(0) - 0.05 * da.var(1))
    assert same(da.evaluate(ez, vals), da.exp(n), 1e-12)
    assert same(da.evaluate(da.sqrt(z), vals), da.sqrt(n), 1e-12)


# --- T6.1 class ----------------------------------------------------------------

def test_constructors_and_parts(env32):
    (a,) = da.symbols("a")
    re = a * da.svar(0) + 1
    z = da.CSDA(re)
    assert coeffs(z.real) == coeffs(re) and coeffs(z.imag) == ["0"]
    assert coeffs(da.CSDA().real) == ["0"]
    z.imag = 2 * re
    assert coeffs(z.imag) == coeffs(2 * re)
    r = z.real                                          # a copy
    r += 1
    assert z.real.con == 1
    assert repr(z) == "CSDA(order=3, nvars=2, nonzero=(2, 2))"
    c = z.copy()
    assert c is not z and coeffs(copy.copy(z).real) == coeffs(copy.deepcopy(z).real)
    c.real = da.SDA(0.0)
    assert z.real.con == 1
    with pytest.raises(TypeError):
        da.CSDA(da.var(0))


def test_parts_must_share_the_env(env32):
    s = da.svar(0)
    z = da.CSDA(s)
    da.clear()
    da.init(3, 2, 100)
    t = da.svar(1)
    with pytest.raises(da.EnvError):
        da.CSDA(t, s)
    with pytest.raises(da.EnvError):
        z + 1.0


# --- T6.2 arithmetic and functions -----------------------------------------------

def test_scalar_operators_match_cnda(env32):
    x = 0.4 + 0.3 * da.var(0) - 0.2 * da.var(1)
    y = 0.1 - 0.5 * da.var(0) * da.var(1)
    n = da.CNDA(x, y)
    z = da.promote(n)
    w = 0.5 - 2.0j
    for fs, fn in [(lambda u: u + w, lambda u: u + w), (lambda u: w - u, lambda u: w - u),
                   (lambda u: u * w, lambda u: u * w), (lambda u: w / u, lambda u: w / u),
                   (lambda u: u / w, lambda u: u / w), (lambda u: 2 + u, lambda u: 2 + u),
                   (lambda u: 2.0 - u, lambda u: 2.0 - u), (lambda u: u * 3, lambda u: u * 3),
                   (lambda u: 2.0 / u, lambda u: 2.0 / u), (lambda u: -u, lambda u: -u),
                   (lambda u: +u, lambda u: +u), (lambda u: u ** 0.5, lambda u: u ** 0.5)]:
        r = fs(z)
        assert isinstance(r, da.CSDA)
        assert same(da.evaluate(r, {}), fn(n), 1e-13)
    with pytest.raises(ValueError):
        z / 0.0


def test_sda_complex_gives_csda(env32):
    x = 0.4 + 0.3 * da.var(0)
    s = da.promote(x)
    w = 0.5 - 2.0j
    for r, n in [(s + w, x + w), (w + s, w + x), (s - w, x - w), (w - s, w - x),
                 (s * w, x * w), (w * s, w * x), (s / w, x / w), (w / s, w / x)]:
        assert isinstance(r, da.CSDA)
        assert same(da.evaluate(r, {}), n, 1e-13)
    assert type(s + 2) is da.SDA and type(s * 2.0) is da.SDA


@pytest.mark.parametrize("rhs", ["csda", 0.5 - 2.0j, 2.0, 3])
def test_inplace_operators(env32, rhs):
    (a,) = da.symbols("a")
    z1 = make_csda(0.8, a, 0.1, 0.3, 0.2, a)
    z2 = make_csda(0.5, 0.3, a, 0.6, 0.1, 0.2)
    other = z2 if rhs == "csda" else rhs
    for op, iop in (("__add__", "__iadd__"), ("__sub__", "__isub__"),
                    ("__mul__", "__imul__"), ("__truediv__", "__itruediv__")):
        c = z1.copy()
        before = id(c)
        expected = getattr(z1, op)(other)
        c = getattr(c, iop)(other)
        assert id(c) == before
        assert same(da.evaluate(c, {a: 0.7}), da.evaluate(expected, {a: 0.7}), 1e-13)


def test_sda_and_expr_operands(env32):
    b1, b2, c1, c2, e1, k = da.symbols("b1 b2 c1 c2 e1 k")
    z = make_csda(0.8, b1, b2, 0.3, c1, c2)
    s = da.SDA(0.5) + e1 * da.svar(0)
    vals = {b1: 0.3, b2: -0.1, c1: 0.2, c2: 0.4, e1: -0.15, k: 1.7}
    zn = make_cnda(0.8, 0.3, -0.1, 0.3, 0.2, 0.4)
    sn = da.NDA(0.5) - 0.15 * da.var(0)
    operands = [(s, sn), (k, 1.7)]
    try:
        import symengine
        operands.append((symengine.Symbol("k"), 1.7))
    except ImportError:
        pass
    for op in ("__add__", "__sub__", "__mul__", "__truediv__"):
        rop = op.replace("__", "__r", 1)
        for sym, num in operands:
            got = getattr(z, op)(sym)
            assert isinstance(got, da.CSDA)
            assert same(da.evaluate(got, vals), getattr(zn, op)(num), 1e-12)
            got = getattr(z, rop)(sym)
            assert same(da.evaluate(got, vals), getattr(zn, rop)(num), 1e-12)


def test_str_prints_both_parts(env32):
    (b,) = da.symbols("b")
    text = str(make_csda(1.0, b, 0.0, 0.5, 0.2, 0.0))
    assert "Real part" in text and "Imaginary part" in text and "b" in text


def test_unsupported_operands(env32):
    z = da.CSDA(da.svar(0))
    for other in (da.var(0), da.CNDA(da.var(0))):
        with pytest.raises(TypeError):
            z + other
        with pytest.raises(TypeError):
            other * z


def test_functions_without_complex_form(env32):
    z = da.CSDA(da.SDA(0.5) + da.svar(0))
    for name in ("sin", "cos", "tan", "sinh", "cosh", "tanh", "erf"):
        with pytest.raises(TypeError):
            getattr(da, name)(z)
    assert isinstance(da.sin(da.svar(0)), da.SDA)        # SDA overloads kept
    assert isinstance(da.exp(da.CNDA(da.var(0))), da.CNDA)


def test_symengine_keys(env32):
    se = pytest.importorskip("symengine")
    a = se.Symbol("a")
    ea = da.Expr.from_symengine(a)
    z = make_csda(0.5, ea, 0.0, 0.3, 0.0, ea)
    n = make_cnda(0.5, 0.25, 0.0, 0.3, 0.0, 0.25)
    assert same(da.evaluate(z, {a: 0.25}), n, 1e-14)
    assert same(da.evaluate(z, [a], [0.25]), n, 1e-14)


# --- T6.3 algorithms -------------------------------------------------------------

def test_evaluate_forms_and_errors(env32):
    a, b = da.symbols("a b")
    z = make_csda(0.5, a, 0.0, 0.3, 0.0, b)
    n = make_cnda(0.5, 0.1, 0.0, 0.3, 0.0, 0.2)
    assert isinstance(da.evaluate(z, {a: 0.1, b: 0.2}), da.CNDA)
    assert same(da.evaluate(z, [a, b], [0.1, 0.2]), n, 1e-15)
    with pytest.raises(ValueError):
        da.evaluate(z, {a: 0.1})
    with pytest.raises(ValueError):
        da.evaluate(z, [a, b], [0.1])


def test_cd_composition_checks(env32):
    z = [da.CSDA(da.svar(0)), da.CSDA(da.svar(1))]
    v = da.CSDAList([da.CSDA(da.svar(1)), da.CSDA(da.svar(0))])
    out = da.CSDAList([da.CSDA(), da.CSDA()])
    with pytest.raises(ValueError):
        da.cd_composition(z, v[:1], out)
    with pytest.raises(ValueError):
        da.cd_composition(z, v, out[:1])
    with pytest.raises(TypeError):
        da.cd_composition(z, v, [da.CSDA(), da.CSDA()])       # output must be a CSDAList
    with pytest.raises(ValueError):
        da.cd_composition(da.SDAList([da.svar(0), da.svar(1)]), v, v)
    da.cd_composition(z, v, out)                             # swaps the variables
    assert same(da.evaluate(out[0], {}), da.CNDA(da.var(1)), 0)
    assert same(da.evaluate(out[1], {}), da.CNDA(da.var(0)), 0)


def test_csdalist_behaves_like_a_list(env32):
    lst = da.CSDAList([da.CSDA(da.svar(0))])
    lst.append(da.CSDA())
    assert len(lst) == 2 and isinstance(lst[0], da.CSDA)
