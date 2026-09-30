"""SDA bindings (plan T5.1-T5.6)."""

import copy
import math

import pytest

import miradac as da

pytestmark = pytest.mark.skipif(not da.HAS_SYMBOLIC, reason="built without symbolic support")


def coeffs(v):
    return [str(e) for e in v.coeffs()]


@pytest.fixture
def ab(env):
    return da.symbols("a b")


# T5.1

def test_constructors(ab):
    a, _ = ab
    assert coeffs(da.SDA()) == ["0"]
    assert da.SDA(1.5).con == 1.5
    two = da.SDA(2)
    assert repr(two.con) == "Expr('2')"            # exact integer, not 2.0
    assert da.SDA(a).con == a
    with pytest.raises(TypeError):
        da.SDA("a")


def test_promote_and_svar(env):
    x = 1.5 + 2 * da.var(0) + 0.25 * da.var(2)
    s = da.promote(x)
    assert isinstance(s, da.SDA)
    assert coeffs(s) == ["1.5", "2", "0", "0.25"]
    assert coeffs(da.svar(1)) == ["0", "0", "1"]
    with pytest.raises(IndexError):
        da.svar(3)


def test_inspection(ab):
    a, b = ab
    s = a * da.svar(0) + b * da.svar(1) * da.svar(1) + 1
    assert s.con == 1
    assert (s.nvars, s.order) == (3, 4)
    assert s.length == len(s.coeffs())
    assert s.n_element == 3
    assert s.element(1) == a
    assert s.element([0, 2, 0]) == b
    assert s.index_element(1) == ([1, 0, 0], a)
    s.set_element([0, 0, 1], b / 2)
    assert s.element([0, 0, 1]) == b / 2
    s.set_element([1, 1, 0], 3)
    assert s.element([1, 1, 0]) == 3
    with pytest.raises(ValueError):
        s.set_element([-1, 0, 0], 1)
    assert repr(s) == "SDA(order=4, nvars=3, nonzero=5)"
    assert "a" in str(s) and "b" in str(s)
    s.con = a                                # as NDA: resets to a constant
    assert coeffs(s) == ["a"]


def test_iszero_clean_reset_copy(ab):
    a, _ = ab
    s = a * da.svar(0) - a * da.svar(0)
    assert s.iszero()
    t = a * da.svar(0)
    assert not t.iszero()
    c = t.copy()
    assert c is not t and c.coeffs() == t.coeffs()
    assert copy.copy(t).coeffs() == t.coeffs() == copy.deepcopy(t).coeffs()
    t.reset()
    assert t.iszero()
    assert c.element(1) == a
    u = a * da.svar(2) - a * da.svar(2) + 1
    u.clean()
    assert u.length == 1


# T5.2

def test_binary_operators(ab):
    a, _ = ab
    s = a + da.svar(0)                       # [a, 1]
    t = da.SDA(2) * da.svar(0)               # [0, 2]
    x = da.var(0)                            # NDA [0, 1]
    for other, want in [
        (t, "3 + a"),                        # SDA
        (0.5, None), (3, None), (a, None),   # float, int, Expr
        (x, "2 + a"),                        # NDA
    ]:
        r = s + other
        assert isinstance(r, da.SDA)
        if want:
            assert str(r.element(1) + r.con) == want
    assert coeffs(s + 3) == ["3 + a", "1"]
    assert coeffs(3 + s) == ["3 + a", "1"]
    assert coeffs(s - 1) == ["-1 + a", "1"]
    assert coeffs(1 - s) == ["1 - a", "-1"]
    assert coeffs(s * 2) == ["2*a", "2"]
    assert coeffs(2 * s) == ["2*a", "2"]
    assert coeffs(s / 2) == ["(1/2)*a", "1/2"]
    assert coeffs(s + 0.5) == ["0.5 + a", "1"]
    assert coeffs(s * 0.5) == ["0.5*a", "0.5"]
    assert coeffs(s * a) == ["a**2", "a"]
    assert coeffs(a * s) == ["a**2", "a"]
    assert coeffs(s - a) == ["0", "1"]
    assert coeffs(a - s) == ["0", "-1"]
    assert coeffs(s / a) == ["1", "a**(-1)"]
    assert coeffs(s + t) == ["a", "3"]
    assert coeffs(s - t) == ["a", "-1"]
    assert coeffs(s + x) == ["a", "2"]
    assert coeffs(x + s) == ["a", "2"]
    assert coeffs(s - x) == ["a", "0"]
    assert coeffs(x - s) == ["-a", "0"]
    assert coeffs((s * x).expand()) == ["0", "a", "0", "0", "1"]
    assert coeffs(-s) == ["-1.0*a", "-1.0"]
    assert coeffs(+s) == ["a", "1"]
    assert coeffs(s ** 2) == coeffs(s * s)
    one = da.SDA(1) + da.svar(0)
    assert coeffs((1 / one).expand()) == coeffs((one ** -1).expand())
    assert coeffs((x / one)) == coeffs(da.promote(x) / one)
    with pytest.raises(ValueError):
        s / 0
    with pytest.raises(ValueError):
        s / 0.0
    with pytest.raises(ValueError):
        s / da.Expr(0)
    with pytest.raises(TypeError):
        s + "a"


def test_nda_expr_gives_sda(ab):
    a, _ = ab
    x = da.var(0)
    assert coeffs(a + x) == ["a", "1"]
    assert coeffs(x + a) == ["a", "1"]
    assert coeffs(x - a) == ["-a", "1"]
    assert coeffs(a - x) == ["a", "-1"]
    assert coeffs(x * a) == ["0", "a"]
    assert coeffs(a * x) == ["0", "a"]
    assert coeffs(x / a) == ["0", "a**(-1)"]
    assert isinstance(a / (1 + x), da.SDA)


@pytest.mark.parametrize("op", ["+=", "-=", "*=", "/="])
def test_inplace_operators(ab, op):
    a, _ = ab
    base = a + da.svar(0)
    for other in (da.SDA(2) + da.svar(1), 2.0, 2, a, 2 + da.var(1)):
        s = base.copy()
        before = id(s)
        want = {"+=": base + other, "-=": base - other,
                "*=": base * other, "/=": base / other}[op]
        exec(f"s {op} other")
        assert id(s) == before
        assert coeffs(s.expand()) == coeffs(want.expand())


def test_operator_env_checks(env):
    s = da.svar(0)
    da.clear()
    da.init(4, 3, 400)
    with pytest.raises(da.EnvError):
        s + 1
    with pytest.raises(da.EnvError):
        da.svar(0) + s


# T5.3

def test_functions(ab):
    a, _ = ab
    s = a + da.svar(0)
    for f in (da.sqrt, da.exp, da.log, da.sin, da.cos, da.tan, da.asin, da.acos, da.atan,
              da.sinh, da.cosh, da.tanh, da.erf):
        assert isinstance(f(s), da.SDA)
    assert coeffs(da.pow(s, 2)) == coeffs(s * s)
    assert isinstance(da.pow(s, 0.5), da.SDA)
    assert coeffs(da.exp(s))[:2] == ["exp(a)", "1.0*exp(a)"]
    for f in (da.asinh, da.acosh, da.atanh, da.abs):
        with pytest.raises(TypeError):
            f(s)
    x = da.exp(da.var(0))                    # the NDA overloads still come first
    assert isinstance(x, da.NDA)


def test_der_int_forms(ab):
    a, _ = ab
    s = a * da.svar(0) * da.svar(0)
    d = da.der(s, 0)
    assert coeffs(d) == ["0", "2.0*a"]
    out = da.SDA()
    da.da_der(s, 0, out)
    assert coeffs(out) == coeffs(d)
    i = da.int_(da.svar(1) * a, 1)
    assert i.element([0, 2, 0]) == a / 2
    da.da_int(s, 1, out)
    assert coeffs(out) == coeffs(da.da_int(s, 1))
    with pytest.raises(IndexError):
        da.da_der(s, 3)
    assert da.der is da.da_der and da.int_ is da.da_int
    assert da.substitute is da.da_substitute and da.compose is da.da_composition


def test_substitute_and_composition_lists(ab):
    a, _ = ab
    s = a * da.svar(0) + da.svar(1)
    out = da.SDA()
    da.substitute(s, 0, 2.0, out)
    assert coeffs(out) == coeffs(2.0 * a + da.svar(1))
    da.da_substitute(s, [0, 1], da.SDAList([da.SDA(1), da.SDA(a)]), out)
    assert coeffs(out) == ["2*a"]
    da.da_substitute(s, [0, 1], [da.SDA(1), da.SDA(a)], out)     # plain list
    assert coeffs(out) == ["2*a"]
    m = [s, s * s]
    pt = [da.SDA(1), da.SDA(2), da.SDA(0)]
    outs = da.SDAList([da.SDA(), da.SDA()])
    da.compose(m, pt, outs)
    assert [coeffs(v) for v in outs] == [["2 + a"], ["4 + 4*a + a**2"]]
    with pytest.raises(TypeError):
        da.compose(m, pt, [da.SDA(), da.SDA()])                  # output must be an SDAList
    with pytest.raises(ValueError):
        da.compose(m, pt[:2], outs)
    with pytest.raises(ValueError):
        ms = da.SDAList(m)
        da.compose(ms, pt, ms)                                   # output is also an input
    with pytest.raises(ValueError):
        da.da_substitute(s, [0, 0], [da.SDA(1), da.SDA(1)], out)
    vals = da.compose(m, [1.0, 2.0, 0.0])
    assert all(isinstance(e, da.Expr) for e in vals)
    assert [float(e.subs({a: 0.5})) for e in vals] == [2.5, 6.25]


# T5.4

def test_evaluate(ab):
    a, b = ab
    s = da.exp(a * da.svar(0) + 1)
    got = da.evaluate(s, {a: 0.3})
    ref = da.exp(0.3 * da.var(0) + 1)
    assert isinstance(got, da.NDA)
    assert got.length == ref.length
    for p, q in zip(got.coeffs(), ref.coeffs()):
        assert abs(p - q) <= 1e-13 * max(1.0, abs(q))
    x = 1.5 + 2.0 * da.var(0) + 0.5 * da.var(1) * da.var(2)
    back = da.evaluate(da.promote(x), {})
    assert (back.coeffs() == x.coeffs()).all()
    two = da.evaluate(a * da.svar(0) + b, [a, b], [2, 3.5])
    assert list(two.coeffs()) == [3.5, 2.0]
    with pytest.raises(ValueError):
        da.evaluate(s, {})                                      # a has no value
    with pytest.raises(ValueError):
        da.evaluate(s, [a, b], [1.0])


# T5.6

def test_simplify_expand_subs(ab):
    a, b = ab
    s = (a + 1) ** 2 * da.svar(0) + (da.sin(da.SDA(a)) ** 2).con + (da.cos(da.SDA(a)) ** 2).con
    e = s.expand()
    assert e is not s
    assert e.element(1) == a ** 2 + 2 * a + 1
    assert s.element(1) == (a + 1) ** 2                        # unchanged
    t = (a * 2 / 2 + b - b) * da.svar(1)
    assert t.simplify().element([0, 1, 0]) == a
    r = s.subs({a: 2})
    assert r.element(1) == 9
    assert s.subs({a: b}).element(1) == (b + 1) ** 2


# T5.5 (shared mode, see test_symengine_interop.py)

def test_symengine_objects(env):
    se = pytest.importorskip("symengine")
    assert da.symengine_interop_status()["mode"] == "shared"
    x, y = se.symbols("x y")
    s = da.SDA(x)
    assert da._rcp_address(s.con_symengine()) == da._rcp_address(x)
    assert isinstance(s.con_symengine(), se.Symbol)
    s.set_element([1, 0, 0], y)
    cs = s.coeffs(as_symengine=True)
    assert [da._rcp_address(c) for c in cs] == [da._rcp_address(x), da._rcp_address(y)]
    assert all(isinstance(c, da.Expr) for c in s.coeffs())
    s.con = y
    assert da._rcp_address(s.con_symengine()) == da._rcp_address(y)

    ex, ey = da.Expr.from_symengine(x), da.Expr.from_symengine(y)
    t = ex + da.svar(0)
    for got, want in [(t + y, t + ey), (y + t, ey + t), (t - y, t - ey), (y - t, ey - t),
                      (t * y, t * ey), (y * t, ey * t), (t / y, t / ey), (y / t, ey / t)]:
        assert isinstance(got, da.SDA)
        assert got.coeffs() == want.coeffs()
    u = t.copy()
    before = id(u)
    u += y
    u *= y
    u -= x
    u /= y
    assert id(u) == before
    assert u.coeffs() == ((((t + ey) * ey) - ex) / ey).coeffs()

    r = da.evaluate(x * da.svar(0) + y, {x: 2.0, y: 0.5})
    assert list(r.coeffs()) == [0.5, 2.0]
    r = da.evaluate(x * da.svar(0) + y, [x, ey], [2.0, 0.5])
    assert list(r.coeffs()) == [0.5, 2.0]
    assert (x * da.svar(0)).subs({x: y}).element(1) == ey
