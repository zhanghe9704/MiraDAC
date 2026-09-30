"""Multi-env API (plan T7.3, T7.4)."""

import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

import miradac as da

DATA = Path(__file__).resolve().parents[2] / "test"
symbolic = pytest.mark.skipif(not da.HAS_SYMBOLIC, reason="built without symbolic support")


@pytest.fixture
def envs():
    """A: the default env (order 2, 2 vars), current. B: a second env (order 4, 3 vars)."""
    da.init(2, 2, 200)
    b = da.Env(4, 3, 4000)
    yield da.Env.default(), b
    b.close()
    da.clear()


# --- T7.3 Env ----------------------------------------------------------------

def test_env_creation_and_properties(envs):
    a, b = envs
    assert da.Env.current() == a and da.Env.current() is a      # creating b kept a current
    assert (b.current_order, b.max_order, b.nvars, b.full_length) == (4, 4, 3, 35)
    assert (b.pool_size, b.count, b.remain, b.live_slots) == (4000, 0, 4000, 0)
    assert b.retired is False and a != b and len({a, b, da.Env.default()}) == 2
    assert (a == None, a != None, None in [a], a == 3) == (False, True, False, False)  # noqa: E711
    assert repr(b) == "Env(order=4, nvars=3, pool_size=4000)"
    x = da.var(0, env=b)
    assert b.count == 1 and x.env is b and da.var(0).env is a
    with pytest.raises(ValueError):
        da.Env(2, 2, 0)
    with pytest.raises(TypeError):
        da.Env()


def test_current_and_default_raise_without_env():
    da.clear()
    with pytest.raises(da.EnvError):
        da.Env.current()
    with pytest.raises(da.EnvError):
        da.Env.default()
    e = da.Env(2, 2, 10)
    with pytest.raises(da.EnvError):
        da.Env.current()                          # creating an env does not select it
    e.select()
    assert da.Env.current() is e
    e.close()
    assert da.Env.current() is e and e.retired    # the closed env stays selected
    with pytest.raises(da.EnvError):
        da.NDA()
    da.init(2, 2, 10)
    da.clear()


def test_with_nests_and_select(envs):
    a, b = envs
    with b:
        assert da.Env.current() is b
        with a:
            with b:
                assert da.var(0).env is b
            assert da.Env.current() is a
        with b:                                   # re-entrant
            assert da.Env.current() is b
        assert da.Env.current() is b
    assert da.Env.current() is a
    with pytest.raises(RuntimeError):
        with b:
            raise RuntimeError
    assert da.Env.current() is a
    b.select()
    assert da.NDA().env is b
    a.select()


def test_env_keyword_on_constructors(envs):
    a, b = envs
    made = [da.NDA(env=b), da.NDA(1.5, env=b), da.NDA.from_coeffs(np.ones(4), env=b),
            da.var(2, env=b), da.CNDA(env=b), da.CNDA(1 + 2j, env=b)]
    if da.HAS_SYMBOLIC:
        x, = da.symbols("x")
        made += [da.SDA(env=b), da.SDA(2.5, env=b), da.SDA(x, env=b), da.svar(2, env=b),
                 da.CSDA(env=b)]
    assert all(v.env is b for v in made)
    assert da.NDA(1.5).env is a
    assert da.exponents(b).shape == (35, 3) and da.exponents().shape == (6, 2)
    with pytest.raises(IndexError):
        da.var(2)                                 # a has 2 variables, b has 3
    with pytest.raises(TypeError):
        da.NDA(b)                                 # env is keyword-only


def test_close_with_live_vectors(envs):
    a, b = envs
    x = da.exp(da.var(0, env=b))
    c = da.CNDA(x, x)
    b.close()
    assert b.retired and repr(b) == "Env(retired)"
    b.close()                                     # idempotent
    for use in (lambda: x + 1.0, lambda: da.exp(x), lambda: x.con, lambda: c * c,
                lambda: da.NDA(env=b), lambda: b.count, lambda: b.select(),
                lambda: a.import_(x)):
        with pytest.raises(da.EnvError):
            use()
    assert x.env is b and da.Env.current() is a


def test_close_default_env_clears_it(envs):
    a, _ = envs
    x = da.var(0)
    a.close()
    assert a.retired
    with pytest.raises(da.EnvError):
        da.Env.default()
    with pytest.raises(da.EnvError):
        x * x
    da.init(2, 2, 10)                             # the fixture's clear() needs a default env


def test_default_env_without_user_vectors_is_retired_not_freed():
    da.init(2, 2, 10)
    a = da.Env.default()
    da.clear()
    assert a.retired
    da.init(3, 3, 20)
    assert da.Env.default() is not a and a.retired
    b = da.Env.default()
    b.close()
    assert b.retired


def test_failed_env_construction_keeps_current_env(envs):
    a, _ = envs
    with pytest.raises(MemoryError):
        da.Env(4, 3, 2**31)
    assert da.Env.current() == a


def test_process_exits_cleanly_after_close():
    code = textwrap.dedent("""
        import miradac as da
        da.init(3, 2, 100)
        e = da.Env(4, 3, 100)
        keep = [da.exp(da.var(0, env=e)), da.CNDA(1j, env=e)]
        if da.HAS_SYMBOLIC:
            keep.append(da.exp(da.svar(1, env=e)))
        e.close()
        try:
            keep[0] + 1.0
        except da.EnvError:
            print("EnvError")
        d = da.var(0)
        da.clear()
    """)
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    # stdout may start with the editable build's rebuild log
    assert r.returncode == 0 and r.stdout.endswith("EnvError\n"), r.stderr
    assert "Traceback" not in r.stderr and "terminate" not in r.stderr


def test_import_and_promote(envs):
    a, b = envs
    b2 = da.Env(4, 3, 100)
    try:
        x = 1.5 + da.var(0, env=b) * da.var(1, env=b) - 2.0 * da.var(2, env=b)
        y = b2.import_(x)
        assert y.env is b2 and b2.count == 1 and da.Env.current() is a
        z = b.import_(y)
        assert z.env is b and (z.coeffs() == x.coeffs()).all()
        c = da.CNDA(x, 2.0 * x)
        cc = b.import_(b2.import_(c))
        assert cc.env is b and (cc.real.coeffs() == x.coeffs()).all()
        assert (cc.imag.coeffs() == (2.0 * x).coeffs()).all()
        with pytest.raises(ValueError):
            a.import_(x)                          # different layout
        if da.HAS_SYMBOLIC:
            s = b2.promote(x)
            assert s.env is b2
            assert [str(e) for e in s.coeffs()] == [str(e) for e in da.promote(x).coeffs()]
            t = b.import_(s)
            assert t.env is b and [str(e) for e in t.coeffs()] == [str(e) for e in s.coeffs()]
            cs = b2.promote(c)
            assert cs.env is b2 and cs.real.env is b2
            cs2 = b.import_(cs)
            assert cs2.env is b
            assert [str(e) for e in cs2.imag.coeffs()] == [str(e) for e in cs.imag.coeffs()]
            with pytest.raises(ValueError):
                a.promote(x)
            with pytest.raises(ValueError):
                a.import_(s)
    finally:
        b2.close()


# --- T7.4 --------------------------------------------------------------------

def test_arithmetic_truncates_at_each_envs_order(envs):
    a, b = envs
    xa, xb = da.var(0), da.var(0, env=b)
    assert (xa * xa * xa).iszero()                # order 2
    assert (xb * xb * xb).element([3, 0, 0]) == 1.0
    assert (xb ** 5).iszero() and not (xb ** 4).iszero()
    assert da.exp(xa).n_element == 3 and da.exp(xb).n_element == 5


def test_mixing_envs_raises(envs):
    a, b = envs
    xa, xb = da.var(0), da.var(0, env=b)
    for op in (lambda: xa + xb, lambda: xb * xa, lambda: xb.__iadd__(xa),
               lambda: da.CNDA(xa, xb), lambda: da.CNDA(xb) * da.CNDA(xa),
               lambda: da.compose([xb, xb, xb], [xa, xa, xa], da.NDAList([da.NDA(env=b)] * 3))):
        with pytest.raises(da.EnvError):
            op()
    if da.HAS_SYMBOLIC:
        sa, sb = da.svar(0), da.svar(0, env=b)
        for op in (lambda: sa + sb, lambda: sb * xa, lambda: da.CSDA(sa, sb)):
            with pytest.raises(da.EnvError):
                op()


def test_env_order_affects_only_that_env(envs):
    a, b = envs
    with b.order(2):
        assert (b.current_order, a.current_order, da.current_order()) == (2, 2, 2)
        xb = da.var(0, env=b)
        assert (xb * xb * xb).iszero()
        with b.order(1):
            assert b.current_order == 1 and (xb * xb).iszero()
        assert b.current_order == 2
        with da.order(1):                         # the current env, a
            assert (a.current_order, b.current_order) == (1, 2)
        assert a.current_order == 2
    assert b.current_order == 4 and (xb * xb * xb).element([3, 0, 0]) == 1.0
    with pytest.raises(ValueError):
        with b.order(5):
            pass
    with pytest.raises(KeyError):
        with b.order(3):
            raise KeyError
    assert b.current_order == 4


def _snap(r, env):
    """Comparable form of a result; checks every DA result lives in env."""
    if isinstance(r, (list, tuple)) or type(r).__name__.endswith("List"):
        return [_snap(v, env) for v in r]
    if isinstance(r, (da.NDA, da.CNDA)) or (da.HAS_SYMBOLIC and isinstance(r, (da.SDA, da.CSDA))):
        assert r.env is env
        if isinstance(r, da.NDA):
            return list(r.coeffs())
        if isinstance(r, da.CNDA):
            return [list(r.real.coeffs()), list(r.imag.coeffs())]
        return str(r.real) + str(r.imag) if isinstance(r, da.CSDA) else [str(c) for c in r.coeffs()]
    if isinstance(r, np.ndarray):
        return r.tolist()
    return str(r) if da.HAS_SYMBOLIC and isinstance(r, da.Expr) else r


def _numeric_cases(b):
    x = 0.3 + 0.5 * da.var(0, env=b) + da.var(1, env=b) - 0.2 * da.var(2, env=b)
    y, z = da.sin(x), da.cos(x)
    xfile = 1.0 + da.var(0, env=b) + 2.0 * da.var(1, env=b) + 5.0 * da.var(2, env=b)
    v0, v1, v2 = (da.var(i, env=b) for i in range(3))
    m = [2.0 * v0 + 0.3 * v1 + 0.1 * v0 * v0, -0.4 * v0 + 1.5 * v1 + 0.2 * v1 * v1,
         v2 + 0.1 * v0 * v2]

    def out():
        return da.NDA(env=b)

    def outs(n=3):
        return da.NDAList([da.NDA(env=b) for _ in range(n)])

    def with_out(f, o):
        f(o)
        return o

    cases = {name: (lambda f=getattr(da, name): f(x))
             for name in ("sqrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan",
                          "sinh", "cosh", "tanh", "asinh", "atanh", "erf", "abs")}
    cases.update({
        "acosh": lambda: da.acosh(x + 1.5),
        "pow_int": lambda: da.pow(x, 3),
        "pow_float": lambda: da.pow(x, 0.3),
        "compare_da_with_file": lambda: da.compare_da_with_file(
            str(DATA / "exp_da.txt"), da.exp(xfile), 1e-14),
        "der": lambda: da.der(y, 1),
        "int_": lambda: da.int_(y, 1),
        "der_out": lambda: with_out(lambda o: da.da_der(y, 1, o), out()),
        "int_out": lambda: with_out(lambda o: da.da_int(y, 1, o), out()),
        "substitute_const": lambda: with_out(lambda o: da.da_substitute_const(y, 0, 0.5, o), out()),
        "substitute_da": lambda: with_out(lambda o: da.substitute(y, 0, z, o), out()),
        "substitute_float": lambda: with_out(lambda o: da.substitute(y, 0, 0.5, o), out()),
        "substitute_multi": lambda: with_out(lambda o: da.substitute(y, [0, 1], [z, x], o), out()),
        "substitute_lists": lambda: with_out(
            lambda o: da.substitute([x, y, z], [0, 1], [z, x], o), outs()),
        "compose_da": lambda: with_out(lambda o: da.compose([x, y, z], [z, x, y], o), outs()),
        "compose_float": lambda: da.compose([x, y, z], [0.1, -0.2, 0.3]),
        "compose_complex": lambda: da.compose([x, y, z], [0.1j, -0.2, 0.3]),
        "inv_map": lambda: with_out(lambda o: da.inv_map(m, 3, o), outs()),
        "devide_by_element": lambda: da.devide_by_element(y, z),
        "read_da_from_file": lambda: (lambda d: (da.read_da_from_file(
            str(DATA / "da_der.txt"), d), d))(out()),
        "compare_da_vectors": lambda: da.compare_da_vectors(y, y.copy()),
        "evaluate_map": lambda: da.evaluate_map(da.NDAList([x, y, z]), np.full((4, 3), 0.1)),
        "call": lambda: y([0.1, 0.2, 0.3]),
    })
    return cases


def _symbolic_cases(b):
    p, q = da.symbols("p q")
    s = p + 0.5 * da.svar(0, env=b) + q * da.svar(1, env=b) * da.svar(2, env=b)
    t = da.svar(1, env=b) * p

    def out():
        return da.SDA(env=b)

    def outs(n=2):
        return da.SDAList([da.SDA(env=b) for _ in range(n)])

    def with_out(f, o):
        f(o)
        return o

    cases = {"s_" + name: (lambda f=getattr(da, name): f(s))
             for name in ("sqrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan",
                          "sinh", "cosh", "tanh", "erf")}
    pt = lambda: [da.SDA(1, env=b), da.SDA(p, env=b), da.SDA(0, env=b)]  # noqa: E731
    cases.update({
        "s_pow_int": lambda: da.pow(s, 2),
        "s_pow_float": lambda: da.pow(s, 0.5),
        "s_der": lambda: da.der(s, 0),
        "s_int_": lambda: da.int_(s, 1),
        "s_der_out": lambda: with_out(lambda o: da.da_der(s, 1, o), out()),
        "s_int_out": lambda: with_out(lambda o: da.da_int(s, 2, o), out()),
        "s_substitute_const": lambda: with_out(lambda o: da.da_substitute_const(s, 0, 2.0, o), out()),
        "s_substitute_da": lambda: with_out(lambda o: da.substitute(s, 1, t, o), out()),
        "s_substitute_float": lambda: with_out(lambda o: da.substitute(s, 0, 2.0, o), out()),
        "s_substitute_multi": lambda: with_out(lambda o: da.substitute(s, [0, 1], [t, s], o), out()),
        "s_substitute_lists": lambda: with_out(
            lambda o: da.substitute([s, t], [0, 1], [t, s], o), outs()),
        "s_compose_da": lambda: with_out(lambda o: da.compose([s, t], pt(), o), outs()),
        "s_compose_float": lambda: da.compose([s, t], [1.0, 2.0, 0.5]),
        "s_evaluate": lambda: da.evaluate(da.exp(s), {p: 0.3, q: -0.2}),
    })
    return cases


def test_functions_on_other_env_vectors(envs):
    """Every T1.5, T2.2 and T5.3 function on env-B vectors, with A current, gives the
    same result as with B current."""
    a, b = envs
    cases = _numeric_cases(b)
    if da.HAS_SYMBOLIC:
        cases.update(_symbolic_cases(b))
    for name, f in cases.items():
        with b:
            ref = _snap(f(), b)
        with a:
            got = _snap(f(), b)
        assert got == ref, name
    assert cases["compare_da_with_file"]() is True


def _inv_map_coeffs():
    x, y = da.var(0), da.var(1)
    m = da.NDAList([2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y,
                    -0.4 * x + 1.5 * y + 0.2 * y * y])
    inv = da.NDAList([da.NDA(), da.NDA()])
    da.inv_map(m, 2, inv)
    return [list(v.coeffs()) for v in inv]


def test_inv_map_in_non_default_env():
    da.init(5, 2, 1000)
    ref = _inv_map_coeffs()
    e = da.Env(5, 2, 1000)
    with e:
        got = _inv_map_coeffs()
    da.clear()
    with e:                                       # also works with no default env
        again = _inv_map_coeffs()
    e.close()
    assert got == ref and again == ref


@symbolic
def test_erf_sda_in_non_default_env():
    # The last case of test/test_symbolic.cc.
    sx, sy = da.symbols("x y")

    def erf_coeffs():
        s = da.erf(sx + da.svar(0) * da.svar(0) + sy * da.svar(1))
        r = da.evaluate(s, [sx, sy], [0.3, -0.7])
        return list(r.coeffs())

    da.init(5, 2, 400, table=True)
    ref = erf_coeffs()
    e = da.Env(5, 2, 400, table=True)
    with e:
        got = erf_coeffs()
    e.close()
    da.clear()
    assert len(got) == len(ref)
    assert all(abs(g - r) < 1e-13 for g, r in zip(got, ref))
    assert ref[0] == pytest.approx(math.erf(0.3))
