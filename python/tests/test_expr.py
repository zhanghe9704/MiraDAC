import math
import random

import pytest

import miradac as da

pytestmark = pytest.mark.skipif(not da.HAS_SYMBOLIC, reason="built without symbolic support")


def test_ctors_and_repr():
    assert str(da.Expr(3)) == "3"
    assert float(da.Expr(0.25)) == 0.25
    assert repr(da.Expr("a + 2*b")) == "Expr('a + 2*b')"
    assert da.Expr(2**100) == da.Expr(str(2**100))
    with pytest.raises(ValueError):
        da.Expr("a +* b")


def test_symbols():
    a, b, c = da.symbols("a b c")
    assert (str(a), str(b), str(c)) == ("a", "b", "c")
    assert da.symbols("x, y") == da.symbols("x y")
    assert len(da.symbols("  z ")) == 1


def test_arithmetic():
    a, b = da.symbols("a b")
    assert a + b == da.Expr("a + b")
    assert 2 * a - b / 3 == da.Expr("2*a - b/3")
    assert 1 + a == a + 1
    assert 1.5 * a == da.Expr("1.5*a")
    assert 1 / a == a**-1
    assert a**b == da.Expr("a**b")
    assert 2**a == da.Expr("2**a")
    assert a**0.5 == da.Expr("a**0.5")
    assert -a == da.Expr("-a")
    assert (a - a) == 0


def test_float_eq_hash():
    a, b = da.symbols("a b")
    assert math.isclose(float(da.Expr("sin(1) + 2")), math.sin(1) + 2)
    with pytest.raises(TypeError):
        float(a + 1)
    assert a + b == b + a
    assert a != b
    assert hash(a + b) == hash(b + a)
    assert len({a, b, a + 0}) == 2
    assert (a == "a") is False


def test_methods():
    a, b = da.symbols("a b")
    e = (a + b) ** 2
    assert e.expand() == a**2 + 2 * a * b + b**2
    assert e.subs({a: 1.0, b: 2}) == da.Expr(9.0)
    assert e.subs({a: b}) == 4 * b**2
    assert e.diff(a) == 2 * (a + b)
    with pytest.raises(ValueError):
        e.diff(a + b)
    assert e.free_symbols() == {a, b}
    assert da.Expr(3).free_symbols() == set()
    assert (e - e.expand()).is_zero()
    assert not e.is_zero()
    s = da.Expr("sin(a)**2 + cos(a)**2")
    assert isinstance(s.simplify(), da.Expr)
    assert (a * 2 / 2).simplify() == a


def _random_expr(rng, syms, depth=3):
    if depth == 0 or rng.random() < 0.3:
        choice = rng.random()
        if choice < 0.6:
            return rng.choice(syms)
        if choice < 0.8:
            return da.Expr(rng.randint(-5, 5))
        return da.Expr(f"{rng.randint(1, 9)}/{rng.randint(2, 9)}")
    op = rng.choice(["add", "mul", "pow", "sin", "exp", "sub", "div"])
    x = _random_expr(rng, syms, depth - 1)
    if op == "sin":
        return da.Expr(f"sin({x})")
    if op == "exp":
        return da.Expr(f"exp({x})")
    if op == "pow":
        return x ** rng.randint(2, 3)
    y = _random_expr(rng, syms, depth - 1)
    if op == "div" and y.is_zero():
        return x
    return {"add": x + y, "mul": x * y, "sub": x - y, "div": x / y}[op]


def test_sympy_methods():
    sympy = pytest.importorskip("sympy")
    a, b = da.symbols("a b")
    s = (a + 2 * b).to_sympy()
    assert s == sympy.Symbol("a") + 2 * sympy.Symbol("b")
    assert da.Expr.from_sympy(sympy.Symbol("a") * 3) == 3 * a
    assert da.from_sympy(da.to_sympy(a)) == a


def test_sympy_round_trip():
    pytest.importorskip("sympy")
    rng = random.Random(1234)
    syms = list(da.symbols("a b c"))
    for _ in range(20):
        e = _random_expr(rng, syms)
        assert da.Expr.from_sympy(e.to_sympy()) == e, str(e)
