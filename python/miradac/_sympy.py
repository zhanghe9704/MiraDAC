"""sympy bridge for Expr, through strings (plan A.7, T4.2). sympy is imported lazily."""

from typing import Any

from ._core import Expr


def to_sympy(e: Expr) -> Any:
    """The sympy expression equal to the Expr ``e``."""
    import sympy

    return sympy.sympify(str(e))


def from_sympy(e: Any) -> Expr:
    """The Expr equal to the sympy expression ``e``."""
    return Expr(str(e))
