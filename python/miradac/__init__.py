from collections.abc import Iterator
from contextlib import AbstractContextManager, contextmanager

from ._core import *
from ._core import __version__, HAS_SYMBOLIC


def _check_libsymengine() -> None:
    """The loaded libsymengine must be the pinned one we were built against (plan A.8)."""
    import hashlib
    import os

    loaded = _core._libsymengine_path()
    expected = _core._LIBSYMENGINE_EXPECTED
    if not expected:
        return  # a portable wheel: its bundled copy is the only one it can load
    if os.path.realpath(loaded) != expected:
        raise ImportError(f"miradac: loaded libsymengine {loaded}, built against {expected}")
    if _core._LIBSYMENGINE_SHA256:
        with open(loaded, "rb") as f:
            digest = hashlib.sha256(f.read()).hexdigest()
        if digest != _core._LIBSYMENGINE_SHA256:
            raise ImportError(
                f"miradac: {loaded} has sha256 {digest}, the pinned build {expected} "
                f"has {_core._LIBSYMENGINE_SHA256}")


if HAS_SYMBOLIC:
    from . import _core
    from ._core import _rcp_address, _rcp_use_count
    from ._sympy import from_sympy, to_sympy

    _check_libsymengine()
    Expr.to_sympy = to_sympy  # type: ignore[method-assign, assignment]
    Expr.from_sympy = staticmethod(from_sympy)  # type: ignore[method-assign]


class _Base:
    """``base[i]`` is ``var(i)`` in the current env."""

    def __getitem__(self, i: int) -> NDA:
        return var(i)

    def __len__(self) -> int:
        return nvars()


base = _Base()


@contextmanager
def _env_order(self: Env, n: int) -> Iterator[None]:
    """Temporarily truncate this env at order ``n``; nests correctly."""
    saved = self.current_order
    if not self.change_order(n):
        raise ValueError(f"order {n} exceeds the env's order {self.max_order}")
    try:
        yield
    finally:
        if saved == self.max_order:
            self.restore_order()
        else:
            self.change_order(saved)


Env.order = _env_order  # type: ignore[method-assign]


def order(n: int) -> AbstractContextManager[None]:
    """Temporarily truncate the current env at order ``n``; nests correctly."""
    return Env.current().order(n)
