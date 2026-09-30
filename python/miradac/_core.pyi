from collections.abc import Iterable, Iterator, Sequence
from typing import Annotated, Any, overload

import numpy
from numpy.typing import NDArray
import symengine


import contextlib
import typing

__version__: str

HAS_SYMBOLIC: bool = True

class EnvError(RuntimeError):
    """
    A DA vector was used after its env was closed or cleared, vectors of
    different envs met in one operation, or there is no current env.
    """

def init(order: int, nvars: int, pool_size: int, table: bool = False) -> None:
    """
    Create the default DA environment and make it current.

    DA vectors have nvars variables and are truncated at order; the pool
    holds pool_size vectors. table=True also builds the index/exponent
    lookup table. A previous default env is cleared first; its remaining
    vectors raise EnvError on use.
    """

def clear() -> None:
    """Clear the default env. Its remaining vectors raise EnvError on use."""

def count() -> int:
    """Number of vectors in use in the current env's pool."""

def remain() -> int:
    """Number of free vectors in the current env's pool."""

def poolsize() -> int:
    """Capacity of the current env's pool."""

def full_length() -> int:
    """Number of monomials up to the current order of the current env."""

def nvars() -> int:
    """Number of variables of the current env."""

def max_order() -> int:
    """Order the current env was created with."""

def current_order() -> int:
    """Order at which the current env truncates now (see change_order())."""

def get_eps() -> float:
    """Coefficients below eps in magnitude are dropped (shared by all envs)."""

def set_eps(x: float) -> None:
    """Set eps (see get_eps()); x must be positive."""

def change_order(n: int) -> bool:
    """
    Truncate the current env at order n <= max_order(); False if n is too
    large. Prefer the order(n) context manager.
    """

def restore_order() -> None:
    """Truncate the current env at its original order (max_order()) again."""

class Env:
    """
    A DA environment: variables, order and a pool of vectors.

    Env(order, nvars, pool_size, table=False) creates one without changing
    the current env. `with env:` makes it current for the block; the
    env= keyword of the constructors creates vectors in it. Every
    operation on a vector runs in the vector's own env.
    """

    def __init__(self, order: int, nvars: int, pool_size: int, table: bool = False) -> None:
        """Create an env; the current env does not change."""

    @staticmethod
    def current() -> Env:
        """The current env; EnvError if there is none."""

    @staticmethod
    def default() -> Env:
        """The env of the last init(); EnvError after clear()."""

    @property
    def current_order(self) -> int:
        """Order at which this env truncates now."""

    @property
    def max_order(self) -> int:
        """Order this env was created with."""

    @property
    def nvars(self) -> int:
        """Number of variables."""

    @property
    def full_length(self) -> int:
        """Number of monomials up to the current order."""

    @property
    def pool_size(self) -> int:
        """Capacity of the pool."""

    @property
    def count(self) -> int:
        """Number of vectors in use."""

    @property
    def remain(self) -> int:
        """Number of free vectors."""

    @property
    def live_slots(self) -> int:
        """Number of vectors in use, numeric and symbolic."""

    @property
    def retired(self) -> bool:
        """True once the env was closed or cleared."""

    def select(self) -> None:
        """Make this env current (until another is selected)."""

    def __enter__(self) -> Env:
        """Make this env current for the `with` block."""

    def __exit__(self, exc_type: object | None, exc: object | None, tb: object | None) -> None:
        """Restore the env that was current before the `with` block."""

    def close(self) -> None:
        """Release the pool; the env's vectors raise EnvError on use. Idempotent."""

    def change_order(self, n: int) -> bool:
        """
        Truncate this env at order n <= max_order; False if n is too large.
        Prefer the order(n) context manager.
        """

    def restore_order(self) -> None:
        """Truncate this env at its original order (max_order) again."""

    @overload
    def import_(self, v: NDA) -> NDA:
        """
        Copy v (NDA, CNDA, SDA or CSDA) into this env. ValueError if the
        two envs have different variables or orders.
        """

    @overload
    def import_(self, v: CNDA) -> CNDA: ...

    @overload
    def import_(self, v: SDA) -> SDA: ...

    @overload
    def import_(self, v: CSDA) -> CSDA: ...

    @overload
    def promote(self, v: NDA) -> SDA:
        """
        promote(v) (NDA to SDA, CNDA to CSDA) with the result in this env.
        ValueError if the two envs have different variables or orders.
        """

    @overload
    def promote(self, v: CNDA) -> CSDA: ...

    def __eq__(self, arg: object | None, /) -> bool:
        """True if value is the same env."""

    def __ne__(self, arg: object | None, /) -> bool:
        """True unless value is the same env."""

    def __hash__(self) -> int:
        """Hash of the env's identity."""

    def __repr__(self) -> str:
        """Env(order=..., nvars=..., pool_size=...) or Env(retired)."""

    def order(self, n: int) -> contextlib.AbstractContextManager[None]:
        """Temporarily truncate this env at order ``n``; nests correctly."""

class NDA:
    """
    A numeric DA vector: a truncated power series in nvars variables with
    float coefficients. It lives in the env it was created in.
    """

    @overload
    def __init__(self, *, env: Env | None = None) -> None:
        """NDA(): zero. NDA(x): the constant x. In the current env, or env."""

    @overload
    def __init__(self, x: float, *, env: Env | None = None) -> None: ...

    @staticmethod
    def from_coeffs(a: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], *, env: Env | None = None) -> NDA:
        """
        An NDA with the coefficients a (float64, at most full_length()
        values, in monomial order; see exponents()).
        """

    def copy(self) -> NDA:
        """A copy."""

    def __copy__(self) -> NDA:
        """A copy."""

    def __deepcopy__(self, memo: object) -> NDA:
        """A copy."""

    @property
    def con(self) -> float:
        """The constant part. Setting it resets the vector to that constant."""

    @con.setter
    def con(self, arg: float, /) -> None: ...

    @property
    def length(self) -> int:
        """Number of stored coefficients (up to the last non-zero one)."""

    @property
    def n_element(self) -> int:
        """Number of non-zero coefficients."""

    @property
    def nvars(self) -> int:
        """Number of variables of the vector's env."""

    @property
    def order(self) -> int:
        """Current truncation order of the vector's env."""

    def norm(self) -> float:
        """The largest coefficient magnitude."""

    def weighted_norm(self, w: float) -> float:
        """The largest magnitude of a coefficient times w**(its order)."""

    def iszero(self, eps: float | None = None) -> bool:
        """
        True if every coefficient is below eps in magnitude (default get_eps()).
        """

    def clean(self, eps: float | None = None) -> None:
        """Set coefficients below eps in magnitude (default get_eps()) to zero."""

    def reset(self) -> None:
        """Set every coefficient to zero."""

    @overload
    def element(self, i: int) -> float:
        """
        element(i): coefficient i in monomial order (see exponents()), 0 past
        length. element(exps): the coefficient of the monomial with the given
        exponent of each variable.
        """

    @overload
    def element(self, exps: Sequence[int]) -> float: ...

    def set_element(self, exps: Sequence[int], value: float) -> None:
        """Set the coefficient of the monomial with exponents exps."""

    def index_element(self, i: int) -> tuple[list[int], float]:
        """(exponents, coefficient) of coefficient i in monomial order."""

    def coeffs(self) -> NDArray[numpy.float64]:
        """A copy of the length stored coefficients, a float64 array."""

    def __call__(self, pt: Sequence[float]) -> float:
        """The value at the point pt (nvars coordinates)."""

    def __repr__(self) -> str:
        """NDA(order=..., nvars=..., nonzero=...)."""

    def __str__(self) -> str:
        """The table of terms that C++ prints."""

    @overload
    def __add__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __add__(self, arg: float, /) -> NDA: ...

    @overload
    def __add__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __add__(self, arg: Expr, /) -> SDA:
        """Return self + value."""

    @overload
    def __sub__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __sub__(self, arg: float, /) -> NDA: ...

    @overload
    def __sub__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __sub__(self, arg: Expr, /) -> SDA:
        """Return self - value."""

    @overload
    def __mul__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __mul__(self, arg: float, /) -> NDA: ...

    @overload
    def __mul__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __mul__(self, arg: Expr, /) -> SDA:
        """Return self * value."""

    @overload
    def __truediv__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __truediv__(self, arg: float, /) -> NDA: ...

    @overload
    def __truediv__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __truediv__(self, arg: Expr, /) -> SDA:
        """Return self / value."""

    @overload
    def __radd__(self, arg: float, /) -> NDA: ...

    @overload
    def __radd__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __radd__(self, arg: Expr, /) -> SDA:
        """Return value + self."""

    @overload
    def __rsub__(self, arg: float, /) -> NDA: ...

    @overload
    def __rsub__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rsub__(self, arg: Expr, /) -> SDA:
        """Return value - self."""

    @overload
    def __rmul__(self, arg: float, /) -> NDA: ...

    @overload
    def __rmul__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rmul__(self, arg: Expr, /) -> SDA:
        """Return value * self."""

    @overload
    def __rtruediv__(self, arg: float, /) -> NDA: ...

    @overload
    def __rtruediv__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rtruediv__(self, arg: Expr, /) -> SDA:
        """Return value / self."""

    @overload
    def __iadd__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __iadd__(self, arg: float, /) -> NDA:
        """self += value, in place."""

    @overload
    def __isub__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __isub__(self, arg: float, /) -> NDA:
        """self -= value, in place."""

    @overload
    def __imul__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __imul__(self, arg: float, /) -> NDA:
        """self *= value, in place."""

    @overload
    def __itruediv__(self, arg: NDA, /) -> NDA: ...

    @overload
    def __itruediv__(self, arg: float, /) -> NDA:
        """self /= value, in place."""

    def __neg__(self) -> NDA:
        """Return -self."""

    def __pos__(self) -> NDA:
        """Return a copy of self."""

    @overload
    def __pow__(self, arg: int, /) -> NDA: ...

    @overload
    def __pow__(self, arg: float, /) -> NDA:
        """Return self ** value."""

    @property
    def env(self) -> Env: ...

def var(i: int, *, env: Env | None = None) -> NDA:
    """
    The base vector of variable i (0 <= i < nvars), in the current env or
    env. base[i] is the same.
    """

@overload
def sqrt(arg: NDA, /) -> NDA:
    """Square root (NDA, CNDA, SDA or CSDA)."""

@overload
def sqrt(arg: CNDA, /) -> CNDA: ...

@overload
def sqrt(arg: SDA, /) -> SDA: ...

@overload
def sqrt(arg: CSDA, /) -> CSDA: ...

@overload
def exp(arg: NDA, /) -> NDA:
    """Exponential (NDA, CNDA, SDA or CSDA)."""

@overload
def exp(arg: CNDA, /) -> CNDA: ...

@overload
def exp(arg: SDA, /) -> SDA: ...

@overload
def exp(arg: CSDA, /) -> CSDA: ...

@overload
def log(arg: NDA, /) -> NDA:
    """Natural logarithm (NDA, CNDA, SDA or CSDA)."""

@overload
def log(arg: CNDA, /) -> CNDA: ...

@overload
def log(arg: SDA, /) -> SDA: ...

@overload
def log(arg: CSDA, /) -> CSDA: ...

@overload
def sin(arg: NDA, /) -> NDA:
    """Sine (NDA or SDA)."""

@overload
def sin(arg: SDA, /) -> SDA: ...

@overload
def cos(arg: NDA, /) -> NDA:
    """Cosine (NDA or SDA)."""

@overload
def cos(arg: SDA, /) -> SDA: ...

@overload
def tan(arg: NDA, /) -> NDA:
    """Tangent (NDA or SDA)."""

@overload
def tan(arg: SDA, /) -> SDA: ...

@overload
def asin(arg: NDA, /) -> NDA:
    """Arcsine (NDA, CNDA, SDA or CSDA)."""

@overload
def asin(arg: CNDA, /) -> CNDA: ...

@overload
def asin(arg: SDA, /) -> SDA: ...

@overload
def asin(arg: CSDA, /) -> CSDA: ...

@overload
def acos(arg: NDA, /) -> NDA:
    """Arccosine (NDA, CNDA, SDA or CSDA)."""

@overload
def acos(arg: CNDA, /) -> CNDA: ...

@overload
def acos(arg: SDA, /) -> SDA: ...

@overload
def acos(arg: CSDA, /) -> CSDA: ...

@overload
def atan(arg: NDA, /) -> NDA:
    """Arctangent (NDA, CNDA, SDA or CSDA)."""

@overload
def atan(arg: CNDA, /) -> CNDA: ...

@overload
def atan(arg: SDA, /) -> SDA: ...

@overload
def atan(arg: CSDA, /) -> CSDA: ...

@overload
def sinh(arg: NDA, /) -> NDA:
    """Hyperbolic sine (NDA or SDA)."""

@overload
def sinh(arg: SDA, /) -> SDA: ...

@overload
def cosh(arg: NDA, /) -> NDA:
    """Hyperbolic cosine (NDA or SDA)."""

@overload
def cosh(arg: SDA, /) -> SDA: ...

@overload
def tanh(arg: NDA, /) -> NDA:
    """Hyperbolic tangent (NDA or SDA)."""

@overload
def tanh(arg: SDA, /) -> SDA: ...

@overload
def asinh(arg: NDA, /) -> NDA:
    """Inverse hyperbolic sine (NDA, CNDA or CSDA)."""

@overload
def asinh(arg: CNDA, /) -> CNDA: ...

@overload
def asinh(arg: CSDA, /) -> CSDA: ...

@overload
def acosh(arg: NDA, /) -> NDA:
    """Inverse hyperbolic cosine (NDA, CNDA or CSDA)."""

@overload
def acosh(arg: CNDA, /) -> CNDA: ...

@overload
def acosh(arg: CSDA, /) -> CSDA: ...

@overload
def atanh(arg: NDA, /) -> NDA:
    """Inverse hyperbolic tangent (NDA, CNDA or CSDA)."""

@overload
def atanh(arg: CNDA, /) -> CNDA: ...

@overload
def atanh(arg: CSDA, /) -> CSDA: ...

@overload
def erf(arg: NDA, /) -> NDA:
    """Error function (NDA or SDA)."""

@overload
def erf(arg: SDA, /) -> SDA: ...

@overload
def abs(arg: NDA, /) -> float:
    """abs(NDA): the norm(), a float."""

@overload
def abs(arg: CNDA, /) -> float:
    """abs(CNDA): the larger of norm() of the real and the imaginary part."""

@overload
def abs(arg: CSDA, /) -> SDA:
    """abs(CSDA): the modulus sqrt(re**2 + im**2), as an SDA."""

@overload
def pow(arg0: NDA, arg1: int, /) -> NDA:
    """
    v raised to an int or float power (NDA, CNDA, SDA or CSDA); the same
    as v ** p.
    """

@overload
def pow(arg0: NDA, arg1: float, /) -> NDA: ...

@overload
def pow(arg0: CNDA, arg1: int, /) -> CNDA: ...

@overload
def pow(arg0: CNDA, arg1: float, /) -> CNDA: ...

@overload
def pow(arg0: SDA, arg1: int, /) -> SDA: ...

@overload
def pow(arg0: SDA, arg1: float, /) -> SDA: ...

@overload
def pow(arg0: CSDA, arg1: int, /) -> CSDA: ...

@overload
def pow(arg0: CSDA, arg1: float, /) -> CSDA: ...

def compare_da_with_file(filename: str, d: NDA, eps: float = 1e-15) -> bool:
    """compare_da_vectors of the NDA stored in filename and d."""

class NDAList:
    """
    A list of NDA held in C++, for the map-level functions (no per-element
    copies). A plain list is accepted as an input too, copied element by
    element; output arguments must be an NDAList.
    """

    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: NDAList, /) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[NDA], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[NDA]: ...

    @overload
    def __getitem__(self, arg: int, /) -> NDA: ...

    @overload
    def __getitem__(self, arg: slice, /) -> NDAList: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: NDA, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: NDA, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> NDA:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: NDAList, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: NDA, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: NDAList, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: NDA, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: NDA, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: NDA, /) -> None:
        """Remove first occurrence of ``arg``."""

@overload
def da_der(v: NDA, base_id: int) -> NDA:
    """
    Derivative of v with respect to variable base_id (alias der). With
    result given, it is written there and None is returned.
    """

@overload
def da_der(v: NDA, base_id: int, result: NDA) -> None: ...

@overload
def da_der(v: SDA, base_id: int) -> SDA: ...

@overload
def da_der(v: SDA, base_id: int, result: SDA) -> None: ...

@overload
def da_int(v: NDA, base_id: int) -> NDA:
    """
    Integral of v with respect to variable base_id (alias int_). With
    result given, it is written there and None is returned.
    """

@overload
def da_int(v: NDA, base_id: int, result: NDA) -> None: ...

@overload
def da_int(v: SDA, base_id: int) -> SDA: ...

@overload
def da_int(v: SDA, base_id: int, result: SDA) -> None: ...

@overload
def da_substitute_const(iv: NDA, base_id: int, x: float, ov: NDA) -> None:
    """ov = iv with the number x substituted for variable base_id."""

@overload
def da_substitute_const(iv: SDA, base_id: int, x: float, ov: SDA) -> None: ...

@overload
def da_substitute(iv: NDA, base_id: int, v: NDA, ov: NDA) -> None:
    """
    Substitute into iv and write the result to ov (alias substitute).

    base_id is one variable and v a DA vector or a number, or base_id is a
    list of variables and v a list of DA vectors, one per variable. iv and
    ov may also be lists (ov then an NDAList or SDAList of len(iv)).
    """

@overload
def da_substitute(iv: NDA, base_id: int, x: float, ov: NDA) -> None: ...

@overload
def da_substitute(iv: NDA, base_id: Sequence[int], v: NDAList, ov: NDA) -> None: ...

@overload
def da_substitute(ivecs: NDAList, base_id: Sequence[int], v: NDAList, ovecs: NDAList) -> None: ...

@overload
def da_substitute(iv: SDA, base_id: int, v: SDA, ov: SDA) -> None: ...

@overload
def da_substitute(iv: SDA, base_id: int, x: float, ov: SDA) -> None: ...

@overload
def da_substitute(iv: SDA, base_id: Sequence[int], v: SDAList, ov: SDA) -> None: ...

@overload
def da_substitute(ivecs: SDAList, base_id: Sequence[int], v: SDAList, ovecs: SDAList) -> None: ...

@overload
def da_composition(ivecs: NDAList, v: NDAList, ovecs: NDAList) -> None:
    """
    Compose a map with v (alias compose): ivecs[i](v[0], ..., v[nvars-1]).

    With DA vectors v, the results go to ovecs (a list of len(ivecs)) and
    None is returned; with numbers (float or complex) v, the list of
    values is returned. v has nvars() entries.
    """

@overload
def da_composition(ivecs: NDAList, v: Sequence[float]) -> list[float]: ...

@overload
def da_composition(ivecs: NDAList, v: Sequence[complex]) -> list[complex]: ...

@overload
def da_composition(ivecs: SDAList, v: SDAList, ovecs: SDAList) -> None: ...

@overload
def da_composition(ivecs: SDAList, v: Sequence[float]) -> list[Expr]: ...

def inv_map(ivecs: NDAList, dim: int, ovecs: NDAList) -> None:
    """
    Inverse of the map given by the first dim vectors of ivecs, written to
    the first dim vectors of ovecs (an NDAList). The map must have zero
    constant parts.
    """

def devide_by_element(t: NDA, b: NDA) -> NDA:
    """Coefficient-wise t / b (a coefficient of t stays where b's is zero)."""

def read_da_from_file(filename: str, d: NDA) -> bool:
    """
    Read an NDA written by str()/print into d; False if the file is
    incomplete.
    """

def compare_da_vectors(a: NDA, b: NDA, eps: float = 1e-15) -> bool:
    """
    True if a == b to eps: every coefficient of a - b is below eps times
    b's (below eps where b's is zero).
    """

der = da_der

int_ = da_int

substitute = da_substitute

compose = da_composition

def evaluate_map(map: NDAList, pts: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)]) -> NDArray[numpy.float64]:
    """
    Evaluate a map at many points in one C++ loop.

    pts is a float64 array of shape (N, nvars()); the result has shape
    (N, len(map)) and equals compose(map, pts[p]) for each point p.
    """

def exponents(env: Env | None = None) -> NDArray[numpy.int32]:
    """
    The exponents of every monomial, an int32 array of shape
    (full_length(), nvars()): row i belongs to coefficient i. For the
    current env, or env.
    """

class CNDA:
    """A complex numeric DA vector: a real and an imaginary NDA."""

    @overload
    def __init__(self, *, env: Env | None = None) -> None:
        """
        CNDA(): zero. CNDA(re, im=None): re + i*im from NDA parts of one env.
        CNDA(z): the constant z. Without parts, in the current env or env.
        """

    @overload
    def __init__(self, re: NDA, im: NDA | None = None) -> None: ...

    @overload
    def __init__(self, z: complex, *, env: Env | None = None) -> None: ...

    def copy(self) -> CNDA:
        """A copy."""

    def __copy__(self) -> CNDA:
        """A copy."""

    def __deepcopy__(self, memo: object) -> CNDA:
        """A copy."""

    @property
    def real(self) -> NDA:
        """The real part (a copy); assigning sets it."""

    @real.setter
    def real(self, arg: NDA, /) -> None: ...

    @property
    def imag(self) -> NDA:
        """The imaginary part (a copy); assigning sets it."""

    @imag.setter
    def imag(self, arg: NDA, /) -> None: ...

    def __repr__(self) -> str:
        """CNDA(order=..., nvars=..., nonzero=(re, im))."""

    def __str__(self) -> str:
        """The table of terms that C++ prints."""

    @overload
    def __add__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __add__(self, arg: float, /) -> CNDA: ...

    @overload
    def __add__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __add__(self, arg: NDA, /) -> CNDA:
        """Return self + value."""

    @overload
    def __sub__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __sub__(self, arg: float, /) -> CNDA: ...

    @overload
    def __sub__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __sub__(self, arg: NDA, /) -> CNDA:
        """Return self - value."""

    @overload
    def __mul__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __mul__(self, arg: float, /) -> CNDA: ...

    @overload
    def __mul__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __mul__(self, arg: NDA, /) -> CNDA:
        """Return self * value."""

    @overload
    def __truediv__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __truediv__(self, arg: float, /) -> CNDA: ...

    @overload
    def __truediv__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __truediv__(self, arg: NDA, /) -> CNDA:
        """Return self / value."""

    @overload
    def __radd__(self, arg: float, /) -> CNDA: ...

    @overload
    def __radd__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __radd__(self, arg: NDA, /) -> CNDA:
        """Return value + self."""

    @overload
    def __rsub__(self, arg: float, /) -> CNDA: ...

    @overload
    def __rsub__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rsub__(self, arg: NDA, /) -> CNDA:
        """Return value - self."""

    @overload
    def __rmul__(self, arg: float, /) -> CNDA: ...

    @overload
    def __rmul__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rmul__(self, arg: NDA, /) -> CNDA:
        """Return value * self."""

    @overload
    def __rtruediv__(self, arg: float, /) -> CNDA: ...

    @overload
    def __rtruediv__(self, arg: complex, /) -> CNDA: ...

    @overload
    def __rtruediv__(self, arg: NDA, /) -> CNDA:
        """Return value / self."""

    @overload
    def __iadd__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __iadd__(self, arg: float, /) -> CNDA: ...

    @overload
    def __iadd__(self, arg: complex, /) -> CNDA:
        """self += value, in place."""

    @overload
    def __isub__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __isub__(self, arg: float, /) -> CNDA: ...

    @overload
    def __isub__(self, arg: complex, /) -> CNDA:
        """self -= value, in place."""

    @overload
    def __imul__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __imul__(self, arg: float, /) -> CNDA: ...

    @overload
    def __imul__(self, arg: complex, /) -> CNDA:
        """self *= value, in place."""

    @overload
    def __itruediv__(self, arg: CNDA, /) -> CNDA: ...

    @overload
    def __itruediv__(self, arg: float, /) -> CNDA: ...

    @overload
    def __itruediv__(self, arg: complex, /) -> CNDA:
        """self /= value, in place."""

    def __neg__(self) -> CNDA:
        """Return -self."""

    def __pos__(self) -> CNDA:
        """Return a copy of self."""

    @overload
    def __pow__(self, arg: int, /) -> CNDA: ...

    @overload
    def __pow__(self, arg: float, /) -> CNDA:
        """Return self ** value."""

    @property
    def env(self) -> Env: ...

class CNDAList:
    """A list of CNDA held in C++; see NDAList."""

    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: CNDAList, /) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[CNDA], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[CNDA]: ...

    @overload
    def __getitem__(self, arg: int, /) -> CNDA: ...

    @overload
    def __getitem__(self, arg: slice, /) -> CNDAList: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: CNDA, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: CNDA, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> CNDA:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: CNDAList, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: CNDA, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: CNDAList, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: CNDA, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: CNDA, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: CNDA, /) -> None:
        """Remove first occurrence of ``arg``."""

@overload
def cd_composition(ivecs: NDAList, v: CNDAList, ovecs: CNDAList) -> None: ...

@overload
def cd_composition(ivecs: CNDAList, v: CNDAList, ovecs: CNDAList) -> None: ...

@overload
def cd_composition(ivecs: CNDAList, v: NDAList, ovecs: CNDAList) -> None: ...

@overload
def cd_composition(ivecs: SDAList, v: CSDAList, ovecs: CSDAList) -> None:
    """
    Compose a map with complex DA vectors: ovecs[i] = ivecs[i](v[0], ..., v[nvars-1]).

    ivecs is a real or complex map, v a real or complex list with nvars()
    entries; the results go to ovecs (a CNDAList or CSDAList of len(ivecs)).
    """

@overload
def cd_composition(ivecs: CSDAList, v: CSDAList, ovecs: CSDAList) -> None: ...

@overload
def cd_composition(ivecs: CSDAList, v: SDAList, ovecs: CSDAList) -> None: ...

def read_cd_from_file(filename: str, cd: CNDA) -> bool:
    """
    Read a CNDA written by str()/print into cd; False if the file is
    incomplete.
    """

def compare_cd_vectors(a: CNDA, b: CNDA, eps: float = 1e-15) -> bool:
    """compare_da_vectors of the real parts and of the imaginary parts."""

def compare_cd_with_file(filename: str, d: CNDA, eps: float = 1e-15) -> bool:
    """compare_cd_vectors of the CNDA stored in filename and d."""

class Expr:
    """
    A symbolic scalar (SymEngine::Expression), the coefficient type of SDA.
    """

    @overload
    def __init__(self, arg: int, /) -> None:
        """
        Expr(int): an exact integer. Expr(float). Expr(str): the parsed
        expression, e.g. Expr('a + 2*b').
        """

    @overload
    def __init__(self, arg: float, /) -> None: ...

    @overload
    def __init__(self, arg: str, /) -> None: ...

    @overload
    def __add__(self, arg: Expr, /) -> Expr: ...

    @overload
    def __add__(self, arg: int, /) -> Expr: ...

    @overload
    def __add__(self, arg: float, /) -> Expr:
        """Return self + value."""

    @overload
    def __radd__(self, arg: int, /) -> Expr: ...

    @overload
    def __radd__(self, arg: float, /) -> Expr:
        """Return value + self."""

    @overload
    def __sub__(self, arg: Expr, /) -> Expr: ...

    @overload
    def __sub__(self, arg: int, /) -> Expr: ...

    @overload
    def __sub__(self, arg: float, /) -> Expr:
        """Return self - value."""

    @overload
    def __rsub__(self, arg: int, /) -> Expr: ...

    @overload
    def __rsub__(self, arg: float, /) -> Expr:
        """Return value - self."""

    @overload
    def __mul__(self, arg: Expr, /) -> Expr: ...

    @overload
    def __mul__(self, arg: int, /) -> Expr: ...

    @overload
    def __mul__(self, arg: float, /) -> Expr:
        """Return self * value."""

    @overload
    def __rmul__(self, arg: int, /) -> Expr: ...

    @overload
    def __rmul__(self, arg: float, /) -> Expr:
        """Return value * self."""

    @overload
    def __truediv__(self, arg: Expr, /) -> Expr: ...

    @overload
    def __truediv__(self, arg: int, /) -> Expr: ...

    @overload
    def __truediv__(self, arg: float, /) -> Expr:
        """Return self / value."""

    @overload
    def __rtruediv__(self, arg: int, /) -> Expr: ...

    @overload
    def __rtruediv__(self, arg: float, /) -> Expr:
        """Return value / self."""

    @overload
    def __pow__(self, arg: Expr, /) -> Expr: ...

    @overload
    def __pow__(self, arg: int, /) -> Expr: ...

    @overload
    def __pow__(self, arg: float, /) -> Expr:
        """Return self ** value."""

    @overload
    def __rpow__(self, arg: int, /) -> Expr: ...

    @overload
    def __rpow__(self, arg: float, /) -> Expr:
        """Return value ** self."""

    def __neg__(self) -> Expr:
        """Return -self."""

    def __pos__(self) -> Expr:
        """Return a copy of self."""

    @overload
    def __eq__(self, arg: Expr, /) -> bool: ...

    @overload
    def __eq__(self, arg: int, /) -> bool: ...

    @overload
    def __eq__(self, arg: float, /) -> bool: ...

    @overload
    def __eq__(self, arg: object, /) -> bool:
        """Structural equality with an Expr, int or float."""

    @overload
    def __ne__(self, arg: Expr, /) -> bool: ...

    @overload
    def __ne__(self, arg: int, /) -> bool: ...

    @overload
    def __ne__(self, arg: float, /) -> bool: ...

    @overload
    def __ne__(self, arg: object, /) -> bool:
        """Negation of __eq__."""

    def __hash__(self) -> int:
        """SymEngine's hash, consistent with __eq__."""

    def __str__(self) -> str:
        """The expression as text, e.g. 'a + 2*b'."""

    def __repr__(self) -> str:
        """Expr('...')."""

    def __float__(self) -> float:
        """The numeric value; TypeError if free symbols remain."""

    def subs(self, mapping: dict[Any, Any]) -> Expr:
        """
        A new Expr with mapping ({symbol: value}, Expr, int or float)
        substituted.
        """

    def expand(self) -> Expr:
        """The expanded expression."""

    def diff(self, sym: Expr) -> Expr:
        """The derivative with respect to the symbol sym."""

    def free_symbols(self) -> set[Expr]:
        """The set of symbols in the expression."""

    def is_zero(self) -> bool:
        """True if the expression is zero (the test SDA uses)."""

    def simplify(self) -> Expr:
        """A simplified copy (the simplification SDA.simplify() uses)."""

    def to_symengine(self) -> object:
        """
        The expression as a symengine.py object: the same C++ object when
        symengine.py shares MiraDAC's libsymengine, else via a string (see
        symengine_interop_status()).
        """

    @staticmethod
    def from_symengine(obj: object) -> Expr:
        """The Expr of a symengine.Basic (see to_symengine())."""

    def to_sympy(self) -> typing.Any:
        """The sympy expression equal to this Expr."""

    @staticmethod
    def from_sympy(e: typing.Any) -> Expr:
        """The Expr equal to the sympy expression ``e``."""

def symbols(names: str) -> tuple[Expr, ...]:
    """
    Symbols named in a string separated by spaces or commas:
    a, b = symbols('a b').
    """

def symengine_interop_status() -> dict[str, Any]:
    """
    Whether symengine.py objects cross by pointer (mode 'shared') or by
    string (mode 'string'), with the result of each check and, in string
    mode, the first failed check (reason). The checks run once.
    """

def _libsymengine_path() -> str:
    """Realpath of the libsymengine this module is linked to."""

_LIBSYMENGINE_EXPECTED: str

_LIBSYMENGINE_SHA256: str

def _rcp_address(obj: object) -> int:
    """Address of the SymEngine object behind obj (for tests)."""

def _rcp_use_count(obj: object) -> int:
    """SymEngine reference count of the object behind obj (for tests)."""

class SDA:
    """
    A symbolic DA vector: a truncated power series whose coefficients are
    symbolic expressions (Expr). Mixed operations with NDA give an SDA.
    """

    @overload
    def __init__(self, *, env: Env | None = None) -> None:
        """
        SDA(): zero. SDA(x): the constant x (float, int as an exact integer,
        Expr or symengine.Basic). In the current env, or env.
        """

    @overload
    def __init__(self, x: float, *, env: Env | None = None) -> None: ...

    @overload
    def __init__(self, x: object, *, env: Env | None = None) -> None: ...

    def copy(self) -> SDA:
        """A copy."""

    def __copy__(self) -> SDA:
        """A copy."""

    def __deepcopy__(self, memo: object) -> SDA:
        """A copy."""

    @property
    def con(self) -> Expr:
        """
        The constant part, an Expr. Setting it (Expr, int, float or
        symengine.Basic) resets the vector to that constant.
        """

    @con.setter
    def con(self, arg: object, /) -> None: ...

    def con_symengine(self) -> object:
        """The constant part as a symengine.py object."""

    @property
    def length(self) -> int:
        """Number of stored coefficients (up to the last non-zero one)."""

    @property
    def n_element(self) -> int:
        """Number of non-zero coefficients."""

    @property
    def nvars(self) -> int:
        """Number of variables of the vector's env."""

    @property
    def order(self) -> int:
        """Current truncation order of the vector's env."""

    def iszero(self) -> bool:
        """True if every coefficient is zero."""

    def clean(self) -> None:
        """Set the coefficients for which Expr.is_zero() holds to exactly zero."""

    def reset(self) -> None:
        """Set every coefficient to zero."""

    @overload
    def element(self, i: int) -> Expr:
        """
        element(i): coefficient i in monomial order (see exponents()), 0 past
        length. element(exps): the coefficient of the monomial with the given
        exponent of each variable. An Expr.
        """

    @overload
    def element(self, exps: Sequence[int]) -> Expr: ...

    def set_element(self, exps: Sequence[int], value: object) -> None:
        """
        Set the coefficient of the monomial with exponents exps (Expr, int,
        float or symengine.Basic).
        """

    def index_element(self, i: int) -> tuple[list[int], Expr]:
        """(exponents, coefficient) of coefficient i in monomial order."""

    def coeffs(self, as_symengine: bool = False) -> list[Any]:
        """
        The length stored coefficients: a list of Expr, or of symengine.py
        objects with as_symengine=True.
        """

    def simplify(self) -> SDA:
        """A new SDA with every coefficient simplified."""

    def expand(self) -> SDA:
        """A new SDA with every coefficient expanded."""

    def subs(self, mapping: dict[Any, Any]) -> SDA:
        """
        A new SDA with mapping ({symbol: value}, Expr, int, float or
        symengine.Basic) substituted in every coefficient.
        """

    def __repr__(self) -> str:
        """SDA(order=..., nvars=..., nonzero=...)."""

    def __str__(self) -> str:
        """The table of terms that C++ prints."""

    @overload
    def __add__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __add__(self, arg: float, /) -> SDA: ...

    @overload
    def __add__(self, arg: int, /) -> SDA: ...

    @overload
    def __add__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __add__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __add__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __add__(self, arg: symengine.Basic, /) -> SDA:
        """Return self + value."""

    @overload
    def __sub__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __sub__(self, arg: float, /) -> SDA: ...

    @overload
    def __sub__(self, arg: int, /) -> SDA: ...

    @overload
    def __sub__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __sub__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __sub__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: symengine.Basic, /) -> SDA:
        """Return self - value."""

    @overload
    def __mul__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __mul__(self, arg: float, /) -> SDA: ...

    @overload
    def __mul__(self, arg: int, /) -> SDA: ...

    @overload
    def __mul__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __mul__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __mul__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: symengine.Basic, /) -> SDA:
        """Return self * value."""

    @overload
    def __truediv__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __truediv__(self, arg: float, /) -> SDA: ...

    @overload
    def __truediv__(self, arg: int, /) -> SDA: ...

    @overload
    def __truediv__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __truediv__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __truediv__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: symengine.Basic, /) -> SDA:
        """Return self / value."""

    @overload
    def __radd__(self, arg: float, /) -> SDA: ...

    @overload
    def __radd__(self, arg: int, /) -> SDA: ...

    @overload
    def __radd__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __radd__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __radd__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __radd__(self, arg: symengine.Basic, /) -> SDA:
        """Return value + self."""

    @overload
    def __rsub__(self, arg: float, /) -> SDA: ...

    @overload
    def __rsub__(self, arg: int, /) -> SDA: ...

    @overload
    def __rsub__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __rsub__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __rsub__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rsub__(self, arg: symengine.Basic, /) -> SDA:
        """Return value - self."""

    @overload
    def __rmul__(self, arg: float, /) -> SDA: ...

    @overload
    def __rmul__(self, arg: int, /) -> SDA: ...

    @overload
    def __rmul__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __rmul__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __rmul__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rmul__(self, arg: symengine.Basic, /) -> SDA:
        """Return value * self."""

    @overload
    def __rtruediv__(self, arg: float, /) -> SDA: ...

    @overload
    def __rtruediv__(self, arg: int, /) -> SDA: ...

    @overload
    def __rtruediv__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __rtruediv__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __rtruediv__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rtruediv__(self, arg: symengine.Basic, /) -> SDA:
        """Return value / self."""

    @overload
    def __iadd__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __iadd__(self, arg: float, /) -> SDA: ...

    @overload
    def __iadd__(self, arg: int, /) -> SDA: ...

    @overload
    def __iadd__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __iadd__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __iadd__(self, arg: symengine.Basic, /) -> SDA:
        """self += value, in place."""

    @overload
    def __isub__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __isub__(self, arg: float, /) -> SDA: ...

    @overload
    def __isub__(self, arg: int, /) -> SDA: ...

    @overload
    def __isub__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __isub__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __isub__(self, arg: symengine.Basic, /) -> SDA:
        """self -= value, in place."""

    @overload
    def __imul__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __imul__(self, arg: float, /) -> SDA: ...

    @overload
    def __imul__(self, arg: int, /) -> SDA: ...

    @overload
    def __imul__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __imul__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __imul__(self, arg: symengine.Basic, /) -> SDA:
        """self *= value, in place."""

    @overload
    def __itruediv__(self, arg: SDA, /) -> SDA: ...

    @overload
    def __itruediv__(self, arg: float, /) -> SDA: ...

    @overload
    def __itruediv__(self, arg: int, /) -> SDA: ...

    @overload
    def __itruediv__(self, arg: Expr, /) -> SDA: ...

    @overload
    def __itruediv__(self, arg: NDA, /) -> SDA: ...

    @overload
    def __itruediv__(self, arg: symengine.Basic, /) -> SDA:
        """self /= value, in place."""

    def __neg__(self) -> SDA:
        """Return -self."""

    def __pos__(self) -> SDA:
        """Return a copy of self."""

    @overload
    def __pow__(self, arg: int, /) -> SDA: ...

    @overload
    def __pow__(self, arg: float, /) -> SDA:
        """Return self ** value."""

    @property
    def env(self) -> Env: ...

@overload
def promote(nda: NDA) -> SDA:
    """
    The NDA as an SDA (a CNDA as a CSDA) in the same env; integer-valued
    coefficients become exact integers.
    """

@overload
def promote(cnda: CNDA) -> CSDA: ...

def svar(i: int, *, env: Env | None = None) -> SDA:
    """promote(var(i)): the symbolic base vector of variable i."""

class SDAList:
    """A list of SDA held in C++; see NDAList."""

    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: SDAList, /) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[SDA], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[SDA]: ...

    @overload
    def __getitem__(self, arg: int, /) -> SDA: ...

    @overload
    def __getitem__(self, arg: slice, /) -> SDAList: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: SDA, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: SDA, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> SDA:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: SDAList, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: SDA, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: SDAList, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: SDA, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: SDA, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: SDA, /) -> None:
        """Remove first occurrence of ``arg``."""

@overload
def evaluate(sda: SDA, values: dict[Any, float]) -> NDA: ...

@overload
def evaluate(sda: SDA, syms: Sequence[Any], vals: Sequence[float]) -> NDA: ...

@overload
def evaluate(sda: CSDA, values: dict[Any, float]) -> CNDA:
    """
    Substitute numbers for the symbols: an SDA gives an NDA, a CSDA a CNDA.

    values is a dict {symbol: float}, or give syms and vals as two
    sequences. Symbols are Expr or symengine.Basic; every symbol in a
    coefficient needs a value (ValueError otherwise).
    """

@overload
def evaluate(sda: CSDA, syms: Sequence[Any], vals: Sequence[float]) -> CNDA: ...

class CSDA:
    """A complex symbolic DA vector: a real and an imaginary SDA."""

    @overload
    def __init__(self, *, env: Env | None = None) -> None:
        """
        CSDA(): zero, in the current env or env. CSDA(re, im=None): re + i*im
        from SDA parts of one env.
        """

    @overload
    def __init__(self, re: SDA, im: SDA | None = None) -> None: ...

    def copy(self) -> CSDA:
        """A copy."""

    def __copy__(self) -> CSDA:
        """A copy."""

    def __deepcopy__(self, memo: object) -> CSDA:
        """A copy."""

    @property
    def real(self) -> SDA:
        """The real part (a copy); assigning sets it."""

    @real.setter
    def real(self, arg: SDA, /) -> None: ...

    @property
    def imag(self) -> SDA:
        """The imaginary part (a copy); assigning sets it."""

    @imag.setter
    def imag(self, arg: SDA, /) -> None: ...

    def __repr__(self) -> str:
        """CSDA(order=..., nvars=..., nonzero=(re, im))."""

    @overload
    def __add__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __add__(self, arg: float, /) -> CSDA: ...

    @overload
    def __add__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __add__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __add__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __add__(self, arg: symengine.Basic, /) -> CSDA:
        """Return self + value."""

    @overload
    def __sub__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: float, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __sub__(self, arg: symengine.Basic, /) -> CSDA:
        """Return self - value."""

    @overload
    def __mul__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: float, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __mul__(self, arg: symengine.Basic, /) -> CSDA:
        """Return self * value."""

    @overload
    def __truediv__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: float, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __truediv__(self, arg: symengine.Basic, /) -> CSDA:
        """Return self / value."""

    @overload
    def __radd__(self, arg: float, /) -> CSDA: ...

    @overload
    def __radd__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __radd__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __radd__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __radd__(self, arg: symengine.Basic, /) -> CSDA:
        """Return value + self."""

    @overload
    def __rsub__(self, arg: float, /) -> CSDA: ...

    @overload
    def __rsub__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rsub__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __rsub__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __rsub__(self, arg: symengine.Basic, /) -> CSDA:
        """Return value - self."""

    @overload
    def __rmul__(self, arg: float, /) -> CSDA: ...

    @overload
    def __rmul__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rmul__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __rmul__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __rmul__(self, arg: symengine.Basic, /) -> CSDA:
        """Return value * self."""

    @overload
    def __rtruediv__(self, arg: float, /) -> CSDA: ...

    @overload
    def __rtruediv__(self, arg: complex, /) -> CSDA: ...

    @overload
    def __rtruediv__(self, arg: SDA, /) -> CSDA: ...

    @overload
    def __rtruediv__(self, arg: Expr, /) -> CSDA: ...

    @overload
    def __rtruediv__(self, arg: symengine.Basic, /) -> CSDA:
        """Return value / self."""

    @overload
    def __iadd__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __iadd__(self, arg: float, /) -> CSDA: ...

    @overload
    def __iadd__(self, arg: complex, /) -> CSDA:
        """self += value, in place."""

    @overload
    def __isub__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __isub__(self, arg: float, /) -> CSDA: ...

    @overload
    def __isub__(self, arg: complex, /) -> CSDA:
        """self -= value, in place."""

    @overload
    def __imul__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __imul__(self, arg: float, /) -> CSDA: ...

    @overload
    def __imul__(self, arg: complex, /) -> CSDA:
        """self *= value, in place."""

    @overload
    def __itruediv__(self, arg: CSDA, /) -> CSDA: ...

    @overload
    def __itruediv__(self, arg: float, /) -> CSDA: ...

    @overload
    def __itruediv__(self, arg: complex, /) -> CSDA:
        """self /= value, in place."""

    def __neg__(self) -> CSDA:
        """Return -self."""

    def __pos__(self) -> CSDA:
        """Return a copy of self."""

    @overload
    def __pow__(self, arg: int, /) -> CSDA: ...

    @overload
    def __pow__(self, arg: float, /) -> CSDA:
        """Return self ** value."""

    @property
    def env(self) -> Env: ...

    def __str__(self) -> str:
        """The real part's table of terms, then the imaginary part's."""

class CSDAList:
    """A list of CSDA held in C++; see NDAList."""

    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: CSDAList, /) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[CSDA], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[CSDA]: ...

    @overload
    def __getitem__(self, arg: int, /) -> CSDA: ...

    @overload
    def __getitem__(self, arg: slice, /) -> CSDAList: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: CSDA, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: CSDA, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> CSDA:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: CSDAList, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: CSDA, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: CSDAList, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: CSDA, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: CSDA, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: CSDA, /) -> None:
        """Remove first occurrence of ``arg``."""
