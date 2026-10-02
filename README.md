# MiraDAC - Mixed-type Integrated RepresentAtion for Differential Algebra in C++

MiraDAC is a C++17 library for Truncated Power Series Algebra (TPSA) / Differential Algebra (DA).
It provides a single templated engine that supports:

- **NDA** (`da::DAVector<double>`) — numerical DA, zero dependencies, fast inner loops.
- **SDA** (`da::DAVector<SymEngine::Expression>`) — symbolic DA with SymEngine coefficients.
- **Interoperability** — mixed `NDA * SDA` expressions (result is SDA) and `da::evaluate(sda, values)` to convert SDA back to NDA.

Both flavors share one memory-pooled, monomial-indexed engine; the only difference is the coefficient type `T`.

## Features

- Templated DA engine monomorphized at compile time — no virtual dispatch, no runtime overhead.
- RAII `Pool<T>` with O(1) free-list; no per-vector heap allocation.
- Multiple independent `DAEnv`s (different order/nv) can coexist in one process.
- Environment guard (`DA_CHECK_ENV`, default on) catches cross-environment operations at runtime.
- `da_clear()` is safe to call while DA vectors are still in scope (see
  [Environment lifetime](#environment-lifetime)).
- Complex DA (`std::complex<NDA>`) and `cd_composition` ported from the numerical reference.
- Math functions: `sqrt, exp, log, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, pow, abs, erf`.
- Composition, substitution, derivative, integration, norm, zero-check.

## Build

### Numeric-only (no dependencies)

```bash
cmake -S . -B build -DWITH_SYMBOLIC=OFF
cmake --build build -j4
cd test && ../build/test/run_tests
```

### Full build with symbolic support

Symbolic support needs the pinned SymEngine build (see [SymEngine version](#symengine-version)).
Build it once, then point CMake at it:

```bash
eval "$(scripts/setup_symengine.sh --print-env)"   # sets SymEngine_DIR
cmake -S . -B build -DWITH_SYMBOLIC=ON -DSymEngine_DIR="$SymEngine_DIR"
cmake --build build -j4
cd test && ../build/test/run_tests
```

### Install

```bash
cmake --install build --prefix /usr/local
# or a custom prefix:
cmake --install build --prefix /opt/da
```

After install, downstream CMake projects can do:

```cmake
find_package(da REQUIRED)
target_link_libraries(myapp PRIVATE da::daStatic)
```

### Build options

| Option | Default | Description |
|--------|---------|-------------|
| `WITH_SYMBOLIC` | `ON` | Enable SDA / SymEngine support |
| `DA_CHECK_ENV` | `1` | Runtime cross-environment guard (set to `0` to disable) |
| `DA_BUILD_TESTS` | `ON` | Build the C++ test suite (forced `OFF` in a Python build) |
| `DA_BUILD_EXAMPLES` | `ON` | Build the examples (forced `OFF` in a Python build) |
| `DA_BUILD_BENCH` | `OFF` | Build `python/bench/bench_cpp`, the C++ baseline for the Python benchmark |
| `DA_IGNORE_SYMENGINE_PIN` | `OFF` | Accept any SymEngine and only warn when it is not the pinned build. Exists solely to benchmark against another SymEngine release (the symbolic benchmark's comparison with 0.15.0); never use it for a real build |

## SymEngine version

MiraDAC is pinned to one SymEngine build, recorded in `cmake/symengine_pin.txt`:

| Item | Value |
|---|---|
| SymEngine commit | `153b7e98f310bccaae586dab6b49284ccd5f4174` (v0.14.0 + 14 commits; reports version `0.14.0`) |
| symengine.py | `0.14.1` (the release that pins this commit) |
| Build options | `BUILD_SHARED_LIBS=ON`, `INTEGER_CLASS=gmp`, `WITH_SYMENGINE_THREAD_SAFE=OFF`, no FLINT/MPFR/MPC/LLVM, `Release` |
| Install prefix | `$HOME/.local/opt/symengine-0.14.0-153b7e98` (dedicated, so it never mixes with another SymEngine) |

**Why.** The Python package hands SymEngine objects to and from symengine.py by pointer (zero
copy). That is only safe when both load one and the same `libsymengine.so`. symengine.py only
supports the SymEngine commit it pins, so MiraDAC uses that commit too, and symengine.py is built
from source against it.

**Install.** Build the pinned SymEngine and build symengine.py against it into `.venv`
(needs `curl`, CMake, a C++ compiler and the GMP headers):

```bash
uv venv .venv --python 3.13                  # once, if .venv does not exist
scripts/setup_symengine.sh                   # [--prefix DIR] [--python PYTHON]
eval "$(scripts/setup_symengine.sh --print-env)"   # before every CMake configure
```

The script verifies both downloads by sha256, writes a stamp
`<prefix>/share/symengine/miradac-pin.txt` (commit, options, library sha256), and checks that the
installed `symengine_wrapper` links `libsymengine.so.0.14` from the prefix. A rerun skips whatever
is already correct. CMake refuses to configure MiraDAC against any other SymEngine (a different
version, or the v0.14.0 tag, which reports the same version but has no stamp), and the installed
`daConfig.cmake` requires `SymEngine 0.14.0 EXACT` from consumers.

**Do not `pip install symengine` from PyPI.** The PyPI wheel compiles its own private SymEngine
into `symengine_wrapper.so`; it can never share objects with MiraDAC, so interop falls back to
(slow) string conversion with a warning. If you install symengine.py yourself instead of using the
script, build it from source against the pin prefix (its build needs `cython` on `PATH`):

```bash
.venv/bin/python -m pip install cython setuptools
PATH=$PWD/.venv/bin:$PATH CMAKE_PREFIX_PATH=$HOME/.local/opt/symengine-0.14.0-153b7e98 \
    .venv/bin/python -m pip install --no-build-isolation --no-binary symengine symengine==0.14.1
```

**Checking interop.** `miradac.symengine_interop_status()` reports whether interop runs in
zero-copy or string mode and the result of each check: the symengine.py version, that its wrapper
links `libsymengine.so.0.14`, that the process maps exactly one `libsymengine` (the same file
MiraDAC uses), and a layout self-test. Set `MIRADAC_REQUIRE_SHARED_SYMENGINE=1` to turn a failed
check into an error instead of a warning.

**Benchmarking against another SymEngine.** `-DDA_IGNORE_SYMENGINE_PIN=ON` accepts any SymEngine
(for example 0.15.0) and turns the pin check into a warning. It exists solely so the symbolic
benchmark can compare releases; such a build must not be used with symengine.py.

**Upgrading the pin.**
1. Wait for a symengine.py release; take the SymEngine commit from its `symengine_version.txt`.
2. Update `cmake/symengine_pin.txt`: `SYMENGINE_COMMIT`, `SYMENGINE_VERSION`,
   `SYMENGINE_PY_VERSION`, and both sha256 values (`SYMENGINE_TARBALL_SHA256` of
   `https://github.com/symengine/symengine/archive/<commit>.tar.gz`, `SYMENGINE_PY_SDIST_SHA256`
   as listed on PyPI).
3. Run `scripts/setup_symengine.sh` on a fresh prefix, then the full C++ suite, the Python suite,
   and the symbolic benchmark against the old and the new pin.
4. Commit the pin change on its own.

## Quick start

### Numerical DA (NDA)

```cpp
#include "da/da.h"
using namespace da;

da_init(6, 6, 2000);   // order 6, 6 variables, pool of 2000 slots

NDA x = base[0];       // first basis variable
NDA y = base[1];

NDA f = sqrt(1.0 + x + 0.5*y);   // Taylor expansion of sqrt around origin
f.print();

std::vector<NDA> ivecs = {x, y};
std::vector<NDA> maps  = {x + 0.1*y, y - 0.05*x};
std::vector<NDA> ovecs(2);
da_composition(ivecs, maps, ovecs);   // compose f o maps

da_clear();
```

### Symbolic DA (SDA, requires WITH_SYMBOLIC=ON)

```cpp
#include "da/da.h"
using namespace da;
using SymEngine::Expression;
using SymEngine::symbol;

da_init(4, 2, 500, true);  // table=true needed for symbolic

Expression a = symbol("a");
Expression b = symbol("b");

SDA sx = base[0];     // symbolic first basis variable
SDA sy = base[1];

SDA g = a * sx + b * sy;
g.print();

// Differentiate
SDA dg_dx = da_der(g, 0);
dg_dx.print();

da_clear();
```

### Mixed NDA x SDA and evaluate

```cpp
#include "da/da.h"
using namespace da;
using SymEngine::Expression;
using SymEngine::symbol;

da_init(4, 2, 500, true);

Expression alpha = symbol("alpha");

NDA nx = base[0];         // NDA basis
SDA sg = alpha * nx;      // NDA * Expression -> SDA via promote()

// Evaluate at alpha=2.5: sg -> NDA
SymEngine::map_basic_basic vals = {{alpha.get_basic(),
                                    Expression(2.5).get_basic()}};
NDA result = evaluate(sg, vals);
result.print();

da_clear();
```

## Symbolic complex DA

`WITH_SYMBOLIC=ON` provides `std::complex<da::SDA>` (aliased as `da::CSDA`), representing a complex
DA map as a **real + imaginary pair of symbolic DA vectors**. This mirrors the numerical
`da::CNDA = std::complex<da::NDA>` but with `SymEngine::Expression` coefficients.

> **`CSDA` is the supported way to do complex arithmetic with symbolic DA.** Keep the imaginary
> unit in the `std::complex` wrapper, never inside a coefficient.
>
> SymEngine has its own imaginary unit (`SymEngine::I`), and putting it in an `SDA` coefficient
> will compile and produce plausible-looking results. It is **not supported**, for three reasons:
>
> - Real and imaginary parts can no longer be separated, which is the whole point of a complex
>   DA map (and the reason this design was chosen over SymEngine-native complex numbers).
> - `evaluate()` cannot convert such a vector back to numbers — `NDA` coefficients are `double`,
>   so SymEngine raises `Not Implemented`.
> - Mixing the two, e.g. an `I` buried in the real part of a `CSDA`, double-counts the imaginary
>   unit. This is untested and produces wrong answers.
>
> Pick one representation per computation, and let it be `CSDA`.

### Quick start

```cpp
#include "da/da.h"
#include <symengine/expression.h>
#include <symengine/symbol.h>
using namespace da;
using SE = SymEngine;

da_init(3, 2, 300);   // order 3, 2 variables

// Symbolic coefficients
SE::Expression a = SE::Expression(SE::symbol("a"));
SE::Expression b = SE::Expression(SE::symbol("b"));

// Build CSDA: z = (1.5 + a*x1) + i*(0.5 + b*x1)
// Use numeric constant parts so scalar DA functions can branch on them.
SDA re = SDA(SE::Expression(1.5)) + a * promote(da::base[0]);
SDA im = SDA(SE::Expression(0.5)) + b * promote(da::base[0]);
CSDA z(re, im);

// Complex DA arithmetic
CSDA z2  = z * z;
CSDA ez  = da::exp(z);     // exp, sqrt, log, asin, acos, atan, ...
CSDA sz  = da::sqrt(z);
SDA  mag = da::abs(z);     // sqrt(re^2 + im^2) as SDA

// Separate real and imaginary parts
const SDA& re_ez = get_real(ez);
const SDA& im_ez = get_imag(ez);

// Evaluate at a=0.3, b=0.2 -> CNDA
SE::map_basic_basic vals;
vals[SE::symbol("a")] = SE::real_double(0.3);
vals[SE::symbol("b")] = SE::real_double(0.2);
CNDA ez_numeric = evaluate(ez, vals);

// Compose complex maps
std::vector<CSDA> maps_in  = { z };
std::vector<CSDA> maps_sub = { CSDA(promote(da::base[0]), SDA(SE::Expression(0.0))) };
std::vector<CSDA> maps_out(1);
da::cd_composition(maps_in, maps_sub, maps_out);
```

### Supported operations

| Operation | Symbolic (`CSDA`) |
|---|---|
| `+`, `-`, `*`, `/` with `CSDA`, `double`, `complex<double>` | Yes |
| `exp`, `sqrt`, `log` | Yes |
| `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh` | Yes |
| `pow(z, int)`, `pow(z, double)` | Yes |
| `abs(CSDA)` → `SDA` | Yes (symbolic magnitude `sqrt(re²+im²)`) |
| `cd_composition` (all 3 overloads) | Yes |
| `promote(CNDA)` → `CSDA` | Yes |
| `evaluate(CSDA, map)` → `CNDA` | Yes |
| `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh` | **No** — see below |
| `SymEngine::I` inside a coefficient | **Not supported** — use `CSDA` |

The trigonometric and hyperbolic functions are not defined for complex DA in either flavor:
`CSDA` here has exactly the same function set as the numerical `CNDA`, which in turn matches
`ref/tpsa`. The scalar versions (`sin(SDA)`, `cos(SDA)`, …) do exist — only the
`std::complex<...>` overloads are absent. If you need them, build them from `exp`:
`sin(z) = (exp(i·z) - exp(-i·z)) / 2i`.

### Correctness model

All symbolic complex DA results are validated by **evaluation cross-check**: compute an
operation symbolically on `CSDA`, substitute concrete symbol values via `evaluate`, and verify
coefficient-by-coefficient agreement with the same operation run numerically on `CNDA` (within 1e-10).

### Note on constant parts

The scalar SDA functions compute the zeroth Taylor coefficient with SymEngine's own functions
(e.g. `exp(c0)`, `sqrt(c0)`, `sin(c0)`, `atan(c0)`, `erf(c0)`), so the constant part may be
**fully symbolic** — for example `exp(a + a*x1)` yields the seed `exp(a)`. No evaluation to a
double is required, for either `SDA` or `CSDA` transcendental functions.

The only constraints are the usual mathematical-domain conditions, checked symbolically:

- `sqrt` / `log` require a non-zero constant part. If it is (symbolically) zero, they emit a
  warning and return a NaN / -inf seed.
- `asin` / `acos` are defined via `atan` / `sqrt` and remain fully symbolic.

(Numeric-bound domain guards — e.g. the numerical `asin` requiring `|c0| <= 1` — apply to the
`NDA` path, not `SDA`.)

## Environment lifetime

A `DAVector` holds a raw pointer to the `DAEnv` it was created in, so the
environment has to stay readable for as long as the vector exists.

`da_clear()` handles this for you. If any DA vector still references the
default environment, the environment is **retired** rather than deleted: its
pools and monomial tables are freed immediately, and the (now empty) `DAEnv`
object is kept alive until the program exits, so the vectors' destructors
remain valid. The same applies when `da_init()` is called a second time.

```cpp
da::da_init(6, 6, 2000);
NDA x = base[0];
NDA y = sqrt(1.0 + x);

da::da_clear();     // fine: x and y are still in scope
                    // their destructors run safely afterwards
```

Computing with a vector after its environment has been cleared is still
invalid — the coefficient storage is gone. Only destruction is guaranteed.

For environments created with `da_make_env()`, use `da_destroy_env()` rather
than `delete`:

```cpp
da::DAEnv& env2 = da::da_make_env(4, 2, 500);
// ... work in env2 ...
da::da_destroy_env(env2);   // deletes it, or retires it if vectors remain
```

`delete &env2` is only correct when nothing references the environment any
more; `da_destroy_env()` checks for you (`env.live_slots()` reports the
count).

## Python

The package `miradac` binds NDA, CNDA, SDA, CSDA and the multi-environment API with
[nanobind](https://nanobind.readthedocs.io). Bound calls add about 50–130 ns to the C++ time, so
Python code driving DA arithmetic runs close to C++ speed (see `python/bench/REPORT.md`).

### Build

Everything Python lives in the virtual env `.venv` at the repository root.

```bash
uv venv .venv --python 3.13                         # once
scripts/setup_symengine.sh                          # pinned SymEngine + symengine.py (see above)
uv pip install --python .venv nanobind scikit-build-core pytest numpy sympy mypy
uv pip install --python .venv --reinstall pip       # a uv venv has no bin/pip
eval "$(scripts/setup_symengine.sh --print-env)"    # sets SymEngine_DIR
.venv/bin/pip install --no-build-isolation -Ceditable.rebuild=true -e .
.venv/bin/pytest python/tests -q
```

#### Numeric-only install

Without symbolic support (no SymEngine, no GMP; NDA, CNDA, `Env` and the map functions only,
`miradac.HAS_SYMBOLIC` is `False`):

```bash
pip install . -Ccmake.define.WITH_SYMBOLIC=OFF
```

No pinned SymEngine or `SymEngine_DIR` is needed, and the build ignores any SymEngine installed
on the machine.

#### Portable wheels

Linux x86_64 wheels (`manylinux_2_28`, CPython 3.10–3.14) need nothing installed on the target
machine: they carry the pinned SymEngine and GMP inside the wheel. Each version tag `v*` attaches
them to a GitHub Release (workflow `.github/workflows/wheels.yml`), and the workflow can also be
run by hand to get them as artifacts. Install one with `pip install miradac-<version>-<tag>.whl`.
To build them locally (needs Docker):

```bash
uvx cibuildwheel --platform linux --output-dir wheelhouse
```

A wheel's SymEngine is its own private copy, so it cannot share expressions with symengine.py
without copying: interop works in string mode (`miradac.symengine_interop_status()["mode"] ==
"string"`, with the reason given). For zero-copy interop, build from source after
`scripts/setup_symengine.sh`, as above.

#### Notes on the build

The editable install rebuilds the extension on `import miradac` whenever a source file changed
(set `SymEngine_DIR` in that shell too). `.venv/bin/pip wheel . --no-build-isolation -w dist`
builds a wheel; it links SymEngine from the pin prefix by RPATH, so it only runs on a machine with
that prefix.

### Example

```python
import numpy
import miradac as da

da.init(order=4, nvars=3, pool_size=1000)       # default env
x = 1.0 + da.var(0) + 2*da.var(1)               # NDA
y = da.exp(x); y += x*x
with da.order(2):                               # temporary truncation, nests correctly
    z = da.sin(y)
c = da.CNDA(x, y); w = da.exp(c)                # complex numeric
a, b = da.symbols("a b")                        # Expr
s = a*da.svar(0) + da.SDA(1.5)                  # SDA
v = da.evaluate(da.exp(s), {a: 0.3})            # -> NDA
m = da.NDAList([x, y, z]); pts = numpy.random.rand(10000, 3)
out = da.evaluate_map(m, pts)                   # (10000, 3) ndarray, one C++ loop

e = da.Env(order=10, nvars=2, pool_size=500)    # second env; does not change the current env
with e:
    q = da.var(0) * da.var(1)                   # lives in e
f = da.Env(order=4, nvars=3, pool_size=100)
x2 = f.import_(x)                               # copy into f (ValueError if layouts differ)
e.close()                                       # q now raises EnvError on use
```

`python/examples/` holds Python ports of the C++ examples. Every binding has a docstring
(`help(da.NDA)`), and the package ships type stubs (`_core.pyi`, `py.typed`) checked with
`mypy --strict`.

### Notes

- Every operation runs in the env of its first DA operand, whichever env is current. Vectors of
  different envs in one operation raise `miradac.EnvError`, as does a vector whose env was closed
  or cleared.
- Map-level functions (`compose`, `substitute`, `inv_map`, `cd_composition`, `evaluate_map`) take
  `NDAList`/`CNDAList`/`SDAList`/`CSDAList`, which hold the vectors in C++. A plain list works as
  an input but is copied element by element; output arguments must be one of these list types.
- CNDA⊕NDA, CSDA⊕SDA and CSDA⊕`Expr` (or `symengine.Basic`) work on both sides; CSDA⊕NDA and
  CSDA⊕CNDA do not (promote first: `z + da.promote(x)`).
- The GIL is never released: DA pools and SymEngine are not thread safe. Use processes for
  parallelism.
- Exceptions: `EnvError` (a `RuntimeError`), `ValueError` for invalid arguments and domain errors,
  `IndexError` for out-of-range indices, `RuntimeError` when a pool runs out.
- The stub `python/miradac/_core.pyi` is generated; after changing a binding, regenerate it with
  the command at the top of `python/stubgen_patterns.txt` (a test fails while it is stale).

## Julia

The package `MiraDAC.jl` (`julia/MiraDAC`) binds NDA, CNDA, SDA, CSDA and the multi-environment
API through a C API library, `libmiradac_c` (`capi/`), called with `ccall`. In-place operations
add a few tens of nanoseconds to the C++ time (see `julia/MiraDAC/bench/REPORT.md`). It needs
Julia 1.10 or later.

### Build

```bash
eval "$(scripts/setup_symengine.sh --print-env)"    # pinned SymEngine (see above)
cmake -S . -B build -G Ninja -DDA_BUILD_CAPI=ON -DSymEngine_DIR="$SymEngine_DIR"
cmake --build build                                  # -> build/capi/libmiradac_c.so
ctest --test-dir build -R capi                       # C API tests
julia julia/dev_setup.jl                             # points MiraDAC.jl at that library
julia --project=julia/MiraDAC -e 'using Pkg; Pkg.test()'
```

`julia/dev_setup.jl` stores the library path as the package's `libmiradac` preference
(`julia/MiraDAC/LocalPreferences.toml`, also in the `bench` and `docs` environments);
`MiraDAC.set_library!(path)` changes it later. A numeric-only library (`-DWITH_SYMBOLIC=OFF`,
no SymEngine) works too: `MiraDAC.HAS_SYMBOLIC` is then `false` and the symbolic tests are
skipped. The documentation (Documenter.jl) builds with
`julia --project=julia/MiraDAC/docs julia/MiraDAC/docs/make.jl` into `julia/MiraDAC/docs/build`.

### Example

```julia
using MiraDAC

init!(4, 3, 10_000)                           # default env: order 4, 3 variables, 10 000 slots
x = 1.0 + davar(1) + 2davar(2)                # NDA; variables are 1-based
y = exp(x); add!(y, y, x * x)                 # in place: y = y + x*x, no new slot
z = with_order(2) do                          # temporary truncation, nests correctly
    sin(y)
end
c = CNDA(x, y); w = exp(c)                    # complex numeric
a, b = dasymbols("a b")                       # SymExpr
s = a * sdavar(1) + SDA(1.5)                  # SDA
v = evaluate(exp(s), Dict(a => 0.3))          # -> NDA
m = NDAList([x, y, z]); pts = rand(3, 10_000) # points as columns
out = evaluate_map(m, pts)                    # 3×10000 Matrix, one C++ loop

e = DAEnv(10, 2, 500)                         # second env; does not change the current env
q = with_env(e) do
    davar(1) * davar(2)                       # lives in e
end
f = DAEnv(4, 3, 100)
x2 = import_vec(f, x)                         # copy into f (ArgumentError if layouts differ)
close(e)                                      # q now throws EnvError on use
```

`julia/MiraDAC/examples/` holds Julia ports of the C++ examples; every exported name has a
docstring (`?NDA`).

### Notes

- **Pool and GC.** A DA vector takes a slot in its env's fixed-size pool, and Julia frees it
  only when the garbage collector finalizes it. When a pool is full, an operation runs the GC
  and retries; `PoolExhaustedError` means more vectors are reachable than the pool holds. Size
  pools generously (each retry costs a young collection, ~100 µs) and use the in-place forms
  `add! sub! mul! div!` and `exp!(out, a)`, `sin!`, … in hot loops: they take no new slot.
- **Threads.** One env may be used by one Julia thread at a time (the C++ pools have no locks);
  different threads may use different envs. Symbolic objects are for one thread at a time.
- Every operation runs in the env of its first DA operand. Vectors of different envs in one
  operation, or of a closed or cleared env, throw `EnvError`.
- Indices are 1-based: `davar(1)` is the first variable. DA types are not `Number`s.
- With SymEngine.jl loaded, `SymExpr(::SymEngine.Basic)` and `SymEngine.Basic(::SymExpr)`
  convert through strings (the two use different SymEngine libraries).

## API summary

| Reference (numerical) | New (da::) |
|---|---|
| `DAVector` | `da::NDA` / `da::DAVector<double>` |
| `SymbDA::DAVector` | `da::SDA` / `da::DAVector<Expression>` |
| `da_init(o,nv,np)` | `da::da_init(o,nv,np)` |
| `da_clear()` | `da::da_clear()` |
| `da[i]` | `da::base[i]` |
| `ad_composition(...)` | `da::da_composition(ivecs, maps, ovecs)` |
| `eval(m)` / `eval_funs(...)` | `da::evaluate(sda, values)` |
| (n/a) | `da::promote(nda)` → SDA |
| (n/a) | `da::convert(vec, env)` / `env.import<T>(src_env, slot)` |

## Documentation

Generate API docs with Doxygen:

```bash
doxygen Doxyfile
# open doc/doxygen/html/index.html
```

## Benchmarks

See `examples/benchmark_composition.cc`. Build and run:

```bash
cmake -S . -B build -DWITH_SYMBOLIC=OFF
cmake --build build --target benchmark_composition -j4
./build/examples/benchmark_composition
```

The benchmark composes one DA vector of 6 bases with 6 DA vectors at orders 2, 4, 6 (and optionally 8, 10),
reporting time for both the MiraDAC engine and the reference `ref/tpsa` engine side by side.

## Acknowledgement

This work is supported by the U.S. Department of Energy, Office of Science, Office of Nuclear
Physics under contract DE-AC05-06OR23177 and contract No. 89243126CSC000213.

## License

MIT — see `LICENSE`.
