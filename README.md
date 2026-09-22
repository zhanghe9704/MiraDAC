# MiraDAC — Unified Numerical and Symbolic Differential Algebra in C++

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

SymEngine must be installed (or built from source). Point CMake at it:

```bash
cmake -S . -B build \
    -DWITH_SYMBOLIC=ON \
    -DSymEngine_DIR=/path/to/symengine/lib/cmake/symengine
cmake --build build -j4
cd test && LD_LIBRARY_PATH=/path/to/symengine/lib ../build/test/run_tests
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

## License

MIT — see `LICENSE`.
