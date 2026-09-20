# MiraDAC — Unified Numerical and Symbolic Differential Algebra

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
