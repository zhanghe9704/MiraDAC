# MiraDAC — Python Binding Development Plan

**Goal:** A Python package `miradac` exposing the full MiraDAC library — NDA, CNDA, SDA, CSDA and
the multi-environment API — with call overhead low enough that Python code driving DA arithmetic
runs close to C++ speed.

This document has two parts:

1. **Part A — Design** (decisions, with the reason for each).
2. **Part B — Staged implementation plan** (ordered tasks, each with files, steps and an
   acceptance check, written so a small model can execute one task without making design
   decisions).

> **Ground rules for any implementer**
> - All design decisions are in Part A. Do not invent new ones. If a task cannot be done as
>   written (an API does not exist, a signature differs), stop and report; do not improvise.
> - Every stage must build and pass its tests before the next stage starts.
> - Do not modify C++ library code except in tasks explicitly marked **[C++]**. Every [C++] task
>   must keep the existing C++ test suite green (`ctest` in the C++ build dir).
> - `../tpsa/python-wrapper/src/tpsa_python.cc` (pybind11 binder of the ancestor library) is a
>   **reference for API coverage and naming only**. Do not copy its patterns blindly: it has no
>   env handling and converts Python lists to `std::vector` by copying every element.
> - Use the virtual env at `.venv` (repo root) for everything Python. Never `pip install` into
>   the system Python.

---

## Part A — Design

### A.1 Facts about the C++ library that drive the design

| Fact | Source | Consequence |
|---|---|---|
| `DAVector<T>` holds `DAEnv* env_` and `unsigned slot_` (public members). Data lives in a fixed pool, pre-allocated at `da_init`. | `davector.h:60-61`, `pool.h:93` | A Python object only needs to store the 16-byte `DAVector`. Slot pointers are stable, so memcpy in/out is safe. |
| Pool is fixed size; exhaustion throws `std::runtime_error("Pool::assign: Run out of vectors")`. | `pool.h:121-124` | Maps to a Python exception when a constructor runs out. CPython refcounting frees temporaries immediately, so pool pressure matches C++. Two C++ defects found in Stage 1 (the `noexcept` move ctor terminated the process when the pool was full; a pool that ran out could not be reused) were fixed on 2026-09-30, see T1.7. |
| Move ctor still calls `pool.assign()` for the moved-from object. | `davector.h:84-88` | Returning results by value costs one extra assign/free pair. Measure in Stage 2; optimize only if the gate fails (T2.6). |
| Binary DA×DA ops call `check_env`, which throws `std::logic_error` on mismatch. | `env.h:262-266` | Map to `miradac.EnvError`. |
| Constructors (`NDA()`, `NDA(1.0)`) take the **current thread-local env**. Ops like `a + 2.0` use `a`'s env. Library functions (`exp`, …) may create temporaries in the current env. | `davector.h:73,93`; `env.h:177-183` | Every bound call on a DA object selects that object's env for the duration of the call (A.4). This makes multi-env correct without auditing every C++ function. |
| `da_clear()` / `da_destroy_env()` **retire** an env that still has live vectors: memory freed, shell kept until exit. | `env.h:106-122, 200+` | `DAEnv*` pointers never dangle. A vector of a retired env must not be used: every bound call checks `env_->retired()`. |
| Global `da::base` belongs to the default env only. `inv_map` and `erf(SDA)` use it internally. | `base.h:44`, `src/base.cpp:191,203`, `src/functions.cpp:790` | Python gets `var(i)`, built per env. `inv_map` and `erf(SDA)` in a non-default env need a [C++] fix (T7.1). |
| `da_change_order(n)` / `da_restore_order()` act on the current env. Restore always returns to the **original** order (not a stack). | `da.h:293-305`, `src/layout.cpp:417-435` | The Python context manager saves the previous order and restores it, giving correct nesting. |
| SymEngine is built **without** `WITH_SYMENGINE_THREAD_SAFE` (non-atomic refcounts). The pinned build (A.8) keeps this setting. | `symengine_config.h` | Never release the GIL around symbolic code. |
| MiraDAC builds and passes its full C++ suite (805 assertions / 71 cases, plus the ASan run) against SymEngine commit `153b7e98` (0.14.0 + 14 commits), not only 0.15.0. | Verified 2026-09-29 in a scratch build; re-verified in Stage P with the suite at 844 assertions / 77 cases | The downgrade in A.8 is safe for the C++ library. |
| Pool and layout have no locks. | `pool.h` | Never release the GIL around any DA code (A.6). |
| Root `CMakeLists.txt` unconditionally adds `test/` and `examples/`. | `CMakeLists.txt:117,120` | Needs build options (T0.1). |

### A.2 Binding library: **nanobind** (≥ 2.0)

Chosen for performance:
- Lower per-call dispatch overhead than pybind11 (typically 2–3× less). This dominates cost at
  low DA orders, where one C++ op takes tens to hundreds of ns.
- The C++ object is stored inline in the Python instance (no separate heap allocation per object).
- First-class `nb::ndarray` for numpy transfer.
- Faster compile, smaller binary.

Build the module with `NOMINSIZE` (nanobind optimizes binding code for size by default; we want
speed) and `NB_STATIC`. LTO is tried in Stage 2 and kept only if the benchmark improves.

### A.3 Package layout

```
pyproject.toml                 # repo root: scikit-build-core backend, package "miradac"
CMakeLists.txt                 # root: adds python/ when SKBUILD is defined
python/
  CMakeLists.txt               # nanobind module target miradac._core
  src/
    module.cpp                 # NB_MODULE(_core, m): calls bind_* in order
    common.h                   # EnvError, EnvGuard, exception translators, helpers
    arith.h                    # template helpers binding operator sets
    bind_env.cpp               # env-level functions, Env class
    bind_nda.cpp               # NDA, NDAList, NDA functions and algorithms
    bind_cnda.cpp              # CNDA, CNDAList
    bind_expr.cpp              # Expr (SymEngine::Expression)       [symbolic only]
    se_bridge.h                # symengine.py zero-copy bridge + checks (A.8, A.9) [symbolic only]
    bind_sda.cpp               # SDA, SDAList, interop               [symbolic only]
    bind_csda.cpp              # CSDA, CSDAList                      [symbolic only]
  miradac/
    __init__.py                # re-exports _core; pure-Python sugar (order(), base proxy)
    _sympy.py                  # optional sympy bridge (string based)
    py.typed, _core.pyi        # stubs (Stage 8)
  stubgen_patterns.txt         # nanobind.stubgen patterns for _core.pyi (Stage 8)
  tests/                       # pytest
  bench/                       # benchmark scripts + C++ baseline
  examples/
```

Symbolic files compile only when `DA_WITH_SYMBOLIC` is defined; `miradac.HAS_SYMBOLIC` reports it.

### A.4 Env handling in every bound call (`EnvGuard`)

Every bound function or method that touches a DA object creates this guard on entry:

```cpp
struct EnvGuard {                       // python/src/common.h
    da::DAEnv* prev;
    explicit EnvGuard(da::DAEnv* e) {
        if (e->retired()) throw EnvError("DA environment has been cleared");
        prev = da::da_exchange_env(e);  // [C++] T0.2: set current, return previous (may be null)
    }
    ~EnvGuard() { da::da_exchange_env(prev); }
};
```

Cost: one thread-local read and write each way, about 1–2 ns. Benefits:
- Ops on vectors of env B work correctly while env A is current (temporaries created inside
  library functions land in the right env).
- Use after `clear()` raises `EnvError` instead of undefined behavior.

For functions with several DA arguments, guard on the first one; `check_env` inside C++ catches
mismatches in the rest. Constructors without a DA argument (`NDA(1.0)`) guard on the `env=`
keyword if given, else on the current env (raise `EnvError` if there is none).

### A.5 Performance rules for binding code

1. **Explicit overloads, no implicit conversions.** Never use `nb::implicitly_convertible`. For
   each operator, bind DA⊕DA first, then DA⊕float, then DA⊕complex/Expr.
2. **In-place operators** (`__iadd__` …) call the C++ compound operators directly: no allocation.
3. **Return by value** with the default policy (move). Never return `reference` to a temporary.
4. **Opaque list types** (`NDAList`, `CNDAList`, `SDAList`, `CSDAList`) via `nb::bind_vector`. All
   map-level algorithms (composition, substitute, inv_map) take these by reference, so there are
   no per-element copies. Plain Python lists are accepted too but copy each element; document this.
5. **Batch loops in C++.** Evaluating a map at N points is one call taking an `(N, nvars)` float64
   array and returning `(N, nmap)`; no per-point Python round trip.
6. **numpy transfer by memcpy** into/out of the pool slot. Return copies, not views: slot length
   changes under arithmetic and `change_order`, so a view would be unsafe.
7. **Base vectors are returned by value** (`var(i)` builds a new vector with `ad_var`). Cheap, and
   safe across `init()`/`clear()` cycles.

### A.6 GIL

Never released in v1. Pools have no locks, and SymEngine refcounts are not atomic. Releasing the
GIL would let another Python thread corrupt a pool. Threads therefore give no speed-up; use
processes. (Deferred: per-env locks plus a thread-safe SymEngine build.)

### A.7 Symbolic coefficients: own `Expr` class plus zero-copy symengine.py interop

`SDA` coefficients are `SymEngine::Expression`.
- **`miradac.Expr`** binds `SymEngine::Expression` directly. It is always available and needs
  nothing beyond MiraDAC's own SymEngine.
- **symengine.py interop** (`Expr.to_symengine()`, `Expr.from_symengine(obj)`, and SDA/CSDA
  coefficient access returning or accepting symengine.py objects). Zero-copy: the C++ pointer
  (`RCP<const Basic>`) is handed across. This is safe only when symengine.py and MiraDAC use **one
  and the same `libsymengine.so`**, which A.8 guarantees and verifies.
- If the verification fails (for example, the PyPI wheel is installed), interop falls back to
  string conversion and warns once. With `MIRADAC_REQUIRE_SHARED_SYMENGINE=1` it raises instead.
- sympy interop stays string based, in pure Python, off the hot path.

### A.8 SymEngine version pin

**Why a pin is needed.** The PyPI `symengine` wheel compiles its own SymEngine copy (0.14.1 wheel:
57 MB `symengine_wrapper.so`, no `libsymengine.so` dependency, bundled FLINT/MPFR/MPC/GMP). A
pointer from that copy used by MiraDAC's copy is undefined behavior. The fix is to build
symengine.py from source against the same shared SymEngine that MiraDAC uses. symengine.py only
supports the SymEngine commit it pins, so MiraDAC must use that commit too.

**The pin (single source of truth: `cmake/symengine_pin.txt`):**

| Item | Value |
|---|---|
| symengine.py version | `0.14.1` (latest on PyPI as of 2026-09-29) |
| SymEngine commit | `153b7e98f310bccaae586dab6b49284ccd5f4174` (from symengine.py's `symengine_version.txt`) |
| SymEngine reported version | `0.14.0` (the commit is v0.14.0 + 14 commits; its CMake still says 0.14.0) |
| SymEngine build options | `BUILD_SHARED_LIBS=ON`, `INTEGER_CLASS=gmp`, `WITH_SYMENGINE_THREAD_SAFE=OFF`, no FLINT/MPFR/MPC/LLVM, `CMAKE_BUILD_TYPE=Release` |
| Install prefix (default) | `$HOME/.local/opt/symengine-0.14.0-153b7e98` — a dedicated prefix, so its headers and CMake files never mix with another SymEngine version |

Because the commit reports the same version string as the v0.14.0 tag, a version check alone
cannot tell them apart. The build script therefore writes a **stamp file**,
`<prefix>/share/symengine/miradac-pin.txt`, one `KEY=value` per line (readable with CMake 3.16
`file(STRINGS)`, no JSON parser needed): `SYMENGINE_COMMIT`, `SYMENGINE_VERSION`, one line per build
option, and `LIB_SHA256` (sha256 of the installed `libsymengine.so.0.14.0`).

**Checks, in the order they run:**

| When | Check | On failure |
|---|---|---|
| MiraDAC CMake configure (C++ and Python builds) | `find_package(SymEngine 0.14.0 EXACT CONFIG)`; the stamp exists; stamp commit equals the pin; stamp options equal the pin | Configure error naming the found prefix and the fix command (only a warning with `DA_IGNORE_SYMENGINE_PIN=ON`, see below) |
| MiraDAC CMake configure | Embed the SymEngine prefix and `lib_sha256` into `_core` as constants | — |
| `import miradac` | `dladdr` on a SymEngine symbol gives the loaded `libsymengine` path; its realpath and sha256 equal the embedded ones (catches `LD_LIBRARY_PATH` swaps) | `ImportError` with both paths |
| First symengine.py interop call | 1. `symengine.__version__ == "0.14.1"`. 2. `symengine_wrapper` has `DT_NEEDED libsymengine.so.0.14` (the PyPI wheel has none). 3. `/proc/self/maps` holds exactly one `libsymengine` mapping, the same realpath as ours. 4. Layout self-test: create an object with symengine.py's exported `c2py` and read back the pointer (A.9). | String fallback plus a one-time warning (or `RuntimeError` if `MIRADAC_REQUIRE_SHARED_SYMENGINE=1`). `miradac.symengine_interop_status()` reports each check. |

**Measures that keep both sides on the same library:**
- `scripts/setup_symengine.sh` is the only supported way to build SymEngine and symengine.py. It
  builds the pinned commit into the pinned prefix, writes the stamp, then builds symengine.py
  0.14.1 from its sdist with `CMAKE_PREFIX_PATH=<prefix>` and verifies the result (check 2 and the
  wrapper's RUNPATH point into the prefix).
- symengine.py's own CMake asks for `SymEngine 0.14.0` **without** `EXACT`, so it would silently
  accept a newer SymEngine found on the default search path (for example 0.15.0 in `~/.local`).
  The script always sets `CMAKE_PREFIX_PATH` and verifies afterwards; never build symengine.py by hand.
- The PyPI wheel is **supported**, in string mode. It can never be zero-copy, whatever the version:
  it compiles SymEngine into `symengine_wrapper.so` with no exported symbols, so MiraDAC cannot link
  to that copy, and two copies in one process must not share objects. Its build options also differ
  (0.14.1 wheel: FLINT, LLVM, MPFR and MPC enabled; ours: none).
- `pyproject.toml`: `[tool.uv] no-binary-package = ["symengine"]`. This affects only this repo's
  own development environment, so developers always test the zero-copy path. It does not restrict
  users. The README tells users who want zero-copy to run `scripts/setup_symengine.sh` (or pip with
  `--no-binary symengine` plus `CMAKE_PREFIX_PATH` set to the pin prefix).
- `miradac` does not declare `symengine` as a hard dependency (interop is optional); it declares
  the extra `miradac[symengine]` with `symengine==0.14.1`.
- The MiraDAC C++ package config (`daConfig.cmake`) also requires `SymEngine 0.14.0 EXACT`, so C++
  consumers cannot mix versions either.
- CMake option `DA_IGNORE_SYMENGINE_PIN` (default `OFF`): when `ON`, any SymEngine version is
  accepted and a failed pin check only warns. It exists solely so the T0.5 symbolic benchmark can
  compare against SymEngine 0.15.0; such a build must never be used with symengine.py.
- Upgrading the pin is a documented procedure (T0.9): change `cmake/symengine_pin.txt` only after a
  new symengine.py release pins a new commit, then rerun the script and all tests.

### A.9 Crossing the symengine.py boundary without Cython

symengine.py ships `symengine_wrapper.pxd`:
- `cdef class Basic` has a single field `rcp_const_basic thisptr` and no `cdef` methods, so its
  instance layout is `PyObject_HEAD` followed by the RCP (one pointer).
- `c2py(rcp_const_basic)` is exported through `symengine_wrapper.__pyx_capi__["c2py"]` (a
  PyCapsule, verified present in 0.14.1).

MiraDAC's `_core` uses them from C++:
- **C++ → Python:** fetch the `c2py` capsule once, call it with the RCP.
- **Python → C++:** check `isinstance(obj, symengine.Basic)`, then read the RCP at offset
  `sizeof(PyObject)`, and copy it (which increments the refcount).
- The layout self-test (A.8, check 4) creates a symbol with `c2py` and confirms the pointer read
  back is the same. It runs once, so a layout change in a future symengine.py is caught before any
  real object crosses.

### A.10 Exceptions

| C++ | Python |
|---|---|
| `EnvError` (wrapper) and plain `std::logic_error` (from `check_env`) | `miradac.EnvError` (subclass of `RuntimeError`) |
| `std::invalid_argument`, `std::domain_error` | `ValueError` |
| `std::out_of_range` | `IndexError` |
| `std::runtime_error` (includes pool exhaustion) | `RuntimeError` |

Register one translator that rethrows and catches in this order: `invalid_argument`,
`domain_error`, `out_of_range`, `logic_error`; anything else falls through to nanobind defaults.

### A.11 Python API summary (target)

```python
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
x2 = e.import_(x)                               # ValueError if layouts differ
e.close()                                       # q now raises EnvError on use
```

---

## Part B — Staged implementation plan

Conventions for every task:
- **Build (editable):** `.venv/bin/pip install --no-build-isolation -Ceditable.rebuild=true -e .`
  (after the first build, the module rebuilds automatically on import).
- **Test:** `.venv/bin/pytest python/tests -q`.
- **SymEngine:** every build uses the pinned SymEngine from Stage P
  (`eval "$(scripts/setup_symengine.sh --print-env)"` before building).
- **C++ tests:** `cmake -S . -B build -G Ninja && cmake --build build && ctest --test-dir build`.
- "Acceptance" lists the exact checks. A task is done only when all of them pass.

### Stage P — SymEngine pin (runs before Stage 0)

Background and all values: A.8. **Do not modify `~/.local` without the user's explicit approval**
(only TP.4 touches it).

**TP.1 Pin file.**
- Files: `cmake/symengine_pin.txt`.
- Content (`KEY=value` lines): `SYMENGINE_COMMIT=153b7e98f310bccaae586dab6b49284ccd5f4174`,
  `SYMENGINE_VERSION=0.14.0`, `SYMENGINE_TARBALL_SHA256=<sha256 of
  https://github.com/symengine/symengine/archive/<commit>.tar.gz>`, `SYMENGINE_PY_VERSION=0.14.1`,
  `SYMENGINE_PY_SDIST_SHA256=<sha256 listed on PyPI for symengine-0.14.1.tar.gz>`, and one line per
  build option from the A.8 table. Compute both sha256 values by downloading the files once.
  *As built:* keys starting with `SYMENGINE_` are metadata; every other key is a SymEngine CMake
  option (`BUILD_SHARED_LIBS=ON`, `INTEGER_CLASS=gmp`, `WITH_SYMENGINE_THREAD_SAFE=OFF`,
  `WITH_FLINT=OFF`, `WITH_MPFR=OFF`, `WITH_MPC=OFF`, `WITH_LLVM=OFF`, `CMAKE_BUILD_TYPE=Release`),
  passed as `-DKEY=value` and copied into the stamp. Lines starting with `#` are comments.
- Acceptance: file present; both sha256 values match fresh downloads.

**TP.2 Setup script.**
- Files: `scripts/setup_symengine.sh`.
- Usage: `scripts/setup_symengine.sh [--prefix DIR] [--python PYTHON] [--print-env]`. Default prefix
  `$HOME/.local/opt/symengine-<version>-<first 8 chars of commit>`; default Python `.venv/bin/python`.
- Steps (fail fast with a clear message at each):
  1. Read every value from `cmake/symengine_pin.txt`.
  2. If `<prefix>/share/symengine/miradac-pin.txt` exists and matches the pin and `LIB_SHA256`
     matches the installed library, skip to step 5.
  3. Download the SymEngine tarball for the pinned commit; verify its sha256; configure with exactly
     the pinned options plus `BUILD_TESTS=OFF BUILD_BENCHMARKS=OFF`; build; install into the prefix.
  4. Write the stamp file (A.8).
  5. Download the symengine.py sdist; verify its sha256; verify its `symengine_version.txt` equals
     the pinned commit.
  6. Build and install it into the given Python with
     `CMAKE_PREFIX_PATH=<prefix> <python> -m pip install --no-build-isolation --no-binary symengine <sdist>`
     (install `cython` and `setuptools` into that Python first).
  7. Verify: `readelf -d` on `symengine_wrapper*.so` shows `NEEDED libsymengine.so.0.14` and a
     RUNPATH inside the prefix; `python -c "import symengine; assert symengine.__version__ == '0.14.1'"`.
  8. With `--print-env`, print `export SymEngine_DIR=<prefix>/lib/cmake/symengine`.
- Acceptance: a run on a clean prefix passes every step; a second run skips the SymEngine build;
  corrupting one byte of the stamp's `LIB_SHA256` makes step 2 rebuild.
- *As built:* steps 5–6 are also skipped when the installed symengine.py already passes step 7, so
  `--print-env` before every build is cheap. All progress goes to stderr; stdout carries only the
  `--print-env` line. SymEngine is installed with `CMAKE_INSTALL_LIBDIR=lib`. Step 6 bootstraps pip
  with `ensurepip` when the Python has none (a `uv venv` has no pip) and puts the Python's `bin/`
  first on `PATH`, because symengine.py's CMake looks for `cython` on `PATH`.

**TP.3 [C++] CMake pin check.**
- Files: `cmake/CheckSymEnginePin.cmake` (new), `CMakeLists.txt`, the package-config template in
  `cmake/` (the file that calls `find_dependency(SymEngine ...)`).
- Do: when `WITH_SYMBOLIC` is ON, `find_package(SymEngine <SYMENGINE_VERSION> EXACT REQUIRED CONFIG)`
  with the version read from `cmake/symengine_pin.txt`; then include `CheckSymEnginePin.cmake`, which
  locates `<SymEngine prefix>/share/symengine/miradac-pin.txt`, reads it with `file(STRINGS)`, and
  compares the commit and every option with the pin. On mismatch or a missing stamp:
  `message(FATAL_ERROR ...)` naming the found `SymEngine_DIR`, the expected commit, and
  `scripts/setup_symengine.sh --print-env` as the fix. Export the requirement to consumers:
  `find_dependency(SymEngine 0.14.0 EXACT)` in the package config.
- *As built:* `CheckSymEnginePin.cmake` reads the pin and does the `find_package` itself, without
  `REQUIRED`, so a version mismatch produces the message above (listing the rejected
  configs and their versions) instead of CMake's generic one. It also sets `DA_SYMENGINE_PREFIX`
  and `DA_SYMENGINE_LIB_SHA256` (from the stamp) for the constants embedded in `_core` (A.8), and
  `DA_SYMENGINE_FIND_VERSION`, which `daConfig.cmake.in` uses for `find_dependency`. Option
  `DA_IGNORE_SYMENGINE_PIN` (default `OFF`): when `ON`, `find_package` takes any version, the
  check only warns, and the package config drops `0.14.0 EXACT`. It exists solely for the T0.5
  benchmark comparison against SymEngine 0.15.0.
- Acceptance:
  1. Configure against the pinned prefix succeeds and `ctest` passes (805 assertions when written;
     844 assertions / 77 cases as of Stage P).
  2. Configure against `~/.local` (0.15.0) fails with the message.
  3. Configure against a SymEngine built from tag `v0.14.0` without a stamp fails with the message.
  4. The `test/consumer` project builds against the installed MiraDAC and the pinned SymEngine.

**TP.4 Migrate this machine (needs user approval before running).**
- Do: run `scripts/setup_symengine.sh`; rebuild MiraDAC against the pin; reinstall it into
  `~/.local` (same install command as before). Leave SymEngine 0.15.0 in `~/.local` untouched; the
  EXACT check keeps MiraDAC from using it.
- Acceptance: `ldd ~/.local/lib/libdaShared.so` shows `libsymengine.so.0.14` from the pin prefix;
  the out-of-tree `find_package(da)` check from the previous install still prints the same values.
- *Done 2026-09-30 (user approved).* Installed from commit `8f9c1cc`. The first install showed
  `libsymengine.so.0.14 => not found` for `libdaShared.so`: CMake strips the build RUNPATH on
  install, and the pin prefix is not on the loader path (0.15.0 had been in `~/.local/lib`, which
  consumers reached only through their own RUNPATH). `daShared` now sets
  `INSTALL_RPATH_USE_LINK_PATH`, so the installed library carries
  `RUNPATH=<pin prefix>/lib`. Both acceptance checks pass: `ldd` resolves
  `libsymengine.so.0.14` from the pin prefix, and the out-of-tree consumer, configured with no
  `SymEngine_DIR`, `CMAKE_PREFIX_PATH` or `LD_LIBRARY_PATH`, prints `1.000000 0.500000 -0.125000`,
  `a` and `(3.933052, 2.148636)`, as before.

**TP.5 Documentation.**
- Files: `README.md` (new section "SymEngine version"), `cmake/symengine_pin.txt` (header comment).
- Content: why the pin exists (A.8, first paragraph); the exact commit and version; how to install
  (`scripts/setup_symengine.sh`); the warning that `pip install symengine` from PyPI installs a
  wheel that cannot share objects, and the `--no-binary symengine` flag; how to read
  `miradac.symengine_interop_status()`; the **upgrade procedure**: wait for a symengine.py release;
  take the commit from its `symengine_version.txt`; update `cmake/symengine_pin.txt` (both sha256
  values too); rerun TP.2 on a fresh prefix, the full C++ suite, the Python suite, and the T0.5
  symbolic benchmark against old and new; commit the pin change on its own.
  *As built:* the section also documents `DA_IGNORE_SYMENGINE_PIN` as existing solely for the T0.5
  benchmark comparison, and the option is in the README's build-option table.
- Acceptance: a reader following only the README reaches a working build on a clean prefix.

### Stage 0 — Tooling and scaffolding

**T0.1 [C++] Build options.**
- Files: `CMakeLists.txt`.
- Do: add `option(DA_BUILD_TESTS "Build C++ tests" ON)` and
  `option(DA_BUILD_EXAMPLES "Build examples" ON)`; wrap `add_subdirectory(test)` and
  `add_subdirectory(examples)` in them. Add `if(SKBUILD) set(CMAKE_POSITION_INDEPENDENT_CODE ON) endif()`
  before the library targets, and at the end `if(SKBUILD) add_subdirectory(python) endif()`.
  When `SKBUILD` is set, force tests and examples OFF.
- Acceptance: default configure + `ctest` unchanged; configure with both options OFF builds only
  `daShared` and `daStatic`.
- *As built:* `DA_BUILD_BENCH` (T0.5) is declared next to the two options; all three are in the
  README's build-option table.

**T0.2 [C++] `da_exchange_env`.**
- Files: `include/da/env.h`, `src/env.cpp`, `test/test_env.cc`.
- Do: add `DAEnv* da_exchange_env(DAEnv* env) noexcept;`. It sets the thread-local current env to
  `env` (may be null) and returns the previous pointer (may be null). Doxygen comment in the style
  of the neighboring functions.
- Acceptance: new Catch2 test: exchange to env B returns A; exchange back returns B;
  after `exchange(nullptr)`, `da_current_env()` throws. Full `ctest` green.
  *As built:* suite at 851 assertions / 78 cases, ASan run green.
- (Done 2026-09-29.) Also fix the comment at `include/da/layout.h:129` to say that `restore_order()` returns to the
  original order given at construction (not to the value before the last `change_order()`), and
  the matching comment on `da_restore_order` in `da.h`.

**T0.3 Python environment.**
- Do: `uv venv .venv --python 3.13`;
  `uv pip install --python .venv nanobind scikit-build-core pytest numpy sympy`.
  Add `.venv/` and `__pycache__/` to `.gitignore`.
- Acceptance: `.venv/bin/python -c "import nanobind, numpy, sympy; print(nanobind.__version__)"`
  prints a version ≥ 2.0.
- *As built:* `.venv` is the one Stage P created (it already holds cython, setuptools and
  symengine.py 0.14.1 built from source); installed versions: nanobind 3.1.0,
  scikit-build-core 1.1.0, pytest 9.1.1, numpy 2.5.3, sympy 1.14.0. A `uv venv` has no
  `bin/pip` entry point (Stage P's `ensurepip` only created `pip3`), so
  `uv pip install --python .venv --reinstall pip` was run to make the Build command's
  `.venv/bin/pip` exist.

**T0.4 Packaging skeleton.**
- Files: `pyproject.toml` (root), `python/CMakeLists.txt`, `python/src/module.cpp`,
  `python/miradac/__init__.py`, `python/tests/test_smoke.py`.
- `pyproject.toml`: backend `scikit_build_core.build`;
  `requires = ["scikit-build-core>=0.10", "nanobind>=2.0"]`; `[tool.scikit-build]` with
  `cmake.build-type = "Release"`, `install.components = ["python"]`,
  `wheel.packages = ["python/miradac"]`, `cmake.define.WITH_SYMBOLIC = "ON"`.
- `python/CMakeLists.txt`: `find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)`;
  find nanobind via `execute_process(COMMAND ${Python_EXECUTABLE} -m nanobind --cmake_dir ...)`;
  `nanobind_add_module(_core NB_STATIC NOMINSIZE src/module.cpp)`;
  `target_link_libraries(_core PRIVATE daStatic)`; set the install RPATH to the SymEngine library
  directory; `install(TARGETS _core LIBRARY DESTINATION miradac COMPONENT python)`.
- `module.cpp`: `NB_MODULE(_core, m)` with `m.attr("__version__") = "0.1.0"` and
  `m.attr("HAS_SYMBOLIC")` set from `DA_WITH_SYMBOLIC`.
- `__init__.py`: `from ._core import *` and `from ._core import __version__, HAS_SYMBOLIC`.
- SymEngine is found through `SymEngine_DIR=<pin prefix>/lib/cmake/symengine` (set by the user or
  by `scripts/setup_symengine.sh --print-env`); the pin check from TP.3 runs automatically. Do not
  hard-code a home path in any committed file.
- Acceptance: the build command succeeds; `test_smoke.py` asserts
  `miradac.__version__ == "0.1.0"` and `miradac.HAS_SYMBOLIC is True`; the built wheel
  (`.venv/bin/pip wheel . --no-build-isolation -w /tmp/whl`) contains only the `miradac` package,
  with no C++ headers or libraries.
- *As built:* `pyproject.toml` also sets `build-dir = "build-py/{wheel_tag}"`, because
  scikit-build-core refuses `editable.rebuild` without a build dir (`build-py/` is covered by
  `.gitignore`'s `build-*/`). It also carries the A.8 items that live in this file: the extra
  `miradac[symengine]` (`symengine==0.14.1`) and `[tool.uv] no-binary-package = ["symengine"]`.
  The install RPATH is `${DA_SYMENGINE_PREFIX}/lib` (from `CheckSymEnginePin.cmake`). `_core` has
  no `NEEDED libsymengine` yet, since nothing in it calls SymEngine; that appears with Stage 4.
  An editable rebuild runs a full `cmake --install` into `build-py/.../install/`, but only the
  `miradac` package is mapped onto `sys.path`; the wheel holds just `miradac/__init__.py` and
  `miradac/_core*.so`.

**T0.5 C++ baseline benchmark.**
- Files: `python/bench/bench_cpp.cpp`, `python/bench/CMakeLists.txt`; root option
  `DA_BUILD_BENCH` (default OFF) adds it.
- Do: for (nvars, order) in {(3,4), (6,6), (6,10)}, time `a+b`, `a*b`, `a+=b`, `a*2.0`, `exp(a)`,
  and `da_composition` of an nvars-map with an nvars-map. Use dense random vectors, a warm-up,
  and the minimum over interleaved rounds (see *As built*). Print JSON `{case: ns_per_op}`.
- Also time symbolic cases at (3,3): `SDA*SDA` and `exp(SDA)` with one symbol per coefficient.
- Acceptance: writes `python/bench/baseline_cpp.json`; two consecutive runs agree within 5%. Run
  the symbolic cases once against the pinned SymEngine and once against 0.15.0 (`~/.local`), and
  record both. If the pin is more than 10% slower on any case, stop and report to the user.
- *As built:* case names are `n<nvars>o<order>/<op>` (`add`, `mul`, `iadd`, `mul_const`, `exp`,
  `composition`) and `sym_n3o3/<op>` (`mul`, `exp`). Every coefficient is set (via
  `da_element_orders`); a symbolic coefficient `i` is its own symbol (`a<i>`, `b<i>`). A run
  first doubles its repetition count until one sample takes ≥ 10 ms (the warm-up). Then 15
  rounds each run every case in turn, giving it 0.2 s of samples (at least one sample), and a
  case reports its **minimum** sample. This replaces "median of ≥ 7 runs of ≥ 0.2 s", which did
  not reproduce on the development machine (a loaded desktop, `powersave` governor, SMT siblings
  shared with other processes): two consecutive runs differed by more than 5% in 7 to 10 of 20
  cases, up to 108%. Load only ever adds time, and back-to-back runs of one case all fall into
  the same burst of load; interleaving spreads each case's samples over the whole run, short
  samples make clean ones likely, and the minimum keeps them. Times are **thread CPU time**
  (`CLOCK_THREAD_CPUTIME_ID`), not wall time. On Linux the process turns transparent huge pages
  off (`prctl(PR_SET_THP_DISABLE)`): with THP `always` and it on, one pair still differed by 6.3%
  on `n3o4/exp` (818 vs 770 ns); with it off, no case moved that much. The Python benchmark (T2.5) must use the same method and clock
  (`time.thread_time_ns`). `bench_cpp --symbolic` runs only the symbolic cases. A full run takes
  about 5 minutes (`n6o10/composition` alone is ~12 s per op).
  Build and run:
  `cmake -S . -B <dir> -G Ninja -DCMAKE_BUILD_TYPE=Release -DDA_BUILD_BENCH=ON -DDA_BUILD_TESTS=OFF -DDA_BUILD_EXAMPLES=OFF`,
  then `taskset -c 5 <dir>/python/bench/bench_cpp > python/bench/baseline_cpp.json`.
  **Reproducibility (2026-09-29):** the baseline was regenerated on an idle machine (load average
  1.4–1.5, no video call); it and the next run agree within 3.8% on every case (worst
  `n3o4/exp`, 742 vs 770 ns; 17 of 20 cases within 1%). The earlier baseline, taken with load
  average up to 17, was up to 14% slower (`n6o10/composition`). A pair taken during a video call
  plus OBS (every CPU 45–62% busy, zoom on cpu 5) differed in 9 of 20 cases by more than 5%, up
  to 29%: regenerate or compare only when the machine is idle.
  The minimum cannot remove load that lasts the whole run on the pinned core's SMT sibling
  (cpu 1 for cpu 5): with a video call there, `sym_n3o3/mul` rose from 35 400 to 57 000–63 000
  ns. Check the sibling is idle (`ps -eo pcpu,psr,comm --sort=-pcpu`) before regenerating the
  baseline or comparing against it.
  **Symbolic comparison** (4 interleaved `--symbolic` run pairs, same source, `taskset -c 5`,
  ns/op, ratio of medians over the pairs): `sym_n3o3/mul` pinned/0.15.0 = 1.00×,
  `sym_n3o3/exp` 0.95×. With the earlier median method: 1.00× and 1.05×. Rerun idle (2 pairs,
  load average ~1.1): `mul` 0.98× (34 128 / 34 390, 33 767 / 34 619 ns), `exp` 0.97× (70 281 /
  73 046, 71 034 / 72 308 ns). The pin is within 10% on both cases.

### Stage 1 — NDA core

**T1.1 Common infrastructure.**
- Files: `python/src/common.h`, `python/src/module.cpp`.
- Do: `struct EnvError : std::runtime_error`; `EnvGuard` exactly as in A.4; helper
  `da::DAEnv* current_env_or_throw()`; register `miradac.EnvError` via `nb::exception<EnvError>`
  and the translator from A.10.
- Acceptance: compiles; T1.7 tests cover each mapping.
- *As built:* `common.h` also declares `bind_env`/`bind_nda`. The translator is one
  `nb::register_exception_translator` whose payload is the `EnvError` type object.
  `current_env_or_throw()` reads the current env with a `da_exchange_env(nullptr)` round trip
  and also raises `EnvError` for a retired env.

**T1.2 Env-level functions (default env).**
- Files: `python/src/bind_env.cpp`.
- Bind: `init(order, nvars, pool_size, table=False)` → `da_init`; `clear()` → `da_clear`;
  `count()`, `remain()`, `poolsize()`, `full_length()`, `nvars()`, `max_order()` of the current env;
  `get_eps()` / `set_eps(x)`; `change_order(n) -> bool` (call `layout().change_order` directly, so
  nothing is printed to stdout); `restore_order()`; `current_order()`.
- Acceptance: tests for init/clear cycles (init, clear, init again with a different order).
- *As built:* `max_order()` is the order the env was built with and `current_order()` the
  present truncation (`layout().max_order()`, which `change_order` lowers). `Layout` keeps the
  original order private, so `max_order()` derives it from the pool's slot length, which is fixed
  at construction to C(nvars + order, order); no C++ change. `set_eps(x)` raises `ValueError`
  for `x <= 0` instead of calling `da_set_eps`, which prints a warning to stdout. All
  query functions raise `EnvError` when there is no current env. Tests: `python/tests/test_env.py`.

**T1.3 `NDA` class — construction and inspection.**
- Files: `python/src/bind_nda.cpp`.
- Bind (every method uses `EnvGuard`):
  - `NDA()`, `NDA(x: float, env=None)`,
    `NDA.from_coeffs(a: ndarray[float64, 1D], env=None)` (memcpy into the slot, set len);
  - `copy()`, `__copy__`, `__deepcopy__`;
  - `con` (property; setter via `reset_const`), `length`, `n_element`, `norm()`, `weighted_norm(w)`,
    `iszero(eps=None)`, `clean(eps=None)`, `reset()`;
  - `element(i: int)`, `element(exps: Sequence[int])`, `set_element(exps, value)`,
    `index_element(i) -> (list[int], float)`, `coeffs() -> ndarray[float64]` (copy of `length` values);
  - `nvars`, `order` (properties read from `env_->layout()`);
  - `__repr__` → `NDA(order=4, nvars=3, nonzero=12)`; `__str__` → the same text as `print()`
    (use `std::ostringstream` with the existing `operator<<` if there is one; otherwise format from
    `index_element`).
- Module level: `var(i, env=None)` returns a new NDA built with `detail::ad_var` for variable `i`;
  it raises `IndexError` if `i >= nvars`. In `__init__.py`, a `base` object whose `__getitem__`
  calls `var(i)` and whose `__len__` returns `nvars()`.
- Acceptance: tests compare `str(var(0))` with the C++ `base[0].print()` output;
  `NDA.from_coeffs(v.coeffs())` equals `v`.
- *As built:* the `env=` keyword on `NDA(...)`, `from_coeffs` and `var` is left to T7.3, which
  adds it together with the `Env` class it needs. `var(i)` calls the public `da::da_base(i)`
  (T7.1) under `EnvGuard`, not `detail::ad_var`. `length`, `n_element`, `nvars`, `order` and
  `con` are properties. `element(exps)` and `set_element` raise `ValueError` for a negative
  exponent (C++ would index out of bounds). `from_coeffs` raises `ValueError` for more than
  `full_length()` values. The printed header holds the slot number (`V [11]`), so the test
  compares with the C++ `base[0].print()` text (captured from a C++ run and kept in the test)
  after replacing the slot number. Tests: `python/tests/test_nda.py`.

**T1.4 Arithmetic.**
- Files: `python/src/arith.h`, `python/src/bind_nda.cpp`.
- Do: a template `bind_arith<DA>(nb::class_<DA>&)` that binds `__add__ __sub__ __mul__ __truediv__`
  for (DA, DA), then (DA, float); `__radd__ __rsub__ __rmul__ __rtruediv__` for float;
  `__iadd__ __isub__ __imul__ __itruediv__` for DA and float; `__neg__ __pos__`; `__pow__` for int,
  then float, calling `da::pow`. All with `nb::is_operator()`.
- Acceptance: tests for every operator against hand-computed coefficients; `x += y` keeps `id(x)`.
- *As built:* `binop(c, name, rname, op)` and `inplace(c, name, op)` in `arith.h` take a generic
  lambda, so each operator is one line. In-place operators return the same object through
  `nb::rv_policy::none`. C++ `operator/(DA, double)` used to call `std::exit(-1)` for a divisor
  below `DBL_MIN`, so `__truediv__(DA, float)` checked first. Fixed in C++ on 2026-09-30: it now
  throws `std::invalid_argument` like `/=`, which maps to `ValueError`, and the binding check was
  removed. `inv_map` of a singular map (a zero row, or a zero pivot, which the reference replaced
  with `DBL_MIN` and so returned huge or `inf` coefficients) now throws `std::invalid_argument`
  instead of exiting or returning garbage.

**T1.5 Math functions.**
- Bind module-level `sqrt exp log sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh erf abs pow`
  for NDA (`functions.h:29-47`), each guarded. Also bind `compare_da_with_file`.
- Acceptance: `exp`, `log`, `sqrt`, `pow(x,3)` and `pow(x,0.3)` match `test/exp_da.txt`,
  `log_da.txt`, `sqrt_da.txt`, `pow3_da.txt` and `pow0p3_da.txt` via `compare_da_with_file`.
- *As built:* `compare_da_with_file` cannot check `pow3_da.txt`: `read_da_from_file` reports
  success only when the file lists terms up to index `full_length() - 1`, and x³ ends at index 19
  of 35 (the C++ test ignores that return value and uses `compare_da_vectors`, a T2.2 binding).
  The test parses `pow3_da.txt` in Python and compares each of its 20 terms to 1e-14 relative;
  the other four files go through `compare_da_with_file` with eps 1e-14. A templated `unary<F>`
  wrapper binds the one-argument functions.

**T1.6 `order()` context manager.**
- Files: `python/miradac/__init__.py`.
- Do: `@contextmanager def order(n)`: save `current_order()`; call `change_order(n)` and raise
  `ValueError` if it returns False; yield; on exit, if the saved order equals `max_order()` call
  `restore_order()`, else call `change_order(saved)`.
- Acceptance: nested `with order(3): with order(2): ...` restores 3, then the original order; the
  order is also restored when the block raises.
- *As built:* as written; `max_order()` is the original order (T1.2 as-built note).

**T1.7 Error-path tests.**
- Acceptance tests: using an NDA after `clear()` raises `EnvError` and does not crash the process;
  exhausting a pool of 10 raises `RuntimeError`; `var(99)` raises `IndexError`.
- *As built:* `python/tests/test_errors.py`. The exhaustion test fills the pool with `NDA(1.0)`
  (constructed in place) and clears the env. Two C++ defects made other exhaustion paths
  unsafe. **Both fixed on 2026-09-30** (C++ tests `[numeric_move]`, `[pool]`, `[multienv][move]`;
  the former `xfail` test now passes, and the exhaustion test frees and reuses a slot):
  1. A result returned by value is moved into the Python object, and the `noexcept` move ctor
     (`davector.h`) calls `pool.assign()` for the source. With one free slot left, `var(0)` or
     `a + b` fills it and the move then throws inside `noexcept`: `std::terminate`. A strict
     `xfail` test (run in a subprocess) records this. T2.6 step 3 (construct in the instance)
     now avoids the move for operators and math functions, but `var()` and the algorithms
     still return by value. *Fix:* the move ctor hands over the slot and leaves the source
     owning none (`DAVector::no_slot`); the destructor skips it, and copy-assignment into such a
     vector allocates first. Move-assignment now swaps the env together with the slot (before,
     a move-assign between two envs mixed slot indices of different pools). `clear()` also
     leaves the vector with no slot instead of taking a dummy one.
  2. Once `Pool::assign` has thrown, `head_` sits on the sentinel and `Pool::free` never moves
     it, so the freed slots are lost and `Pool::remain()` reads past `free_` (ASan
     heap-buffer-overflow, reproduced in plain C++). *Fix:* `Pool::free` restarts the free list
     when it is empty.
  The `logic_error` mapping is covered by a live vector meeting one of a cleared env (`y + x`
  after `clear(); init(...)`), `domain_error` by `asin(NDA(2.0))`, `invalid_argument` by
  `var(0) / 0.0`.

### Stage 2 — NDA algorithms, bulk API, performance gate

**T2.1 Opaque lists.**
- Do: `nb::bind_vector<std::vector<NDA>>(m, "NDAList")`. Check whether nanobind's `bind_vector`
  accepts a Python `list` where `std::vector<NDA>&` is expected. If it does not, add explicit
  overloads that build a temporary `std::vector<NDA>` from a list. Record the result in a comment.
- Acceptance: tests pass both an `NDAList` and a plain list to one algorithm from T2.2.
- *As built:* `bind_vector` registers `implicitly_convertible<iterable, NDAList>` itself, so a
  plain list is accepted wherever an `NDAList&` is expected (copied element by element into a
  temporary); no extra overloads. A temporary cannot carry results back, so every list **output**
  argument is `.noconvert()` and a plain list there raises `TypeError`. `__getitem__` returns a
  copy (nanobind's default policy); `lst[i] += x` still works through `__setitem__`. Known gap:
  `NDAList` copies elements through the C++ copy ctor outside any `EnvGuard`, so copying a vector
  of a cleared env (`lst[0]`, `append`, list conversion) reads a released pool's free list
  (undefined behavior). Fixing it needs a **[C++]** guard in `Pool::assign` for a released pool
  (not scheduled).

**T2.2 Algorithms.**
- Bind both the output-argument and the returning forms where C++ has both: `da_der`, `da_int`,
  `da_substitute_const`, `da_substitute` (all four NDA overloads, `da.h:343-390`), `da_composition`
  (with an NDAList, a list of float or a list of complex; `da.h:394-430`), `inv_map`,
  `devide_by_element`, `read_da_from_file`, `compare_da_vectors`, `compare_da_with_file`.
- Python names drop the `da_` prefix (`der`, `int_`, `substitute`, `compose`, `inv_map`); keep the
  C++ names as aliases.
- Acceptance: a test for each against the `test/*.txt` files (`substitute_*.txt`, `da_der.txt`,
  `da_int.txt`, `da_composition_*.txt`, `bunch_substitution_*.txt`), mirroring
  `../tpsa/python-wrapper/tests/tests.py`.
- *As built:* the functions are bound under the C++ names with overloads; `der`, `int_`,
  `substitute`, `compose` are the same function objects (`da.der is da.da_der`). Only
  `da_der`/`da_int` have both forms in C++. The output-argument forms are kept for NDA outputs;
  `da_composition` with a list of float or complex **returns** the list of values, since a Python
  list cannot be an output argument. The bindings check what C++ only `assert`s (off in Release):
  base id `< nvars` (`IndexError`), list lengths, duplicate base ids, `inv_map`'s `dim` range and
  zero constant parts (`ValueError`), and that every vector of a call lives in one env
  (`EnvError`; the engine indexes one pool with all slots). An output list that shares a vector
  with an input list raises `ValueError` (the engine resets outputs before reading inputs).
  Tests: `python/tests/test_algorithms.py`.

**T2.3 Batch evaluation.**
- Bind `evaluate_map(m: NDAList, pts: ndarray[float64, (N, nvars)]) -> ndarray (N, len(m))` and
  `NDA.__call__(pt: Sequence[float]) -> float`. Both run one C++ loop over the points, calling
  `da_composition(ivecs, point, out)`.
- Acceptance: results equal a Python loop over `da_composition` to 1e-15; for N = 10⁵ it is at
  least 20× faster than the Python loop.
- *As built:* `NDA.__call__` calls `ad_composition`. `evaluate_map` does **not**: a loop calling
  it (which rebuilds its power tables with heap allocations on every call) was only ~7× faster
  than the Python loop `[compose(m, p) for p in pts.tolist()]`. It evaluates the same products in
  the same summation order as the double form of `ad_composition`, in blocks of 128 points with
  the point index innermost, so the compiler vectorizes over points; on x86-64 Linux with GCC the
  kernel is also cloned for AVX2 (`target_clones`, no FMA, so rounding is unchanged). Results are
  bitwise equal to `da_composition` at finite points; ~31 ns per point at (3,4) with a 3-vector
  map, about 30× faster than that Python loop (the test asserts ≥ 20×).

**T2.4 Exponent table.**
- Bind `exponents(env=None) -> ndarray[int32, (full_length, nvars)]` built from `orders_ref`.
- Acceptance: row `i` equals `index_element(i)[0]` for random `i`.
- *As built:* `exponents()` covers the current env, or the `env=` one (added with T7.3, like
  the other `env=` keywords). Rows cover `full_length()` at the current order. The test uses a dense vector,
  since `index_element(i)` returns zero exponents past the vector's length.

**T2.5 Python benchmark and gate.**
- Files: `python/bench/bench_ops.py` (the T0.5 cases, measured as in T0.5: `time.thread_time_ns`, interleaved rounds, minimum sample). It
  prints a table with columns: case, C++ ns, Python ns, overhead ns, ratio.
- Gate (must pass before leaving Stage 2). Each case gets exactly one rule, chosen by its C++ time:
  1. C++ time < 1 µs: overhead (Python ns − C++ ns) ≤ 150 ns for an op that returns a new object,
     and ≤ 100 ns for an in-place op.
  2. C++ time ≥ 1 µs: ratio (Python / C++) ≤ 1.10. No absolute overhead limit.
- The C++ numbers come from a **fresh** `bench_cpp` run made by `bench_ops.py` itself, interleaved
  with the Python runs under the same pinning, not from the stored `baseline_cpp.json`. (The stored
  file stays as a record of T0.5; machine load and clock drift made it 3–12% off within one day.)
- (Revised twice on 2026-09-29, both approved by the user. Version 1 required ratio ≤ 1.10 for every
  (6,6) case, which allows only 19 ns of overhead for `a += b` (188 ns in C++), less than one
  Python-level call (≈ 50 ns in place, ≈ 90–120 ns for a new object). Version 2 applied the absolute
  limit to every case, which cannot be measured on ops that take 10 µs to 11 s.)
- Also try `LTO` on `nanobind_add_module` together with `INTERPROCEDURAL_OPTIMIZATION` on
  `daStatic`; keep it only if the median improves by ≥ 3%.
- *As built:* `python/bench/bench_ops.py [--rounds N] [--cpp PATH]` runs the numeric T0.5 cases
  (the symbolic ones need SDA, Stage 5), pinned as the C++ run (`taskset -c 5`), with THP off via
  `prctl`. The data is dense random (`from_coeffs`), `c = a + b` etc. rebind a new object each rep.
  The C++ side is a child `bench_cpp --numeric --driven` (default path
  `build-bench/python/bench/bench_cpp`, built as in T0.5), which inherits the pinning and runs one
  round of its own, with its warm-up kept, for each line on stdin; `bench_ops.py` sends one line
  before each Python round, so the two alternate round by round, and reads the child's JSON at the
  end. `--numeric` and `--driven` were added to `bench_cpp` for this; its default run is unchanged.
  `baseline_cpp.json` is no longer read. Each case gets one rule: C++ < 1 µs, overhead ≤ 150 ns
  (`add`, `mul`, `mul_const`, `exp` return a new object) or ≤ 100 ns (`iadd` and
  `compose(m, n, out)` count as in-place); C++ ≥ 1 µs, ratio ≤ 1.10. The table's gate column
  names the rule applied. It exits 1 when the gate fails.
  Under version 2 of the gate (absolute limit on every case) the run of 2026-09-29 failed on 8 cases, all with C++ ≥ 1 µs, only
  because the absolute limit also applied to them (e.g. `n6o10/exp`: +209 µs against the stored
  baseline, +9 µs against a fresh C++ run).
  Result under the current gate (second revision, 2026-09-29, `taskset -c 5`, load 2.1 at start, 3.2 at end from a
  desktop browser on other cores, 15 rounds): **gate PASS**, all 18 cases. Sub-µs cases: overhead
  50–89 ns (worst `n3o4/add` and `n3o4/mul` +89 ns against 150; `n6o6/iadd` +54 ns and
  `n3o4/iadd` +51 ns against 100). Cases ≥ 1 µs: ratios 0.86–1.03 (worst `n6o6/mul` 1.031).
  LTO (`LTO` on the module plus `INTERPROCEDURAL_OPTIMIZATION` on `daStatic`) was tried after
  T2.6 step 3: median over cases 0.993× (slightly slower), so it was dropped.

**T2.6 Optimize (only if T2.5 fails).** In this order, re-running T2.5 after each step:
1. Profile with `perf record` on the bench script; attach the top 10 symbols to the report.
2. Reorder overloads so the most common signature comes first.
3. Construct results directly in the nanobind instance (`nb::inst_alloc`, placement-new copy of the
   left operand, compound op, `nb::inst_mark_ready`) to skip the move constructor's extra slot
   assign/free.
Stop and report if the gate still fails.
- *As built (2026-09-29).* Run while the gate was at version 1, which these steps did not meet;
  the gate was then revised and passes (T2.5):
  1. `perf` is not installed on the development machine; no profile. The cost was analysed
     instead: `Pool::free` and `Pool::alloc` zero the whole slot, so the move ctor's extra
     assign/free cost a full-length memset per result (64 KB at (6,10)).
  2. Nothing to reorder: every operator already has DA⊕DA first, then DA⊕float, and each
     benchmark case hits its first matching overload.
  3. `make_da<DA>(f)` in `arith.h`: `nb::inst_alloc`, placement-new of `DA(f())` (the prvalue
     initializes the instance directly, guaranteed elision, rather than a copy of the left
     operand plus a compound op, which would add a copy pass), `nb::inst_mark_ready`. Used by
     the binary, unary and power operators and the math functions. Overhead at (6,10) `add`
     went from +1483 ns to −640 ns, at (6,6) `add` from +283 to +118 ns.
  The failures left are all at (6,6): `add` 1.20, `iadd` 1.27, `mul_const` 1.24. Their overhead
  (50–120 ns) equals the (3,4) overhead, i.e. the fixed cost of one Python-level call (loop
  iteration, nanobind dispatch, `EnvGuard`, object alloc/free); `iadd` (188 ns in C++) would need
  ≤ 19 ns. Resolved by revising the gate (T2.5), not by further optimization.

### Stage 3 — CNDA (`std::complex<NDA>`)

**T3.1 Class `CNDA`.** Ctors: `CNDA()`, `CNDA(re: NDA, im: NDA | None = None)`,
`CNDA(z: complex, env=None)`. Properties `real` and `imag` (get returns a copy; set assigns).
`conj()` only if the C++ side has it (check `da.h`; skip if not). `__repr__`, `__str__`.

- *As built:* `CNDA(z)` also takes a float. The `env=` keyword is left to T7.3, as for `NDA`.
  `CNDA(re, im)` and the `real`/`imag` setters raise `EnvError` when the parts live in different
  envs. `conj()` is not bound: `da.h` has none. `copy()`, `__copy__` and `__deepcopy__` are bound
  as for `NDA`. `__repr__` gives `CNDA(order=4, nvars=3, nonzero=(re, im))` (non-zero terms of each
  part); `__str__` is the C++ `operator<<` text. Tests: `python/tests/test_cnda.py`.

**T3.2 Arithmetic.** Reuse `arith.h` with the scalar set {CNDA, NDA, complex, float}, both sides,
and in-place. Bind only the combinations that the C++ operators support (`da.h:98-230`); list them
in a comment.

- *As built:* `da.h` defines CNDA⊕CNDA, CNDA⊕complex, CNDA⊕float (both sides) and NDA⊕complex
  (both sides, returning a CNDA), and these are bound; it has no CNDA⊕NDA, so `c + x` raises
  `TypeError` (write `c + CNDA(x)`). `arith.h` now takes the scalar type as a template argument
  (`binops<S>`, `inplace<S>`, `unops`) and deduces the result type, and `env_of`/`same_env` for
  NDA and CNDA moved to `common.h`, with the `NB_MAKE_OPAQUE` declarations (for `CNDA` it also keeps
  nanobind's `std::complex` caster off `std::complex<NDA>`). `da.h` has no compound operators for
  CNDA: `+=`/`-=` update the parts in place, `*=`/`/=` assign the result of `*`/`/`; the object
  keeps its `id`. `__neg__`, `__pos__`, `__pow__` (int, float) as for NDA. Division by a float
  below `DBL_MIN` raises `ValueError` as for NDA. `bind_nda` returns the `NDA` class so that
  `bind_cnda` can add the NDA⊕complex operators after NDA's own. The CNDA overloads are bound in
  the A.5 rule 1 order CNDA, float, complex (operators and in-place), so a Python `int` reaches the
  float overload, not the slower complex one.

**T3.3 Functions.** The same names as T1.5 for CNDA, added as overloads to the same Python names
(NDA overload first).

- *As built:* C++ has complex versions only of `sqrt exp log asin acos atan asinh acosh atanh
  pow abs` (`functions.h:55-68`); these are bound after the NDA overloads. `sin cos tan sinh cosh
  tanh erf` of a CNDA raise `TypeError`. `abs(CNDA)` is the C++ one: the larger of the two parts'
  `norm()`, not a modulus.

**T3.4 Algorithms.** `CNDAList`; `cd_composition` (three overloads, `da.h:438-452`);
`da_composition` with complex points; `read_cd_from_file`, `compare_cd_vectors`,
`compare_cd_with_file`.
- Acceptance: tests against `cd_calculation_*.txt`, `cd_composition_*.txt` and
  `da_composition_cd_*.txt`; the T2.5 gate applied to `CNDA*CNDA` and `exp(CNDA)`.
- *As built:* `CNDAList` behaves like `NDAList` (T2.1): plain lists convert, output arguments
  are `.noconvert()`. `cd_composition` has the three C++ overloads under that one name (no alias).
  "`da_composition` with complex points" is the T2.2 overload (NDA map at a list of complex),
  already bound; a test checks it against `cd_composition` with constant CNDA points. The bindings
  check what C++ only asserts: `len(v) == nvars()`, `len(ivecs) == len(ovecs)` (`ValueError`), one
  env for every vector (`EnvError`). Only the NDA-map form reads `v` after resetting `ovecs`, so
  only there an output that is also in `v` raises `ValueError`; the CNDA-map forms work on copies
  and allow it. Overload resolution tries each list type in turn for a plain list, and nanobind
  printed a warning to stderr for each failed try, so `module.cpp` turns those warnings off
  (`nb::set_implicit_cast_warnings(false)`).
  Benchmark: `bench_cpp` and `bench_ops.py` gained the cases `cmul` (`c = a * b`, CNDA) and `cexp`
  (`c = exp(a)`), under the T2.5 gate (new object). Run of 2026-09-29
  (`taskset -c 5`, 15 rounds, load 2.8 at start, 1.6 at end): **gate PASS**, all 24 cases.
  CNDA cases: `n3o4/cmul` +122 ns against 150 (C++ 831 ns); `n3o4/cexp` 1.076, `n6o6/cmul`
  1.021, `n6o6/cexp` 1.046, `n6o10/cmul` 1.002, `n6o10/cexp` 0.994 (ratio ≤ 1.10).

### Stage 4 — `Expr` (symbolic scalar)

All symbolic stages compile only under `DA_WITH_SYMBOLIC`.

**T4.1 Class `Expr`** bound to `SymEngine::Expression`.
- Ctors from `int`, `float` and `str` (via `SymEngine::parse`). `symbols("a b c") -> tuple[Expr, ...]`
  (also accepts commas).
- Arithmetic with Expr/int/float on both sides, `__pow__`, `__neg__`.
- `__str__`; `__repr__` (`Expr('a + 2*b')`); `__float__` (via `eval_double`; raises `TypeError` if
  free symbols remain); `__eq__` (structural); `__hash__` (SymEngine `hash()`).
- `subs(dict[Expr, Expr|float])`, `expand()`, `diff(sym)`, `free_symbols() -> set[Expr]`,
  `is_zero()` (`symbolic_ops.h:27`), `simplify()` (wraps `simplified_expr`, returns a new Expr).
- *As built:* `Expr(int)` takes any Python int (beyond `long long` it goes through the parser).
  A parse error raises `ValueError`. Arithmetic is bound in the order Expr, int, float (A.5 rule 1),
  both sides, plus `__truediv__`, `__pos__` and `__rpow__`. `__eq__`/`__ne__` also compare with an
  int or float (`Expr(3) == 3`); any other type gives `False`. `subs` takes a dict whose keys and
  values are Expr, int or float. `diff` of a non-symbol raises `ValueError`. `symbols` splits on
  whitespace and commas. Tests: `python/tests/test_expr.py`.

**T4.2 sympy bridge** in `python/miradac/_sympy.py`: `to_sympy(e)` and `from_sympy(e)` via `str`;
import sympy lazily; attach them as `Expr.to_sympy` and `Expr.from_sympy`.
- Acceptance: a test for each method; `from_sympy(to_sympy(e)) == e` for 20 random expressions.
- *As built:* `to_sympy` and `from_sympy` are also module functions (`miradac.to_sympy`,
  `miradac.from_sympy`); `Expr.from_sympy` is a static method. The 20 random expressions use
  symbols, integers, rationals, `+ - * /`, integer powers, `sin` and `exp`. Floats are left out:
  sympy prints a float with 15 significant digits, so a double does not always survive the trip.

**T4.3 symengine.py interop checks.**
- Files: `python/src/bind_expr.cpp`, `python/src/se_bridge.h` (new), `python/miradac/__init__.py`.
- Do: implement the four "first interop call" checks from A.8 and the `import miradac` check
  (`dladdr` + realpath + sha256 against the constants embedded at configure time). Expose
  `symengine_interop_status() -> dict` with keys `mode` (`"shared"` or `"string"`), `symengine_version`,
  `wrapper_path`, `wrapper_needed`, `loaded_libsymengine` (list of realpaths from `/proc/self/maps`),
  `expected_libsymengine`, `layout_selftest` (bool), and `reason` (first failed check, or `None`).
  Checks run once, lazily, on the first interop call; the result is cached.
- Acceptance: with the source-built symengine.py, `mode == "shared"` and every check passes.
  Checks run once, lazily, on the first interop call; the result is cached.
- Acceptance: with the source-built symengine.py, `mode == "shared"` and every check passes.
- *As built:* `python/CMakeLists.txt` embeds four constants in `_core`: the realpath of the
  `symengine` target's library (`DA_SYMENGINE_LIB`, e.g. `<prefix>/lib/libsymengine.so.0.14.0`),
  `DA_SYMENGINE_LIB_SHA256` (from the stamp), the soname `libsymengine.so.<major>.<minor>`, and the
  pinned symengine.py version. The import check is split: `_core._libsymengine_path()` does the
  `dladdr` (on `SymEngine::eval_double`), and `__init__.py` compares the realpath and computes the
  sha256 with `hashlib` (C++ has no hash library here); the sha256 is skipped when the stamp gave
  none (`DA_IGNORE_SYMENGINE_PIN` builds). Check 2 parses the wrapper's ELF `.dynamic` section
  (`<elf.h>`); `wrapper_needed` lists all its `DT_NEEDED` entries. Check 3 collects every mapped
  file whose name starts with `libsymengine`. The layout self-test (check 4) runs only when checks
  1–3 pass (calling a wheel's `c2py` with our RCP would itself be undefined behavior); it also
  checks that `c2py` took one reference. If symengine.py is not importable, `mode` is `"string"`,
  `symengine_version` and `wrapper_path` are `None`, and `reason` says so. Tests:
  `python/tests/test_symengine_interop.py` (the import check is tested by preloading
  `~/.local/lib/libsymengine.so.0.15*` when present: `import miradac` raises `ImportError`).

**T4.4 Zero-copy conversion.**
- Files: `python/src/se_bridge.h`, `python/src/bind_expr.cpp`.
- Do (A.9): `Expr.to_symengine()` and `Expr.from_symengine(obj)`. In `"shared"` mode, pass the RCP
  (capsule `c2py` one way, read at offset `sizeof(PyObject)` the other way, after an `isinstance`
  check against `symengine.Basic`). In `"string"` mode, go through `str()` and `SymEngine::parse`,
  and warn once (`RuntimeWarning` naming the failed check and TP.2 as the fix); raise instead when
  `MIRADAC_REQUIRE_SHARED_SYMENGINE=1`. Add a private `_rcp_address(obj)` for tests.
- Acceptance tests:
  1. Shared mode: `from_symengine(x)` and back keeps the same C++ address (`_rcp_address`);
     1000 round trips of `expand((a+b+c)**20)` leave refcounts balanced (address stays alive, no leak
     under the ASan build).
  2. Shared mode is at least 100× faster than string mode on `expand((a+b+c)**20)`.
  3. In a second venv with the PyPI wheel (`uv pip install symengine==0.14.1`), `mode == "string"`,
     conversions still give equal expressions, a `RuntimeWarning` is emitted once, and with
     `MIRADAC_REQUIRE_SHARED_SYMENGINE=1` the call raises `RuntimeError`.
     conversions still give equal expressions, a `RuntimeWarning` is emitted once, and with
     `MIRADAC_REQUIRE_SHARED_SYMENGINE=1` the call raises `RuntimeError`.
- *As built:* `from_symengine` of a non-`symengine.Basic` raises `TypeError` in both modes.
  `MIRADAC_REQUIRE_SHARED_SYMENGINE` is read on every string-mode call. String mode uses
  `symengine.sympify(str(e))` one way and `SymEngine::parse(str(obj))` the other. Besides
  `_rcp_address(obj)` there is a private `_rcp_use_count(obj)` (the SymEngine refcount), which
  acceptance 1 uses to check that the refcount is balanced after 1000 round trips; the ASan part
  was run on an ASan-instrumented wheel in a scratch dir (`PYTHONMALLOC=malloc`,
  `LD_PRELOAD=libasan.so`): the leak summary is identical, byte for byte, for 0 and 1000 round
  trips (the remaining reports are interpreter and one-time module allocations). Acceptance 2 is a
  pytest comparing shared mode with the string-mode operations (`str` plus parse); measured
  ~270 ns against ~5.4 ms. Acceptance 3 is a manual check (it needs network access and its own
  venv): build the wheel with `pip wheel --no-build-isolation -Cbuild-dir=<tmp>`, install it with
  `symengine==0.14.1` from PyPI into a throwaway venv, and run the string-mode checks there.

### Stage 5 — SDA

**T5.1 Class `SDA`** (`DAVector<SymEngine::Expression>`).
- Ctors: `SDA()`, `SDA(x: Expr | float | int, env=None)`. Module functions
  `promote(nda) -> SDA` and `svar(i, env=None)` (= `promote(var(i))`).
- Methods mirroring NDA where they compile for `T = Expression`: `con` (→ Expr), `length`,
  `n_element`, `element`, `set_element(exps, Expr|float)`, `index_element`, `iszero`, `reset`,
  `coeffs() -> list[Expr]`. Try each NDA method and bind the ones that compile; list the skipped
  ones in a comment.
- *As built:* `python/src/bind_sda.cpp`. The `env=` keyword is left to T7.3, as for `NDA`.
  `SDA(x)` takes a float, an int (exact integer: `SDA(2).con` is `Expr('2')`), an `Expr` or a
  `symengine.Basic`; anything else raises `TypeError`. `con` (property; the setter resets the vector
  to that constant, as for NDA), `length`, `n_element`, `nvars`, `order` are properties;
  `element(i)` and `element(exps)` return an `Expr`; `index_element(i) -> (list[int], Expr)`;
  `copy()`, `__copy__`, `__deepcopy__`, `iszero()`, `clean()`, `reset()`, `coeffs()`, `__repr__`
  (`SDA(order=4, nvars=3, nonzero=5)`), `__str__` (C++ `operator<<`). `DAVector<Expression>` is
  explicitly instantiated, so every NDA method compiles; not bound (listed in the file's header
  comment): `norm()` and `weighted_norm()` (the kernels count a symbolic coefficient as 0),
  `from_coeffs` (use `promote(NDA.from_coeffs(a))`), `__call__` (use `evaluate`), and the `eps`
  argument of `iszero`/`clean` (the symbolic kernels ignore it and use `is_zero`). `svar(i)`
  raises `IndexError` like `var(i)`. `env_of`/`same_env`/`check_base`/`check_disjoint`/`exponents`
  moved to `common.h` as templates over `DAVector<T>`; `from_int` moved to `se_bridge.h`.

**T5.2 Arithmetic.** Scalar set {SDA, NDA (mixed ops from `interop.h:89-116`), Expr, float, int},
both sides, and in-place.
- *As built:* operand order SDA, float, int, Expr, NDA, `symengine.Basic` (A.5 rule 1), both
  sides and in place (`+= -= *= /=` keep `id`). An int becomes an exact `Expr` integer, not a
  float. C++ has no compound SDA⊕NDA operators, so the in-place forms promote the NDA first.
  `interop.h` also defines NDA⊕Expr (both sides, giving an SDA); these are bound on `NDA`, so
  `a + var(0)*var(0)` works as in `test_symbolic.cc`. Division by `0`, `0.0` or `Expr(0)` raises
  `ValueError`. `__pow__` takes int or float. A cleared env raises `EnvError`.

**T5.3 Functions and algorithms.** Math functions for SDA (`functions.h:84+`); `da_der`, `da_int`,
`da_substitute_const`, `da_substitute` and `da_composition` for SDA (`da.h:476-600`); `SDAList`.
- *As built:* C++ has SDA versions of `sqrt exp log sin cos tan asin acos atan sinh cosh tanh
  erf pow`, bound after the NDA and CNDA overloads; `asinh acosh atanh abs` of an SDA raise
  `TypeError`. `SDAList` behaves like `NDAList` (plain lists convert, outputs `.noconvert()`).
  `da_der`/`da_int` in both forms; `da_substitute_const`; `da_substitute` with an SDA, a float
  (C++ has only `da_substitute_const` for that, which it calls), several bases, and lists;
  `da_composition` with SDA points (output list) and with float points (returns `list[Expr]`).
  The checks are those of T2.2. nanobind replaces a function object when an overload is added,
  so `bind_sda` sets the `der`/`int_`/`substitute`/`compose` aliases again.

**T5.4 Evaluation.** `evaluate(sda, values: dict[Expr, float]) -> NDA`, and
`evaluate(sda, syms: Sequence[Expr], vals: Sequence[float]) -> NDA`. Implement both with the
vector overload (`interop.h:175`), which is faster; for a dict, split it into keys and values.
- *As built:* keys and `syms` entries may be `Expr` or `symengine.Basic`, values any real number.
  A coefficient with a symbol that has no value raises `ValueError` (SymEngine's "Symbol not in
  the symbols vector"), as do `syms`/`vals` of different lengths.

**T5.5 symengine.py objects in SDA.** Every SDA entry point that takes an `Expr` also takes a
`symengine.Basic` (converted via T4.4): the `SDA(...)` constructor, `set_element`, arithmetic
scalars, and `evaluate` dict keys. `coeffs(as_symengine=True)` and `con_symengine()` return
symengine.py objects. Tests run in shared mode and check `_rcp_address` identity for coefficients.
- *As built:* also `con` (setter), `subs` keys and values, and the `evaluate` sequence form.
  `se_bridge::is_symengine(obj)` tests `isinstance(obj, symengine.Basic)` without importing
  symengine.py or running the A.8 checks (an object cannot be one if the module was never
  imported); the arithmetic overload for it goes to the next overload otherwise, so other types
  still give `TypeError`. Tests: `test_symengine_objects` in `python/tests/test_sda.py`
  (`_rcp_address` identity for `SDA(x)`, `set_element`, `con`, `coeffs(as_symengine=True)`,
  `con_symengine()`).

**T5.6 Per-coefficient helpers.** `SDA.simplify()`, `SDA.expand()` and `SDA.subs(dict)`, each
returning a new SDA (loop over the slot's coefficients in C++).
- Acceptance: port the cases in `test/test_symbolic.cc` to pytest (same inputs, same expected
  values); `evaluate(exp(a*svar(0) + 1), {a: 0.3})` equals `exp(0.3*var(0) + 1)` to 1e-13;
  `evaluate(promote(x), {}) == x`. Benchmark: Python `exp(SDA)` at (3,3) takes at most 1.05× the
  C++ time (add this case to T0.5 and T2.5).
- *As built:* `simplify()` skips coefficients that are zero (as `DAVector::simplify` does);
  `subs` takes Expr, int, float or `symengine.Basic` keys and values. The copy is built in the
  new instance (`make_da`), then its slot is rewritten, so the source is unchanged. Tests:
  `python/tests/test_symbolic.py` (the port: every case of `test_symbolic.cc` except the
  non-default-env `erf` case, moved to T7.4 since it needs `Env`) and `python/tests/test_sda.py`
  (T5.1-T5.6). T0.5 already timed `sym_n3o3/exp`; `bench_ops.py` now runs the symbolic cases too
  (`bench_cpp --driven`, no longer `--numeric`) with ratio ≤ 1.05 for `sym_n3o3/exp` and the T2.5
  rules for `sym_n3o3/mul`. Run of 2026-09-29 (`taskset -c 5`, 15 rounds, load 2.1 at end):
  **gate PASS**, all 26 cases; `sym_n3o3/exp` 1.014 (73 773 vs 74 783 ns), `sym_n3o3/mul` 1.016.
  A run just before, under load 6–8 (a video call), failed only on `n3o4/cexp` (1.106, a CNDA
  case this stage does not touch; 1.091 in the rerun).

### Stage 6 — CSDA (`std::complex<SDA>`)

**T6.1 Class `CSDA`.** Ctor `CSDA(re: SDA, im: SDA | None = None)`; `promote(cnda) -> CSDA`
(`interop.h:317`); `real` and `imag` properties.
- *As built:* `python/src/bind_csda.cpp`. Also `CSDA()`; `im` defaults to an exact zero SDA.
  `CSDA(re, im)` and the `real`/`imag` setters raise `EnvError` for parts in different envs; the
  getters return copies. `copy()`, `__copy__`, `__deepcopy__` and `__repr__`
  (`CSDA(order=3, nvars=2, nonzero=(2, 2))`) as for CNDA. `__str__` was missing at first because
  C++ had `operator<<` for `std::complex<DAVector<T>>` only for `T = double`; since 2026-09-30 C++
  prints a CSDA as the real part's table, then the imaginary part's, and `__str__` uses it.

**T6.2 Arithmetic and functions.** Scalar set {CSDA, SDA, complex, float, Expr}, as supported by the
`is_da_coeff` templates in `da.h:98-230` (bind only what compiles; list the skipped combinations);
math functions; `abs` (`functions.h:73`). Accept `symengine.Basic` wherever `Expr` is accepted (as in T5.5).
- *As built:* the `is_da_coeff` templates define CSDA⊕CSDA, CSDA⊕complex, CSDA⊕float (both
  sides) and SDA⊕complex (both sides, giving a CSDA); these are bound. Since 2026-09-30 also:
  CSDA⊕SDA and CNDA⊕NDA (both sides; `std::complex`'s own templates already provided them, they
  had only been left unbound), and CSDA⊕Expr / CSDA⊕`symengine.Basic` (both sides; C++ operators
  added in `da.h`). CSDA⊕NDA and CSDA⊕CNDA still raise `TypeError`. In-place operators, `__neg__`, `__pos__`, `__pow__` and division by a float below
  `DBL_MIN` behave as for CNDA; the CNDA operator block moved into `arith.h` as
  `bind_complex_arith<C>`, used by both. SDA⊕complex is bound after SDA's own operand types, so
  `s + 2` stays an SDA. Functions: `sqrt exp log asin acos atan asinh acosh atanh pow` (int, float)
  and `abs` (an SDA, `sqrt(re^2 + im^2)`), added after the other overloads; `sin cos tan sinh
  cosh tanh erf` of a CSDA raise `TypeError`. The complex functions and the three
  `cd_composition` forms of CNDA and CSDA come from one template, `bind_complex_functions<T>`
  (`arith.h`). `symengine.Basic` enters CSDA only through `evaluate` keys (T6.3), since no CSDA
  entry point takes an `Expr`.

**T6.3 Algorithms.** `CSDAList`; `cd_composition` for `T = Expression`;
`evaluate(csda, …) -> CNDA` (`interop.h:337,347`).
- Acceptance: port `test/test_symbolic_cd.cc` and the checks in
  `examples/example_complex_symbolic.cc` (for example, the constant parts of `exp(z)` equal
  `exp(1.5)*cos(0.5)` and `exp(1.5)*sin(0.5)`).
- *As built:* `CSDAList` behaves like `SDAList` (plain lists convert, outputs `.noconvert()`).
  `cd_composition` has the three C++ forms for `T = Expression` with the T3.4 checks.
  `evaluate(csda, dict)` and `evaluate(csda, syms, vals)` both use the vector overload
  (`interop.h:347`), with the T5.4 key types and errors; the SDA and CSDA forms share
  `se_bridge::bind_evaluate<S>`, and `to_expr` moved from `bind_sda.cpp` to `se_bridge.h`.
  Tests: `python/tests/test_csda.py` (every case of `test_symbolic_cd.cc`, the checks of
  `example_complex_symbolic.cc`, and T6.1-T6.3).

### Stage 7 — Multi-env API

**T7.1 [C++] Per-env base in `inv_map` and `erf(SDA)`.**
- **Status: done (2026-09-29).** Implemented with a new public `da::da_base(i)` (`base.h`), which
  builds the i-th base vector in the current env; Python's `var(i)` (T1.3) should call it. The same
  change fixed two order bugs: `inv_map` and multi-base `da_substitute` (`engine.cpp`) called
  `restore_order()` and so reset an order the caller had lowered; both now return to the caller's
  order. Regression tests: `[inv_map]`, `[numeric_order]`, `[symbolic][multienv]`.
- Files: `src/base.cpp` (`inv_map`, lines 191 and 203), `src/functions.cpp` (`erf(const SDA&)`,
  line 790), `test/test_multienv.cc`.
- Do: replace each use of the global `da::base[j]` with a base vector built locally, in the input
  vector's env, via `detail::ad_var`. Build each needed base vector once per call, not inside loops.
- Acceptance: new C++ tests show that `inv_map` of the same map, and `erf` of the same SDA, give
  identical coefficients in a `da_make_env` env and in the default env; both also work after the
  default env was cleared with `da_clear()`. Full `ctest` green.

**T7.2 [C++] Cross-env copy helpers.**
- Files: `include/da/env.h` (or `davector.h`), `test/test_multienv.cc`.
- Do: add `template<class T> DAVector<T> import_to(DAEnv& dst, const DAVector<T>& src);`, which
  wraps `DAEnv::import<T>` and returns a vector that owns the new slot. Add
  `promote_to(DAEnv& dst, const NDA& src) -> SDA` the same way, wrapping `DAEnv::promote`.
- Acceptance: C++ tests show a round trip keeps the data equal, and different layouts throw
  `std::invalid_argument`.
- *As built:* both in `davector.h` (`promote_to` under `DA_WITH_SYMBOLIC`), each one line: the
  result adopts the slot through a new constructor `DAVector(DAEnv& env, unsigned slot)` (takes
  ownership of a slot already taken from `env`'s pool), so neither needs or changes the current
  env. `DAEnv::promote` had no caller and turned every coefficient into
  `Expression(double)` (`0.0`, `1.0`), unlike `promote(const NDA&)`, which makes an
  integer value an exact integer (`0`, `1`); it now uses the same rule, so `Env.promote(x)` and
  `promote(x)` give the same coefficients. Tests: `[import_to]` in `test/test_multienv.cc`
  (suite at 917 assertions / 80 cases, ASan run green).

**T7.3 Class `Env`.**
- Files: `python/src/bind_env.cpp`, `python/miradac/__init__.py`.
- Bind `DAEnv` as `Env`, always with `nb::rv_policy::reference`. The C++ side owns every env, and
  retired shells live until exit, so the pointer never dangles. `__eq__` and `__hash__` compare the
  pointer.
- `Env(order, nvars, pool_size, table=False)`: call `da_make_env`, then restore the previous
  current env with `da_exchange_env`. Creating an env must not change the current env.
- `Env.current()` and `Env.default()` (raise `EnvError` if there is none).
- Properties: `order`, `max_order`, `nvars`, `full_length`, `pool_size`, `count`, `remain`,
  `live_slots`, `retired`.
- `select()`; `__enter__` (exchange, remember the previous env) and `__exit__` (restore it);
  re-entrant via a stack.
- `close()`: calls `da_destroy_env` (idempotent; for the default env, `da_clear`).
- `order(n)`: a context manager scoped to this env (same algorithm as T1.6, under `EnvGuard`).
- `import_(v)` for NDA/SDA/CNDA/CSDA (for CNDA/CSDA, import the real and imag parts) and
  `promote(v)`, both via T7.2.
- An `env=` keyword on all constructors and on `var`/`svar` (implemented through `EnvGuard`).
- A `.env` property on NDA/CNDA/SDA/CSDA returning the owning `Env`.
- *As built:* `python/src/bind_env.cpp`, via `nb::new_` with `rv_policy::reference`.
  **No env is ever deleted.** `da_clear()` and `da_destroy_env()` delete an env that has no
  live vector (only one with live vectors is retired), which would leave an `Env` object
  dangling, and a new env at the same address would reuse the old Python object. So
  `clear()`, `init()` (which drops the previous default env) and `Env.close()` first take
  one pool slot, never freed, from an env with no live slot; C++ then always retires it, and
  its shell (a few hundred bytes) stays until exit. For the default env they first clear the
  `da::base` vectors (which `da_clear()`/`da_init()` would free anyway before deleting), so
  the check sees only user vectors. `Env()` restores the previous current env in a
  destructor, so it stays current even when construction throws (e.g. `MemoryError`). `init()` and `Env()` reject
  `pool_size=0` (`ValueError`), since then there is no slot to take. C++ has no accessor for
  the default env, so the bindings record the env of the last `init()` until `clear()`
  (`Env.default()`). The property `order` would clash with the method `order(n)`, so it is
  **`current_order`** (the module function's name); `max_order` is the construction order.
  Properties other than `retired` raise `EnvError` on a closed env; `repr` gives
  `Env(order=4, nvars=3, pool_size=1000)` or `Env(retired)`. `Env.current()` returns a closed
  env if it is still selected (its use raises). `with` blocks nest, so one thread-local stack
  of previous envs serves every `Env`. Also bound: `change_order(n) -> bool` and
  `restore_order()` on `Env` (under `EnvGuard`); `Env.order(n)` is in `__init__.py`, and the
  module-level `order(n)` (T1.6) is now `Env.current().order(n)`, so it restores the env it
  changed. `import_` takes NDA, CNDA, SDA, CSDA; `promote` takes NDA and CNDA (the parts
  separately); both raise `ValueError` for different layouts and `EnvError` for a closed env
  on either side. `env=` is keyword-only on `NDA()`, `NDA(x)`, `NDA.from_coeffs`, `var`,
  `CNDA()`, `CNDA(z)`, `SDA()`, `SDA(x)`, `svar`, `CSDA()`; `exponents(env=None)` (T2.4) takes
  it too. Constructors with a DA argument take the env of that argument.

**T7.4 Tests.**
- Two envs with different orders; arithmetic in each truncates at the right order.
- Ops on env-B vectors while env A is current give correct results (every function in T1.5,
  T2.2 and T5.3).
- Mixing envs in one binary op raises `EnvError`.
- `close()` with live vectors: later use raises `EnvError`, and the process exits cleanly.
- `with env.order(n)` affects only that env.
- `inv_map` in a non-default env matches the default-env result.
- `erf(SDA)` in a non-default env matches the default-env result (the last case of
  `test/test_symbolic.cc`, left out of the T5.6 port because it needs `Env`).
- Benchmark: rerun T2.5; the gate must still pass.
- *As built:* `python/tests/test_multienv.py`. "Every function" is one test that runs each
  function of T1.5, T2.2 and T5.3 (plus `evaluate_map`, `NDA.__call__` and `evaluate`) on
  vectors of env B (order 4, 3 vars) once with B current and once with the default env A
  (order 2, 2 vars) current, and requires equal results that all live in B. The `close()` exit
  check runs in a subprocess. inv_map and `erf(SDA)` results are compared between the default
  env and an `Env` of the same layout (inv_map also after `clear()`).
  Benchmark (2026-09-30, `taskset -c 5`, 15 rounds): the first run (load 3.0 at start)
  failed only on `n3o4/cexp` at 1.105, the CNDA case that also failed under load in Stage 5
  (its path is unchanged here). The rerun (load 1.25 at start, 1.74 at end) gave
  **gate: PASS** on all 26 cases. Sub-µs overheads were 49–122 ns (worst `n3o4/cmul` +122
  against 150; `n3o4/iadd` +49 and `n6o6/iadd` +49 against 100). Cases ≥ 1 µs had ratios
  0.89–1.09 (worst `n3o4/cexp` 1.086), and `sym_n3o3/exp` was 1.000 against 1.05.
  Rerun after the `Env` lifetime fixes (2026-09-30, load 1.16 at start, 1.17 at end): **gate:
  PASS** on all 26 cases; sub-µs overheads 49–126 ns (worst `n3o4/cmul` +126), ratios ≥ 1 µs
  0.90–1.09 (worst `n3o4/cexp` 1.085), `sym_n3o3/exp` 1.001.

### Stage 8 — Finish

**T8.1 Stubs.** Generate `python/miradac/_core.pyi` with `python -m nanobind.stubgen`; add
`py.typed`; run `mypy --strict` on the examples.
- *As built:* the stub is generated with `-P` (so `int_`, `Env.import_` and the private test
  helpers are in it; stubgen treats a trailing underscore as private) and the pattern file
  `python/stubgen_patterns.txt`, which holds the regeneration command. The patterns add
  `__version__` and the names `__init__.py` attaches at import (`Env.order`, `Expr.to_sympy`,
  `Expr.from_sympy`), and drop the values of `_LIBSYMENGINE_EXPECTED`/`_SHA256` (a home path).
  So that the stub names real types, not `object`: `make_da` returns
  `nb::typed<nb::object, DA>` (its callers return `auto`); dicts, sets, tuples, lists and
  sequences are `nb::typed` (`dict[Any, Any]` for `subs`, `dict[Any, float]` for `evaluate`,
  `tuple[Expr, ...]` for `symbols`, `list[Any]` for `SDA.coeffs`); the symengine.py operand of
  the SDA operators is a new `se_bridge::SymengineBasic` (an `nb::object` whose check is
  `is_symengine`, named `symengine.Basic`; it replaces the `next_overload` lambda). SDA⊕complex
  moved from `bind_csda.cpp` to `bind_sda.cpp`, before the `symengine.Basic` overload, because to
  a type checker the untyped `symengine.Basic` is `Any` and would shadow `complex`
  (`bind_csda` no longer takes the SDA class). `Env.__eq__`/`__ne__` take any object
  including `None` (`nb::arg().none()`), `Expr` gained last `__eq__`/`__ne__` overloads for
  any object (unequal, as before), and
  `Env.__exit__` takes three named arguments instead of `*args`, so the stub is compatible with
  `object`. `__init__.py` and `_sympy.py` are annotated. `pyproject.toml` holds the mypy
  settings: `ignore_missing_imports` for symengine and sympy (no type information), and, for
  `miradac._core` only, `overload-cannot-match` and `misc` disabled (the stub mirrors the
  runtime overloads: int after float, the A.5 order, and in-place operators that take fewer
  types than the binary ones). Tests (`python/tests/test_stubs.py`): the committed stub equals
  a fresh stubgen run; `py.typed` exists; `mypy --strict` passes on `python/examples` and on the
  package itself (`MYPYPATH=python`, otherwise mypy finds the editable install through its
  `.pth` and silences the package as installed). mypy 2.3.1 was installed into `.venv`.

**T8.2 Examples.** Port `examples/examples.cc`, `example_interop.cc`, `example_complex_da.cc`,
`example_1_symbolic.cc` and `example_complex_symbolic.cc` to `python/examples/*.py`; each one runs
as a test.
- *As built:* same names with `.py`. Each script calls `init()`, runs `main()` (so its vectors
  are freed before `clear()`), then `clear()`. `x1 + x2 * 1i` becomes `CNDA(x1, x2)`, since there
  is no NDA⊕CNDA operator (T3.2). `example_complex_symbolic.py` exits 1 when a cross-check fails
  (the C++ one only prints FAIL). Their output matches the C++ examples' up to float formatting
  (Python prints 17 significant digits) and slot numbers. `python/tests/test_examples.py` runs
  each in a subprocess in a temporary directory and requires exit code 0 and no "FAIL".

**T8.3 Docs.** Docstrings on every binding; a "Python" section in the README with build
instructions and the A.11 example.
- *As built:* every function, class, method and property of `_core` has a docstring; for an
  overloaded name the first overload carries it (nanobind joins them). Operators get theirs from
  `op_doc(name)` in `arith.h` ("Return self + value."). The `*List` classes have a class
  docstring; their methods are nanobind's. `test_every_binding_has_a_docstring` in
  `test_stubs.py` checks this. The README example is A.11 with one change: `e.import_(x)` would
  raise `ValueError` (`e` has another layout than the default env), so it imports into a third
  env `f` of the default env's layout. The example was run and passes `mypy --strict`.

**T8.4 Performance report.** `python/bench/REPORT.md`: the final table from T2.5 plus the symbolic
cases, with machine and compiler info.
- *As built:* a fresh `bench_ops.py` run after the Stage 8 binding changes (2026-09-30,
  `taskset -c 5`, 15 rounds, load 1.31 at start, 1.52 at end): **gate: PASS** on all 26
  cases (numeric, CNDA and symbolic). Sub-µs overheads 49–118 ns (worst `n3o4/cmul` +118
  against 150); ratios ≥ 1 µs 0.79–1.09 (worst `n3o4/cexp` 1.088); `sym_n3o3/exp` 1.000. The
  report holds the table, the method, the machine and toolchain, and notes on the numbers.

### Deferred (not in this plan)

- GIL release (needs per-env locks and a thread-safe SymEngine build); free-threaded Python.
- Pickling (needs a portable format that records order and nvars).
- Fast repeated SDA evaluation through SymEngine `Lambdify` (a big win for parameter scans).
- Stable-ABI builds; wheels for other platforms (aarch64, macOS, Windows).
- *Done 2026-09-30:* portable Linux wheels. `[tool.cibuildwheel]` in `pyproject.toml` builds
  manylinux_2_28 x86_64 wheels for CPython 3.10–3.14 (no 3.9: current manylinux images dropped
  it; no free-threaded builds, see A.6). The pinned SymEngine is built in the image
  (`setup_symengine.sh --no-python`) and bundled with GMP by auditwheel; each wheel runs the
  Python suite. With `MIRADAC_BUNDLED_SYMENGINE` the import-time libsymengine check is skipped
  and symengine.py interop reports string mode with that reason (zero-copy needs one shared
  libsymengine, which a wheel's private copy can never be). `.github/workflows/wheels.yml` builds
  them on `v*` tags (attached to a GitHub Release, with the sdist) or by hand. Stored-reference
  comparisons in `test_symbolic.py` use 1e-12: the manylinux build differs by up to 1.6e-13
  relative on an order-5 `tanh` coefficient.
