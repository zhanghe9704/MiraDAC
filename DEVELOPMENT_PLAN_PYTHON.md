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
| Pool is fixed size; exhaustion throws `std::runtime_error("Pool::assign: Run out of vectors")`. | `pool.h:121-124` | Maps to a Python exception, no crash. CPython refcounting frees temporaries immediately, so pool pressure matches C++. |
| Move ctor still calls `pool.assign()` for the moved-from object. | `davector.h:84-88` | Returning results by value costs one extra assign/free pair. Measure in Stage 2; optimize only if the gate fails (T2.6). |
| Binary DA×DA ops call `check_env`, which throws `std::logic_error` on mismatch. | `env.h:262-266` | Map to `miradac.EnvError`. |
| Constructors (`NDA()`, `NDA(1.0)`) take the **current thread-local env**. Ops like `a + 2.0` use `a`'s env. Library functions (`exp`, …) may create temporaries in the current env. | `davector.h:73,93`; `env.h:177-183` | Every bound call on a DA object selects that object's env for the duration of the call (A.4). This makes multi-env correct without auditing every C++ function. |
| `da_clear()` / `da_destroy_env()` **retire** an env that still has live vectors: memory freed, shell kept until exit. | `env.h:106-122, 200+` | `DAEnv*` pointers never dangle. A vector of a retired env must not be used: every bound call checks `env_->retired()`. |
| Global `da::base` belongs to the default env only. `inv_map` and `erf(SDA)` use it internally. | `base.h:44`, `src/base.cpp:191,203`, `src/functions.cpp:790` | Python gets `var(i)`, built per env. `inv_map` and `erf(SDA)` in a non-default env need a [C++] fix (T7.1). |
| `da_change_order(n)` / `da_restore_order()` act on the current env. Restore always returns to the **original** order (not a stack). | `da.h:293-305`, `src/layout.cpp:417-435` | The Python context manager saves the previous order and restores it, giving correct nesting. |
| SymEngine is built **without** `WITH_SYMENGINE_THREAD_SAFE` (non-atomic refcounts). The pinned build (A.8) keeps this setting. | `symengine_config.h` | Never release the GIL around symbolic code. |
| MiraDAC builds and passes its full C++ suite (805 assertions / 71 cases, plus the ASan run) against SymEngine commit `153b7e98` (0.14.0 + 14 commits), not only 0.15.0. | Verified 2026-09-29 in a scratch build | The downgrade in A.8 is safe for the C++ library. |
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
| MiraDAC CMake configure (C++ and Python builds) | `find_package(SymEngine 0.14.0 EXACT CONFIG)`; the stamp exists; stamp commit equals the pin; stamp options equal the pin | Configure error naming the found prefix and the fix command |
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
- Acceptance:
  1. Configure against the pinned prefix succeeds and `ctest` passes (805 assertions).
  2. Configure against `~/.local` (0.15.0) fails with the message.
  3. Configure against a SymEngine built from tag `v0.14.0` without a stamp fails with the message.
  4. The `test/consumer` project builds against the installed MiraDAC and the pinned SymEngine.

**TP.4 Migrate this machine (needs user approval before running).**
- Do: run `scripts/setup_symengine.sh`; rebuild MiraDAC against the pin; reinstall it into
  `~/.local` (same install command as before). Leave SymEngine 0.15.0 in `~/.local` untouched; the
  EXACT check keeps MiraDAC from using it.
- Acceptance: `ldd ~/.local/lib/libdaShared.so` shows `libsymengine.so.0.14` from the pin prefix;
  the out-of-tree `find_package(da)` check from the previous install still prints the same values.

**TP.5 Documentation.**
- Files: `README.md` (new section "SymEngine version"), `cmake/symengine_pin.txt` (header comment).
- Content: why the pin exists (A.8, first paragraph); the exact commit and version; how to install
  (`scripts/setup_symengine.sh`); the warning that `pip install symengine` from PyPI installs a
  wheel that cannot share objects, and the `--no-binary symengine` flag; how to read
  `miradac.symengine_interop_status()`; the **upgrade procedure**: wait for a symengine.py release;
  take the commit from its `symengine_version.txt`; update `cmake/symengine_pin.txt` (both sha256
  values too); rerun TP.2 on a fresh prefix, the full C++ suite, the Python suite, and the T0.5
  symbolic benchmark against old and new; commit the pin change on its own.
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

**T0.2 [C++] `da_exchange_env`.**
- Files: `include/da/env.h`, `src/env.cpp`, `test/test_env.cc`.
- Do: add `DAEnv* da_exchange_env(DAEnv* env) noexcept;`. It sets the thread-local current env to
  `env` (may be null) and returns the previous pointer (may be null). Doxygen comment in the style
  of the neighboring functions.
- Acceptance: new Catch2 test: exchange to env B returns A; exchange back returns B;
  after `exchange(nullptr)`, `da_current_env()` throws. Full `ctest` green.
- (Done 2026-09-29.) Also fix the comment at `include/da/layout.h:129` to say that `restore_order()` returns to the
  original order given at construction (not to the value before the last `change_order()`), and
  the matching comment on `da_restore_order` in `da.h`.

**T0.3 Python environment.**
- Do: `uv venv .venv --python 3.13`;
  `uv pip install --python .venv nanobind scikit-build-core pytest numpy sympy`.
  Add `.venv/` and `__pycache__/` to `.gitignore`.
- Acceptance: `.venv/bin/python -c "import nanobind, numpy, sympy; print(nanobind.__version__)"`
  prints a version ≥ 2.0.

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

**T0.5 C++ baseline benchmark.**
- Files: `python/bench/bench_cpp.cpp`, `python/bench/CMakeLists.txt`; root option
  `DA_BUILD_BENCH` (default OFF) adds it.
- Do: for (nvars, order) in {(3,4), (6,6), (6,10)}, time `a+b`, `a*b`, `a+=b`, `a*2.0`, `exp(a)`,
  and `da_composition` of an nvars-map with an nvars-map. Use dense random vectors, a warm-up,
  and the median of ≥ 7 runs of ≥ 0.2 s each. Print JSON `{case: ns_per_op}`.
- Also time symbolic cases at (3,3): `SDA*SDA` and `exp(SDA)` with one symbol per coefficient.
- Acceptance: writes `python/bench/baseline_cpp.json`; two consecutive runs agree within 5%. Run
  the symbolic cases once against the pinned SymEngine and once against 0.15.0 (`~/.local`), and
  record both. If the pin is more than 10% slower on any case, stop and report to the user.

### Stage 1 — NDA core

**T1.1 Common infrastructure.**
- Files: `python/src/common.h`, `python/src/module.cpp`.
- Do: `struct EnvError : std::runtime_error`; `EnvGuard` exactly as in A.4; helper
  `da::DAEnv* current_env_or_throw()`; register `miradac.EnvError` via `nb::exception<EnvError>`
  and the translator from A.10.
- Acceptance: compiles; T1.7 tests cover each mapping.

**T1.2 Env-level functions (default env).**
- Files: `python/src/bind_env.cpp`.
- Bind: `init(order, nvars, pool_size, table=False)` → `da_init`; `clear()` → `da_clear`;
  `count()`, `remain()`, `poolsize()`, `full_length()`, `nvars()`, `max_order()` of the current env;
  `get_eps()` / `set_eps(x)`; `change_order(n) -> bool` (call `layout().change_order` directly, so
  nothing is printed to stdout); `restore_order()`; `current_order()`.
- Acceptance: tests for init/clear cycles (init, clear, init again with a different order).

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

**T1.4 Arithmetic.**
- Files: `python/src/arith.h`, `python/src/bind_nda.cpp`.
- Do: a template `bind_arith<DA>(nb::class_<DA>&)` that binds `__add__ __sub__ __mul__ __truediv__`
  for (DA, DA), then (DA, float); `__radd__ __rsub__ __rmul__ __rtruediv__` for float;
  `__iadd__ __isub__ __imul__ __itruediv__` for DA and float; `__neg__ __pos__`; `__pow__` for int,
  then float, calling `da::pow`. All with `nb::is_operator()`.
- Acceptance: tests for every operator against hand-computed coefficients; `x += y` keeps `id(x)`.

**T1.5 Math functions.**
- Bind module-level `sqrt exp log sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh erf abs pow`
  for NDA (`functions.h:29-47`), each guarded. Also bind `compare_da_with_file`.
- Acceptance: `exp`, `log`, `sqrt`, `pow(x,3)` and `pow(x,0.3)` match `test/exp_da.txt`,
  `log_da.txt`, `sqrt_da.txt`, `pow3_da.txt` and `pow0p3_da.txt` via `compare_da_with_file`.

**T1.6 `order()` context manager.**
- Files: `python/miradac/__init__.py`.
- Do: `@contextmanager def order(n)`: save `current_order()`; call `change_order(n)` and raise
  `ValueError` if it returns False; yield; on exit, if the saved order equals `max_order()` call
  `restore_order()`, else call `change_order(saved)`.
- Acceptance: nested `with order(3): with order(2): ...` restores 3, then the original order; the
  order is also restored when the block raises.

**T1.7 Error-path tests.**
- Acceptance tests: using an NDA after `clear()` raises `EnvError` and does not crash the process;
  exhausting a pool of 10 raises `RuntimeError`; `var(99)` raises `IndexError`.

### Stage 2 — NDA algorithms, bulk API, performance gate

**T2.1 Opaque lists.**
- Do: `nb::bind_vector<std::vector<NDA>>(m, "NDAList")`. Check whether nanobind's `bind_vector`
  accepts a Python `list` where `std::vector<NDA>&` is expected. If it does not, add explicit
  overloads that build a temporary `std::vector<NDA>` from a list. Record the result in a comment.
- Acceptance: tests pass both an `NDAList` and a plain list to one algorithm from T2.2.

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

**T2.3 Batch evaluation.**
- Bind `evaluate_map(m: NDAList, pts: ndarray[float64, (N, nvars)]) -> ndarray (N, len(m))` and
  `NDA.__call__(pt: Sequence[float]) -> float`. Both run one C++ loop over the points, calling
  `da_composition(ivecs, point, out)`.
- Acceptance: results equal a Python loop over `da_composition` to 1e-15; for N = 10⁵ it is at
  least 20× faster than the Python loop.

**T2.4 Exponent table.**
- Bind `exponents(env=None) -> ndarray[int32, (full_length, nvars)]` built from `orders_ref`.
- Acceptance: row `i` equals `index_element(i)[0]` for random `i`.

**T2.5 Python benchmark and gate.**
- Files: `python/bench/bench_ops.py` (the T0.5 cases, `time.perf_counter_ns`, median of 7). It
  prints a table with columns: case, C++ ns, Python ns, overhead ns, ratio.
- Gate (must pass before leaving Stage 2): at (3,4), overhead ≤ 150 ns per binary op and ≤ 100 ns
  per in-place op; at (6,6) and (6,10), ratio ≤ 1.10.
- Also try `LTO` on `nanobind_add_module` together with `INTERPROCEDURAL_OPTIMIZATION` on
  `daStatic`; keep it only if the median improves by ≥ 3%.

**T2.6 Optimize (only if T2.5 fails).** In this order, re-running T2.5 after each step:
1. Profile with `perf record` on the bench script; attach the top 10 symbols to the report.
2. Reorder overloads so the most common signature comes first.
3. Construct results directly in the nanobind instance (`nb::inst_alloc`, placement-new copy of the
   left operand, compound op, `nb::inst_mark_ready`) to skip the move constructor's extra slot
   assign/free.
Stop and report if the gate still fails.

### Stage 3 — CNDA (`std::complex<NDA>`)

**T3.1 Class `CNDA`.** Ctors: `CNDA()`, `CNDA(re: NDA, im: NDA | None = None)`,
`CNDA(z: complex, env=None)`. Properties `real` and `imag` (get returns a copy; set assigns).
`conj()` only if the C++ side has it (check `da.h`; skip if not). `__repr__`, `__str__`.

**T3.2 Arithmetic.** Reuse `arith.h` with the scalar set {CNDA, NDA, complex, float}, both sides,
and in-place. Bind only the combinations that the C++ operators support (`da.h:98-230`); list them
in a comment.

**T3.3 Functions.** The same names as T1.5 for CNDA, added as overloads to the same Python names
(NDA overload first).

**T3.4 Algorithms.** `CNDAList`; `cd_composition` (three overloads, `da.h:438-452`);
`da_composition` with complex points; `read_cd_from_file`, `compare_cd_vectors`,
`compare_cd_with_file`.
- Acceptance: tests against `cd_calculation_*.txt`, `cd_composition_*.txt` and
  `da_composition_cd_*.txt`; the T2.5 gate applied to `CNDA*CNDA` and `exp(CNDA)`.

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

**T4.2 sympy bridge** in `python/miradac/_sympy.py`: `to_sympy(e)` and `from_sympy(e)` via `str`;
import sympy lazily; attach them as `Expr.to_sympy` and `Expr.from_sympy`.
- Acceptance: a test for each method; `from_sympy(to_sympy(e)) == e` for 20 random expressions.

**T4.3 symengine.py interop checks.**
- Files: `python/src/bind_expr.cpp`, `python/src/se_bridge.h` (new), `python/miradac/__init__.py`.
- Do: implement the four "first interop call" checks from A.8 and the `import miradac` check
  (`dladdr` + realpath + sha256 against the constants embedded at configure time). Expose
  `symengine_interop_status() -> dict` with keys `mode` (`"shared"` or `"string"`), `symengine_version`,
  `wrapper_path`, `wrapper_needed`, `loaded_libsymengine` (list of realpaths from `/proc/self/maps`),
  `expected_libsymengine`, `layout_selftest` (bool), and `reason` (first failed check, or `None`).
  Checks run once, lazily, on the first interop call; the result is cached.
- Acceptance: with the source-built symengine.py, `mode == "shared"` and every check passes.

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

### Stage 5 — SDA

**T5.1 Class `SDA`** (`DAVector<SymEngine::Expression>`).
- Ctors: `SDA()`, `SDA(x: Expr | float | int, env=None)`. Module functions
  `promote(nda) -> SDA` and `svar(i, env=None)` (= `promote(var(i))`).
- Methods mirroring NDA where they compile for `T = Expression`: `con` (→ Expr), `length`,
  `n_element`, `element`, `set_element(exps, Expr|float)`, `index_element`, `iszero`, `reset`,
  `coeffs() -> list[Expr]`. Try each NDA method and bind the ones that compile; list the skipped
  ones in a comment.

**T5.2 Arithmetic.** Scalar set {SDA, NDA (mixed ops from `interop.h:89-116`), Expr, float, int},
both sides, and in-place.

**T5.3 Functions and algorithms.** Math functions for SDA (`functions.h:84+`); `da_der`, `da_int`,
`da_substitute_const`, `da_substitute` and `da_composition` for SDA (`da.h:476-600`); `SDAList`.

**T5.4 Evaluation.** `evaluate(sda, values: dict[Expr, float]) -> NDA`, and
`evaluate(sda, syms: Sequence[Expr], vals: Sequence[float]) -> NDA`. Implement both with the
vector overload (`interop.h:175`), which is faster; for a dict, split it into keys and values.

**T5.5 symengine.py objects in SDA.** Every SDA entry point that takes an `Expr` also takes a
`symengine.Basic` (converted via T4.4): the `SDA(...)` constructor, `set_element`, arithmetic
scalars, and `evaluate` dict keys. `coeffs(as_symengine=True)` and `con_symengine()` return
symengine.py objects. Tests run in shared mode and check `_rcp_address` identity for coefficients.

**T5.6 Per-coefficient helpers.** `SDA.simplify()`, `SDA.expand()` and `SDA.subs(dict)`, each
returning a new SDA (loop over the slot's coefficients in C++).
- Acceptance: port the cases in `test/test_symbolic.cc` to pytest (same inputs, same expected
  values); `evaluate(exp(a*svar(0) + 1), {a: 0.3})` equals `exp(0.3*var(0) + 1)` to 1e-13;
  `evaluate(promote(x), {}) == x`. Benchmark: Python `exp(SDA)` at (3,3) takes at most 1.05× the
  C++ time (add this case to T0.5 and T2.5).

### Stage 6 — CSDA (`std::complex<SDA>`)

**T6.1 Class `CSDA`.** Ctor `CSDA(re: SDA, im: SDA | None = None)`; `promote(cnda) -> CSDA`
(`interop.h:317`); `real` and `imag` properties.

**T6.2 Arithmetic and functions.** Scalar set {CSDA, SDA, complex, float, Expr}, as supported by the
`is_da_coeff` templates in `da.h:98-230` (bind only what compiles; list the skipped combinations);
math functions; `abs` (`functions.h:73`). Accept `symengine.Basic` wherever `Expr` is accepted (as in T5.5).

**T6.3 Algorithms.** `CSDAList`; `cd_composition` for `T = Expression`;
`evaluate(csda, …) -> CNDA` (`interop.h:337,347`).
- Acceptance: port `test/test_symbolic_cd.cc` and the checks in
  `examples/example_complex_symbolic.cc` (for example, the constant parts of `exp(z)` equal
  `exp(1.5)*cos(0.5)` and `exp(1.5)*sin(0.5)`).

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

**T7.4 Tests.**
- Two envs with different orders; arithmetic in each truncates at the right order.
- Ops on env-B vectors while env A is current give correct results (every function in T1.5,
  T2.2 and T5.3).
- Mixing envs in one binary op raises `EnvError`.
- `close()` with live vectors: later use raises `EnvError`, and the process exits cleanly.
- `with env.order(n)` affects only that env.
- `inv_map` in a non-default env matches the default-env result.
- Benchmark: rerun T2.5; the gate must still pass.

### Stage 8 — Finish

**T8.1 Stubs.** Generate `python/miradac/_core.pyi` with `python -m nanobind.stubgen`; add
`py.typed`; run `mypy --strict` on the examples.

**T8.2 Examples.** Port `examples/examples.cc`, `example_interop.cc`, `example_complex_da.cc`,
`example_1_symbolic.cc` and `example_complex_symbolic.cc` to `python/examples/*.py`; each one runs
as a test.

**T8.3 Docs.** Docstrings on every binding; a "Python" section in the README with build
instructions and the A.11 example.

**T8.4 Performance report.** `python/bench/REPORT.md`: the final table from T2.5 plus the symbolic
cases, with machine and compiler info.

### Deferred (not in this plan)

- GIL release (needs per-env locks and a thread-safe SymEngine build); free-threaded Python.
- Pickling (needs a portable format that records order and nvars).
- Fast repeated SDA evaluation through SymEngine `Lambdify` (a big win for parameter scans).
- Portable wheels (`auditwheel` bundling SymEngine and GMP) and stable-ABI builds.
