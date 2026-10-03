# MiraDAC — Julia Binding Development Plan

**Goal:** A Julia package `MiraDAC.jl` exposing the MiraDAC library — NDA, CNDA, SDA, CSDA and the
multi-environment API — with call overhead low enough that Julia code driving DA arithmetic runs
close to C++ speed, and an in-place API that avoids allocation in hot loops.

This document has two parts:

1. **Part A — Design** (decisions, with the reason for each).
2. **Part B — Staged implementation plan** (ordered tasks, each with files, steps and an
   acceptance check, written so a small model can execute one task without making design
   decisions).

> **Ground rules for any implementer**
> - All design decisions are in Part A. Do not invent new ones. If a task cannot be done as
>   written (an API does not exist, a signature differs), make the smallest faithful adaptation,
>   record it as an *As built* note under the task, and continue; if it cannot be done at all,
>   stop and report.
> - Every stage must build and pass its tests before the next stage starts.
> - Do not modify C++ library code (`include/`, `src/`) except in tasks marked **[C++]**. Every
>   [C++] task keeps the C++ suite (`build/test/run_tests`, run from `test/`, and
>   `run_tests_asan`) and the Python suite (`.venv/bin/pytest python/tests -q`) green.
> - Do not change the Python binding (`python/`) except where a task says so.
> - `DEVELOPMENT_PLAN_PYTHON.md` and `python/src/*` are the reference for API coverage, naming
>   and for solutions already found (env guard, exception mapping, benchmark protocol). Reuse
>   the ideas; the Julia code is separate.
> - Do not commit, push or tag unless the user asks. Put throwaway files under a scratch
>   directory, never in the repo.

---

## Part A — Design

### A.1 Facts that drive the design

| Fact | Source | Consequence |
|---|---|---|
| Julia 1.13.1 is installed (`~/.local/bin/julia`). | `julia --version` | Develop and test on 1.13; support Julia ≥ 1.10 (the LTS). |
| CxxWrap.jl 0.17.5 needs `libcxxwrap_julia_jll` 0.14 and Julia ≥ 1.10; a CxxWrap glue library is compiled against one libcxxwrap-julia version. | General registry `C/CxxWrap/Compat.toml` | CxxWrap couples our binary to a third-party C++ ABI (A.2). |
| SymEngine.jl 0.13.x uses `SymEngine_jll` ≤ 0.12 (SymEngine C++ 0.12); MiraDAC pins SymEngine 0.14.0 (`cmake/symengine_pin.txt`). | General registry `S/SymEngine/Compat.toml`, `jll/S/SymEngine_jll/Versions.toml` | Zero-copy interop with SymEngine.jl is impossible (two different SymEngine copies, as with the PyPI wheel in the Python plan A.8). Interop is by string. Loading both in one session risks C++ symbol clashes; T5.6 tests it and has a fallback. |
| Julia's GC is not reference counting; a finalizer runs at some later GC and "must not cause a task switch" (`@doc finalizer`); finalizers may run on another thread than the one using the object. | Julia docs | DA temporaries are not freed promptly, so the fixed-size pool can run out even when little is live. Frees must be deferred, and exhaustion must trigger a GC and a retry (A.5). |
| Pool exhaustion throws a plain `std::runtime_error("Pool::assign: Run out of vectors")`. | `include/da/pool.h` | A C API cannot tell exhaustion from other runtime errors reliably; add a dedicated exception type (T1.1 [C++]). |
| Every C++ constructor uses the current thread-local env; operations use their operands' env. `da_exchange_env` swaps the current env and returns the previous one. | `include/da/env.h`, Python plan A.4 | The C API uses the same `EnvGuard` idea as the Python binding (A.4). |
| Retired envs keep their shell; a DA vector of a retired env must not be used, but its destructor is safe. | `include/da/env.h` | Julia finalizers may run after `close(env)`: freeing is always safe, using is checked (A.4). |
| The Python binding reached ≈ 50 ns (in place) and 75–130 ns (new object) overhead per call. | `python/bench/REPORT.md` | Same performance gate for Julia (A.9). |

### A.2 Binding technology: a C API (`libmiradac_c`) plus `ccall`

| | CxxWrap.jl | C API + `ccall` (**chosen**) |
|---|---|---|
| Code to write | Less: C++ types and templates are wrapped directly. | More, but mechanical: one small `extern "C"` function per operation. |
| Binary coupling | The glue library must be rebuilt for each libcxxwrap-julia ABI; the JLL must track CxxWrap's compat bounds. | Plain C ABI; no third-party C++ ABI. |
| Memory control | CxxWrap boxes every returned object and frees it from its own finalizer, which calls the C++ destructor on whatever thread runs the finalizer, racing with the pool (A.5). | Full control: deferred frees, retry on exhaustion, in-place operations. |
| Call cost | CxxWrap dispatch plus a heap-allocated box per result. | `ccall` is a direct C call (a few ns). |
| Other languages | Julia only. | The same C API can serve Fortran, Rust, MATLAB/Octave, etc. |
| Packaging (BinaryBuilder) | Needs `libcxxwrap_julia_jll`, per Julia version. | One `MiraDAC_jll` with no Julia-specific dependency. |

The C API lives in this repository (`capi/`), is built by CMake when `DA_BUILD_CAPI=ON`, and is
installed next to the C++ library. The Julia package lives in `julia/MiraDAC/` in this repository
(the General registry supports packages in a subdirectory).

### A.3 C API conventions (`capi/include/miradac.h`)

- Prefix `mdac_`. Pure C header (C99), usable from C and C++; `extern "C"` guarded.
- Opaque handles, each a pointer to a heap-allocated C++ object:

  | Handle | C++ object |
  |---|---|
  | `mdac_env*` | `da::DAEnv*` (not owned: envs are owned by the C++ env registry) |
  | `mdac_nda*` | `new da::NDA` |
  | `mdac_cnda*` | `new std::complex<da::NDA>` |
  | `mdac_expr*` | `new SymEngine::Expression` |
  | `mdac_sda*` | `new da::SDA` |
  | `mdac_csda*` | `new std::complex<da::SDA>` |
  | `mdac_ndalist*` (and `cnda`/`sda`/`csda` lists) | `new std::vector<...>` |

- **Status codes.** Every function that can fail returns `mdac_status` (an `int`); results go
  through out-parameters. `MDAC_OK = 0`, `MDAC_ERR_ENV = 1` (env mismatch, cleared env, no env),
  `MDAC_ERR_VALUE = 2` (`std::invalid_argument`, `std::domain_error`), `MDAC_ERR_INDEX = 3`
  (`std::out_of_range`), `MDAC_ERR_POOL = 4` (`da::PoolExhausted`, T1.1), `MDAC_ERR_RUNTIME = 5`
  (anything else), `MDAC_ERR_UNSUPPORTED = 6` (a symbolic function in a numeric-only build).
- `const char* mdac_last_error(void)` returns the message of the last failure on the calling
  thread (thread-local `std::string`).
- **No exception crosses the C boundary.** Every function body is wrapped by one macro pair,
  `MDAC_TRY { ... } MDAC_CATCH`, which maps exceptions in this order: `da::PoolExhausted`,
  `EnvError` (a C API type, like the Python one), `std::invalid_argument`, `std::domain_error`,
  `std::out_of_range`, `std::logic_error` (from `check_env` → `MDAC_ERR_ENV`), `std::exception`,
  `...`.
- **Allocating form and in-place form** for every arithmetic operation and math function:
  `mdac_nda_add(const mdac_nda* a, const mdac_nda* b, mdac_nda** out)` allocates a new result;
  `mdac_nda_add_into(mdac_nda* out, const mdac_nda* a, const mdac_nda* b)` writes into an existing
  vector (copy-assign then compound op; no new slot). Scalar variants take `double` (`_d`
  suffix); complex variants take `double re, double im` (`_z` suffix).
- **Free functions never fail:** `void mdac_nda_free(mdac_nda*)` (one per handle type) is safe on
  a null pointer and on a vector whose env was cleared.
- **Size-query convention** for variable-length outputs: `f(..., T* buf, size_t cap, size_t* n)`
  writes at most `cap` items and always sets `*n` to the full size, so the caller can size the
  buffer and call again.
- **Indices are 0-based in C.** The Julia layer converts (A.6).

### A.4 Env handling in the C API

Every C API function that touches a DA object creates, on entry, an `EnvGuard` on the env of its
first DA argument, exactly as `python/src/common.h` does:

```cpp
struct EnvGuard {                       // capi/src/common.h
    da::DAEnv* prev;
    explicit EnvGuard(da::DAEnv* e) {
        if (e->retired()) throw EnvError("DA environment has been cleared");
        prev = da::da_exchange_env(e);
    }
    ~EnvGuard() { da::da_exchange_env(prev); }
};
```

Functions without a DA argument that create a vector take an explicit `mdac_env*` (never null
in the C API; the Julia layer passes the current env).

### A.5 Memory: deferred frees, retry on exhaustion, in-place API

1. **Julia objects.** `mutable struct NDA; ptr::Ptr{Cvoid}; end` (the same shape for every handle
   type). The inner constructor registers a finalizer.
2. **Deferred free.** The finalizer never calls C. It pushes the pointer onto a per-type free
   queue `FREE_QUEUE::Vector{Ptr{Cvoid}}` protected by a `Threads.SpinLock`, using the safe
   pattern `if trylock(lk); push!(q, p); unlock(lk); else finalizer(f, obj); end` (re-register
   and retry at the next GC). The queue is drained, calling `mdac_*_free` for each pointer, at the
   start of every allocating API call; an `isempty` check comes first, so the common case costs
   one load.
3. **Retry on exhaustion.** Every allocating call goes through one helper. If the status is
   `MDAC_ERR_POOL`, it drains the queue, runs `GC.gc(false)`, drains again and retries; if the
   pool is still full, it runs `GC.gc(true)`, drains and retries once more; then it throws
   `PoolExhaustedError`.
4. **In-place API** for hot loops, which allocates nothing: `add!(out, a, b)`, `sub!`, `mul!`
   (extends `LinearAlgebra.mul!`), `div!`, and `exp!(out, a)`-style forms for the math functions.
5. **Threads.** One env may be used by one Julia thread at a time (the C++ pools have no locks);
   different threads may use different envs. Documented, not enforced.

### A.5b Scopes (`dascope`, `keep!`) — branch `julia-alloc-scope`

**Why.** Outside a scope an allocating operation (`a + b`, `exp(a)`) pays (1) a Julia object, a
finalizer registration and the deferred free (~150–350 ns), and (2) a cold pool slot: Julia frees
a result only after a GC, so a new result never reuses the warm slot just freed, as C++ does since
`efc0f70` (most recently freed slot first); at (6,10) a 64 KB slot is zeroed cold, ~1 µs.

**Design (fixed by the user).**

1. `dascope(f)` / `dascope() do … end` runs `f` and returns its value.
2. Every `NDA`, `CNDA`, `NDAList`, `CNDAList` created while a scope is active on the current task
   is registered in the innermost scope (SDA/CSDA join when they exist). The scope stack is
   task-local (`task_local_storage`, works on Julia ≥ 1.10).
3. Such objects get no finalizer. At scope exit — return or exception — every registered object
   that is not kept is freed at once through `mdac_*_free` and invalidated (`ptr = C_NULL`).
4. Kept: objects reachable from the return value (itself, or inside `Tuple`, `NamedTuple`,
   `AbstractArray`, `AbstractDict` values, recursively) and objects passed to `keep!(x)`. A kept
   object moves to the enclosing scope; at the outermost scope it gets its deferred-free
   finalizer and is left to the GC.
5. Using a freed object throws `MiraDAC.FreedObjectError` ("return it from the dascope block or
   call keep!(x)"). The check is one pointer comparison in `unsafe_convert(Handle, x)`, which
   every `ccall` and the list constructors use; C never sees `C_NULL`.
6. Scopes nest; an exception frees everything not kept and is rethrown.
7. Outside any scope nothing changes (finalizer, deferred free, retry). In-place operations
   create no objects.
8. Exported: `dascope`, `keep!`, `FreedObjectError`; user docs in `julia/MiraDAC/README.md`
   ("Scopes").

**As built** (`julia/MiraDAC/src/scope.jl`, no C API change):

- Every inner constructor ends in `adopt!(new(p))`: if no scope is active it attaches the
  finalizer, otherwise it pushes the object onto the scope's vector for its type (one typed
  `Vector` per type, so freeing needs no dynamic dispatch).
- A global `ACTIVE_SCOPES` atomic counts open scopes over all tasks; while it is 0, `adopt!`
  skips the task-local lookup, so the unscoped path costs one load more than before.
- The per-task `ScopeStack` keeps its `Scope` objects (and their vectors' capacity) for reuse:
  entering a scope allocates nothing. The kept set (`Base.IdSet`, qualified: Base exports
  `IdSet` only from Julia 1.11) is created only when something is kept; containers are recorded
  in it too, which guards against cycles. Arrays of bits types are skipped, and `#undef` slots
  of a partly filled array (`Vector{NDA}(undef, n)`) are skipped with `isassigned`. A lazy
  array wrapper (`view`, `v'`, `reshape`, ...: any non-`Array` with `parent(x) !== x`) keeps its
  parent instead of being indexed: indexing a view of an `NDAList` would copy the elements (new
  pool objects) and leave the list itself unkept, and indexing `v'` calls `adjoint(::NDA)`,
  which does not exist.
  `test_scope.jl` passes on Julia 1.10.10 and 1.13.1.
- Objects are freed in creation order, so the newest (warmest) slot is freed last and handed
  out first by the pool.
- Freeing at exit does not go through the env-aware drain of A.5: the objects belong to envs
  this task is using, which (A.5 item 5) no other thread uses at the same time.
- A scope's objects stay allocated until it ends, so a scope around a loop of N operations needs
  N slots: put the scope inside the loop (or around batches) and carry results in place
  (`add!(x, x * 0.5, davar(1))`). `test_scope.jl` checks 100 000 operations in scopes with a pool
  of 64: no retry, `count` back to baseline.

**Benchmark** (`bench_ops.jl --scoped`): the allocating cases (`add`, `mul`, `mul_const`, `exp`,
`cmul`, `cexp`, `composition`) also run as `scoped_*`, one `dascope` per batch of 16 operations
(`--batch N`). Strict rule for them: overhead ≤ 150 ns when C++ < 1 µs, ratio ≤ 1.10 when
C++ ≥ 1 µs. Results: `julia/MiraDAC/bench/REPORT.md`, section "A.5b". There are 21 scoped
cases (7 operations × 3 sizes). All 21 pass the strict rule in both recorded runs (core 3, load
average 1.8–4.3) and in both of an independent verifier's runs. The script prints the A.5b
(scoped) and A.9 (unscoped) gates on separate lines and exits 1 if either fails.

The unscoped `n6o6` allocating cases sit at the A.9 limit (+294 to +439 ns against 400), and
whether they pass depends on load: one recorded run passes, one fails on `n6o6/mul_const` (+424 ns),
and both verifier runs fail on `n6o6/add`. An A/B against HEAD shows the same times. The cause is
not in the binding: raw C calls show that a result in a cold slot (unscoped Julia frees only after
a GC) costs 250–300 ns more than one in a warm slot at n6o6. Item 7 keeps that behaviour outside a
scope, and the fix is `dascope` itself. The A.9 allowance for unscoped allocating operations is a
user decision and is not changed here.

### A.6 Julia API (target)

- **1-based indices** everywhere in Julia: `davar(1)` is the first variable, coefficient index 1
  is the constant term. The Julia layer converts to 0-based C calls.
- DA types are **not** subtypes of `Number`: they are not field elements with Base's numeric
  promotion rules, and subtyping `Number` would pull in generic methods that allocate or give
  wrong results. Operators are defined explicitly.
- Names chosen not to clash with `Base`, `Statistics`, `LinearAlgebra`, `SpecialFunctions` or
  SymEngine.jl:

| Julia | Meaning |
|---|---|
| `DAEnv(order, nvars, poolsize; table=false)`, `close(env)`, `with_env(f, env)` / `with_env(env) do … end`, `current_env()`, `default_env()` | multi-env API |
| `init!(order, nvars, poolsize)`, `clear!()` | default env (C++ `da_init` / `da_clear`) |
| `with_order(f, n)` | temporary truncation, nests correctly (Python plan T1.6 algorithm) |
| `NDA(x::Real)`, `NDA(coeffs::Vector{Float64})`, `davar(i)` | numeric DA vectors |
| `CNDA(re::NDA, im::NDA)`, `CNDA(z::Complex)` | complex numeric |
| `SymExpr(x)`, `dasymbols("a b")` | symbolic scalar (not `Expr`, which is Julia's AST type) |
| `SDA(x)`, `sdavar(i)`, `promote_sda(v::NDA)` | symbolic DA vectors |
| `CSDA(re::SDA, im::SDA)` | complex symbolic |
| `+ - * / ^`, `exp log sqrt sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh abs` (Base methods) and `MiraDAC.erf` (own function, exported; no SpecialFunctions dependency) | arithmetic and functions |
| `add! sub! mul! div!`, `exp!` … | in-place forms (A.5) |
| `con(v)`, `coeff(v, exps)`, `coeffs(v)`, `nterms(v)`, `norm(v)` (`LinearAlgebra.norm`), `exponents(env)` | inspection |
| `der(v, i)`, `integ(v, i)`, `substitute(v, i, x)`, `compose(f, g)`, `inv_map(m)`, `evaluate(v, dict)` | algorithms |
| `NDAList`, `CNDAList`, `SDAList`, `CSDAList` | maps held in C++ (no per-element copy), convertible from `Vector` |
| `evaluate_map(m, pts::Matrix{Float64})` | points as **columns** (`nvars × N`, column-major, contiguous); result `length(m) × N` |
| `show(io, v)` compact; `show(io, MIME"text/plain"(), v)` full table (C++ `operator<<`) | display |

- Exceptions: `EnvError` and `PoolExhaustedError` (both `<: MiraDACError <: Exception`);
  `MDAC_ERR_VALUE` → `ArgumentError`, `MDAC_ERR_INDEX` → `BoundsError`, and `DomainError` for
  math functions outside their domain (for example `asin` of a constant part above 1).

### A.7 Loading the library

- `const libmiradac = @load_preference("libmiradac", "")` (Preferences.jl), so `ccall` sees a
  constant (fast) and changing the path triggers recompilation.
- `MiraDAC.set_library!(path)` writes the preference to `LocalPreferences.toml`; the dev script
  `julia/dev_setup.jl` sets it to the absolute path of `build/capi/libmiradac_c.so`.
- `__init__` checks the library: `mdac_abi_version()` must equal the package's expected ABI
  version (an integer constant bumped whenever the C API changes incompatibly), else a clear
  error naming the path and both versions.
- Stage 9 adds `MiraDAC_jll` as the default when no preference is set.

### A.8 Symbolic interop with SymEngine.jl

- Through strings only (A.1): a package extension `ext/MiraDACSymEngineExt.jl`, loaded when both
  MiraDAC and SymEngine are loaded, defines `SymEngine.Basic(x::SymExpr)` and `SymExpr(x::Basic)`
  via `string` and parsing.
- Same-process coexistence of the two SymEngine libraries is tested in T5.6; if it fails, the
  fallback there links SymEngine statically with hidden symbols.

### A.9 Performance gate

Same as the Python plan's revised T2.5 gate. A Julia script runs each case round by round,
interleaved with a fresh `bench_cpp --driven` child (`python/bench/bench_cpp.cpp`), both pinned
to one core, and takes the minimum per case:

Each case gets exactly one rule, chosen by its C++ time:

| Julia operation | C++ time < 1 µs | C++ time ≥ 1 µs |
|---|---|---|
| In-place (`add!`, `mul!`, `exp!`, …, and `+=`-style compound forms) — the fast path | overhead ≤ 60 ns | ratio ≤ 1.10 |
| Allocating (`a + b`, `exp(a)`, …) | overhead ≤ 400 ns **or** ratio ≤ 1.5 | overhead ≤ 400 ns **or** ratio ≤ 1.5 |

The benchmark uses a realistic pool: 10 000 slots at (nvars, order) = (3,4) and (6,6), 2 000 at
(6,10) (each slot there is 64 KB).

*Revised 2026-09-30, approved by the user.* The first version asked ≤ 150 ns for allocating
operations, like Python. Julia cannot reach that: it has no reference counting, so every new DA
object pays for a finalizer registration and the deferred free (~200 ns), plus its share of the
GC that the pool-exhaustion retry triggers (one young collection, ~100 µs, per pool-full of
operations: ~290 ns per operation with a pool of 400, ~10 ns with 10 000). The in-place API is
the performance path and has the tight limit; the documentation (T8.2) tells users to size the
pool generously and to use the in-place forms in hot loops. The in-place limit applies only
below 1 µs: 60 ns is below the noise of operations that take milliseconds.

*Revised again 2026-10-01, user decision.* Allocating operations of ≥ 1 µs could not meet
ratio ≤ 1.10 once the C++ pool hands out the most recently freed slot first (commit `efc0f70`):
C++ then reuses a warm slot for every result, while Julia, which frees only after a GC, gets a
cold one and zeroes it (about 1 µs per result at (6,10), where a slot is 64 KB). The rule for
allocating operations is now "overhead ≤ 400 ns or ratio ≤ 1.5" at every size; in-place rules are
unchanged. An explicit scope that frees temporaries promptly (`dascope`) is being developed
separately (branch `julia-alloc-scope`) to remove this cost.

*2026-10-01, later:* even with that rule, `a * 2.0` at (6,6) sits at the 400 ns limit (+319 to
+423 ns over 4 runs; C++ itself dropped from ~390 to ~300 ns with the slot stack) — the same
cold-slot cost below 1 µs. Following the user's direction to keep option 1 and continue while
`dascope` is developed, **allocating cases are informational until `dascope` lands**: the
benchmark still measures and records them, but they do not block a stage. In-place cases stay
blocking under their strict rules. Once `dascope` is merged, scoped allocating cases become
blocking under the strict rule (overhead ≤ 150 ns below 1 µs, ratio ≤ 1.10 above).

*2026-10-02: `dascope` merged* (branch `julia-alloc-scope`, see A.5b). `bench_ops.jl --scoped`
now has three classes: in-place cases (blocking, strict), scoped allocating cases `scoped_*`
(blocking, overhead ≤ 150 ns below 1 µs, ratio ≤ 1.10 above), and unscoped allocating cases
(informational, overhead ≤ 400 ns or ratio ≤ 1.5); it exits non-zero only when a blocking case
fails. Unscoped allocating operations stay informational: their extra cost is inherent to
freeing by GC, and `dascope` is the documented way to avoid it. After the merge, `SDA`, `CSDA`
and their lists are scoped like `NDA` and `CNDA`, and `show` of a freed object prints
`T(freed)` instead of throwing, so an error raised inside a scope keeps its message (both found
in the final dascope verification).

---

## Part B — Staged implementation plan

Conventions for every task:
- **SymEngine:** `eval "$(scripts/setup_symengine.sh --print-env)"` before any symbolic build.
- **C++ and C API build:**
  `cmake -S . -B build -G Ninja -DDA_BUILD_CAPI=ON -DSymEngine_DIR="$SymEngine_DIR" && cmake --build build`
  (library: `build/capi/libmiradac_c.so`).
- **C++ tests:** from `test/`: `../build/test/run_tests` and `../build/test/run_tests_asan`.
- **C API tests:** `ctest --test-dir build -R capi`.
- **Julia tests:** `julia --project=julia/MiraDAC -e 'using Pkg; Pkg.test()'`.
- **Python tests (must stay green):** `.venv/bin/pytest python/tests -q`.
- "Acceptance" lists the exact checks. A task is done only when all pass.

### Stage 0 — Scaffolding

**T0.1 C API build target.**
- Files: `CMakeLists.txt`, `capi/CMakeLists.txt` (new), `capi/include/miradac.h` (new),
  `capi/src/capi_core.cpp` (new).
- Do: add `option(DA_BUILD_CAPI "Build the C API library miradac_c" OFF)`; when ON,
  `add_subdirectory(capi)`. In `capi/CMakeLists.txt`: shared library `miradac_c` from
  `capi/src/*.cpp`, linking `daStatic` privately, C++17, public include dir `capi/include`,
  `-fvisibility=hidden` with only `mdac_*` symbols exported (a `MDAC_API` macro using
  `__attribute__((visibility("default")))`), `INSTALL_RPATH_USE_LINK_PATH ON`. Install the
  library and `miradac.h` (to `include/miradac.h`). First functions: `mdac_abi_version()`
  returning `MDAC_ABI_VERSION` (1) and `mdac_version()` returning the project version string (a
  compile definition from `PROJECT_VERSION`).
- Acceptance: with `DA_BUILD_CAPI=OFF` nothing changes; with ON, `libmiradac_c.so` builds;
  `nm -D --defined-only build/capi/libmiradac_c.so` lists only `mdac_*` symbols (plus
  toolchain-defined ones such as `_init`/`_fini`); C++ and Python suites green.
- *As built:* `daStatic` goes into a shared library, so `DA_BUILD_CAPI=ON` also sets
  `CMAKE_POSITION_INDEPENDENT_CODE ON` (as the `SKBUILD` build does). `daStatic` is compiled with
  default visibility, so `miradac_c` also links with `-Wl,--exclude-libs,ALL` to keep its symbols
  out of the export table. `add_subdirectory(capi)` comes after the tests block, so its
  `add_test` sees `enable_testing()`; the C test is built only when `DA_BUILD_TESTS` is ON.
  std:: template instantiations (e.g. `std::vector<unsigned>` members from T1.3) keep default
  visibility under `-fvisibility=hidden` and leak in unoptimized builds, so the library also links
  with the version script `capi/src/exports.map` (`global: mdac_*; local: *;`).

**T0.2 Pure-C header check.**
- Files: `capi/test/test_c_header.c` (new), `capi/CMakeLists.txt`.
- Do: a C (not C++) program that includes `miradac.h`, calls `mdac_abi_version()` and
  `mdac_version()`, and returns 0 if they match the expected values. Register it with `add_test`
  as `capi_c_header`. Enable the C language only in `capi/` (`enable_language(C)`).
- Acceptance: `ctest --test-dir build -R capi` passes; the compile command for the test uses the
  C compiler (check `build/compile_commands.json` or the verbose build log).

**T0.3 Julia package skeleton.**
- Files: `julia/MiraDAC/Project.toml`, `julia/MiraDAC/src/MiraDAC.jl`,
  `julia/MiraDAC/test/runtests.jl`, `julia/dev_setup.jl` (all new); `.gitignore` (add
  `julia/MiraDAC/LocalPreferences.toml` and `julia/**/Manifest.toml`).
- Do: `Project.toml` with a fresh UUID (`using UUIDs; uuid4()`), `version` equal to the C++
  project version, `[deps]` Preferences, LinearAlgebra, Libdl; `[compat]` julia = "1.10",
  Preferences = "1". `src/MiraDAC.jl`: the module with `const libmiradac =
  @load_preference("libmiradac", "")`, `set_library!(path)`, `c_version()`, and `__init__`
  performing the check of A.7 (clear error if the preference is empty or the ABI differs).
  `julia/dev_setup.jl` sets the preference to the absolute path of
  `build/capi/libmiradac_c.so`.
- Acceptance: `julia --project=julia/MiraDAC julia/dev_setup.jl`, then
  `julia --project=julia/MiraDAC -e 'using MiraDAC; println(MiraDAC.c_version())'` prints the
  project version; `Pkg.test()` runs a first test checking the same.
- *As built:* with an empty preference `__init__` throws, so `using MiraDAC` fails and
  `set_library!` cannot be reached; it serves to change a working setup.
  `julia/dev_setup.jl` therefore activates `julia/MiraDAC`, instantiates it and writes the
  preference with `Preferences.set_preferences!(uuid, ...)` directly. `Test` is a test-only
  dependency (`[extras]`/`[targets]`).

**T0.4 C++ baseline for the Julia bench.**
- Do: reuse `python/bench/bench_cpp.cpp` (`--driven` mode: one benchmark round per stdin line,
  answered with `done`, JSON at end of input); build it with `-DDA_BUILD_BENCH=ON`. Write
  `julia/MiraDAC/bench/README.md` with the command and the protocol, as read from
  `python/bench/bench_ops.py`.
- Acceptance: `build/python/bench/bench_cpp --driven` answers `done` per input line and prints
  JSON at end of input.

**T0.5 CI job.**
- Files: `.github/workflows/ci.yml`.
- Do: a job `julia` (ubuntu-latest): the pinned SymEngine as in the symbolic job (same cache
  step), the C++ build with `-DDA_BUILD_CAPI=ON`, `julia-actions/setup-julia@v2` with a matrix of
  `1.10` and `1`, `julia julia/dev_setup.jl`, `Pkg.test()`.
- Acceptance: the job is green on push (report the run link).
- *As built:* the job runs `scripts/setup_symengine.sh --no-python` (symengine.py is not
  needed) with the symbolic job's cache step, and runs `ctest --test-dir build -R capi` before
  the Julia tests. The real run is pending a push; the job's commands were run locally in order
  with Julia 1.13.1 only (1.10 is not installed locally).

### Stage 1 — C API: core and NDA

**T1.1 [C++] `da::PoolExhausted`.**
- Files: `include/da/pool.h`, `test/test_pool.cc`.
- Do: `struct PoolExhausted : std::runtime_error { using std::runtime_error::runtime_error; };`
  in namespace `da`; `Pool::assign` throws it with the same message.
- Acceptance: a new Catch2 test: exhausting a pool throws `da::PoolExhausted`, which is still
  caught as `std::runtime_error`; the Python `test_pool_exhaustion_raises` is unchanged and green.

**T1.2 Error handling and env functions.**
- Files: `capi/src/common.h` (new), `capi/src/capi_core.cpp`, `capi/include/miradac.h`,
  `capi/test/test_capi.cc` (new, Catch2, reusing `test/catch.hpp`; registered as `capi_core`).
- Do: `EnvError`, `EnvGuard` (A.4), the status enum, the thread-local last error,
  `MDAC_TRY`/`MDAC_CATCH` (A.3). Functions: `mdac_last_error`,
  `mdac_init(order, nvars, poolsize, table)`, `mdac_clear()`, `mdac_env_current(mdac_env** out)`
  (`MDAC_ERR_ENV` if none), `mdac_env_default`,
  `mdac_env_make(order, nvars, poolsize, table, mdac_env** out)` (does **not** change the
  current env, like the Python `Env()`), `mdac_env_select`,
  `mdac_env_exchange(mdac_env* e, mdac_env** prev)`, `mdac_env_close` (`da_destroy_env`;
  `da_clear` for the default env), and the queries `order`, `max_order`, `nvars`,
  `full_length`, `poolsize`, `count`, `remain`, `retired`, `change_order(env, n, int* ok)`,
  `restore_order(env)`, `get_eps`/`set_eps`.
- Acceptance: one Catch2 case per function; a failing call leaves a non-empty
  `mdac_last_error()`.
- *As built:* the queries are `mdac_env_<query>(mdac_env*, unsigned* out)`; `order` is the
  current truncation order and `max_order` the order at creation (the Python `current_order` /
  `max_order`). Functions that cannot fail return their value directly: `mdac_env_retired`
  (`int`), `mdac_get_eps` (`double`), `mdac_nda_env` (`mdac_env*`), `mdac_last_error`.
  `mdac_env_exchange` accepts a null or retired env (it only swaps, so restoring never fails).
  Env handles are raw pointers, so, as in `python/src/bind_env.cpp`, `mdac_init`,
  `mdac_clear` and `mdac_env_close` keep every env's shell (one never-freed slot makes C++
  retire instead of delete). `MDAC_CATCH` calls one out-of-line function that rethrows and maps
  the exception. The test executable links its own `daStatic` for the C++ reference results
  (the C API library's copy is hidden, so the two have separate envs); it is built twice,
  `capi_core` and `capi_core_asan` (against `miradac_c_asan`, an ASan build of the library).

**T1.3 NDA lifecycle and inspection.**
- Files: `capi/src/capi_nda.cpp` (new), the header, `capi/test/test_capi.cc`.
- Do: `mdac_nda_new(env, double x, out)`, `mdac_nda_from_coeffs(env, const double*, size_t n,
  out)`, `mdac_nda_var(env, unsigned i, out)` (via `da::da_base`), `mdac_nda_copy`,
  `mdac_nda_free`, `mdac_nda_env`, `mdac_nda_con`/`set_con`, `mdac_nda_length`,
  `mdac_nda_nterms`, `mdac_nda_norm`, `mdac_nda_coeffs` (size-query convention, memcpy),
  `mdac_nda_coeff(v, const int* exps, size_t k, double*)`, `mdac_nda_set_coeff`,
  `mdac_nda_index_term(v, i, int* exps, double*)`, `mdac_nda_iszero(v, eps, int*)`,
  `mdac_nda_clean`, `mdac_nda_reset`, `mdac_nda_to_string` (the C++ `operator<<` text,
  size-query convention).
- Acceptance: a Catch2 case for each, including: freeing a vector after `mdac_clear()` is safe
  (ASan clean); using it returns `MDAC_ERR_ENV`.
- *As built:* `mdac_nda_set_con` resets the vector to the constant (the Python `con` setter);
  lengths and counts are `size_t`; `mdac_nda_to_string` counts the terminating NUL in `*n` and
  always NUL-terminates when `cap > 0`; `mdac_nda_index_term` returns `MDAC_ERR_INDEX` for
  `i >= full_length`; negative exponents give `MDAC_ERR_VALUE`.

**T1.4 NDA arithmetic, allocating and in-place.**
- Files: `capi/src/capi_nda.cpp`, the header, the test.
- Do: for `add sub mul div`: `mdac_nda_<op>(a, b, out)`, `mdac_nda_<op>_d(a, double, out)`,
  `mdac_nda_d<op>(double, a, out)` (scalar on the left), and in-place
  `mdac_nda_<op>_into(out, a, b)`, `mdac_nda_<op>_d_into(out, a, double)`,
  `mdac_nda_d<op>_into(out, double, a)`; plus `neg`, `pow_i`, `pow_d` in both forms. In-place
  forms: `*out = *a; *out op= *b;` when `out` is neither `a` nor `b`; handle aliasing (`out == a`
  → `*out op= *b`; `out == b` → compute into a temporary, then assign). One macro per shape, so
  each operation is one line.
- Acceptance: Catch2 cases comparing every C API result with the direct C++ result (coefficient
  equality); in-place forms keep the output's slot (`count()` unchanged); aliasing cases give the
  right result.
- *As built:* `NDA::operator*=` and `/=` assign a new product (`*this = *this * other`), which
  moves another slot into the output, so `mul`/`div` `_into` call the engine's three-slot
  kernels (`ad_mult`, `ad_div`, as `operator*`/`operator/` do) when `out` is neither operand;
  `add`/`sub` follow the copy-then-compound form, also for `out == a`. Every other aliasing case,
  `ddiv_into`, `neg`/`pow` and the math functions (T1.5) compute a temporary and copy-assign it
  (the slot is kept; `count()` is unchanged after the call). `div_d_into` repeats the zero check
  of `operator/(NDA, double)`. The header declares the four operator families and the math
  functions with the macros `MDAC_NDA_BINOP_DECL` / `MDAC_NDA_FUNC_DECL`.

**T1.5 NDA math functions.**
- Files: `capi/src/capi_nda.cpp`, the header, the test.
- Do: `sqrt exp log sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh erf` in both
  forms (`mdac_nda_exp(a, out)`, `mdac_nda_exp_into(out, a)`), `mdac_nda_abs(a, double*)`.
- Acceptance: results equal the C++ functions; `asin` of a vector with constant part 2.0 returns
  `MDAC_ERR_VALUE`.

### Stage 2 — Julia: NDA core

**T2.1 Low-level layer.**
- Files: `julia/MiraDAC/src/capi.jl` (new), `src/MiraDAC.jl`.
- Do: one Julia function per C function used so far,
  `ccall((:mdac_..., libmiradac), ...)`, grouped as in the header; `check(st)` turning a status
  into the exception of A.6 with the text of `mdac_last_error()`; the free queues and `drain!`
  (A.5 item 2); `alloc_call(f)` implementing the retry of A.5 item 3, counting retries in
  `POOL_RETRIES::Threads.Atomic{Int}`.
- Acceptance: unit tests: `alloc_call` with a closure returning `MDAC_ERR_POOL` once, then
  `MDAC_OK`, retries once and increments `POOL_RETRIES`; `check` maps each status to the right
  exception type.
- *As built:* `check(st, x)` maps `MDAC_ERR_VALUE` to `DomainError(x, msg)` (math functions),
  `check(st)` to `ArgumentError`; other runtime statuses throw `ErrorException`. A pool has no
  lock, so a drain frees only the pointers of the env the calling operation uses (and of
  retired envs); the others are parked per env until a drain for their env. The finalizers'
  `SpinLock` is held only to swap the queue vector for one allocated beforehand (no GC can run
  under it; otherwise the finalizers of a GC on another thread all fail `trylock` and wait for a
  later GC); the parked table has its own `SpinLock` (neither lock yields, so a task never moves
  to another thread, and away from its thread-local current env, inside an operation).
  `alloc_call(f, x, d)`: `x` gives the env, `d` the `DomainError` value. With more than one
  thread, `GC.gc` returns without collecting when another thread is collecting (possibly only
  the young generation), and the finalizers then run on that thread, so the retry does
  `retry_rounds()` = 12 rounds (incremental, then full collections with a 1 ms
  `Libc.systemsleep`, which does not yield) instead of 2; with one thread it is exactly A.5
  item 3. `POOL_RETRIES` counts rounds. `MiraDAC.drain!()` frees every queue where no other
  thread uses an env.

**T2.2 Types and construction.**
- Files: `julia/MiraDAC/src/nda.jl` (new).
- Do: `mutable struct NDA; ptr::Ptr{Cvoid}; end` (not `<: Number`, A.6) with an inner
  constructor that registers the deferred-free finalizer. The constructors of A.6, `davar(i)`
  (1-based; `BoundsError` outside `1:nvars`), `Base.copy`, `Base.deepcopy_internal`, `con`,
  `coeffs`, `coeff(v, exps)`, `nterms`, `LinearAlgebra.norm`, `Base.iszero(v)` and
  `iszero(v; eps)`, `show` in both forms, `env(v)`.
- Acceptance: tests: `coeffs(NDA(coeffs(v))) == coeffs(v)`; the full `show` of `davar(1)` equals
  the C++ `da::base[0]` text from `mdac_nda_to_string`; `davar(0)` throws `BoundsError`.
- *As built:* `NDA(coeffs)` accepts any `AbstractVector{<:Real}`; `coeffs(v)` returns the
  stored coefficients up to the last non-zero one; `env(v)` is not exported (`MiraDAC.env`).
  The compact `show` is `NDA(order=…, nvars=…, nonzero=…)`, and `NDA(cleared env)` for a vector
  of a cleared env.

**T2.3 Operators and functions.**
- Files: `julia/MiraDAC/src/nda.jl`.
- Do: `Base.:+ - * /` for (NDA, NDA), (NDA, Real), (Real, NDA); unary `-`; `^` for `Integer`
  and `Real`; `Base.exp` … `Base.atanh`, `MiraDAC.erf`, `Base.abs` returning `Float64`. In-place
  `add! sub! mul! div!` (`mul!` extends `LinearAlgebra.mul!`) and `exp!` … for every function.
- Acceptance: tests of every operator against hand-computed coefficients; `@inferred a + b`;
  `@allocated add!(c, a, b) == 0` after a warm-up call; the reference files `test/exp_da.txt`,
  `log_da.txt`, `sqrt_da.txt`, `pow0p3_da.txt` matched through a Julia port of
  `compare_da_with_file` (read the file; compare with the tolerance the C++ test uses).
- *As built:* the port reads the tables with `mdac_nda_index_term` and uses the C++ `eps`
  (`1e-14`); `pow3_da.txt` is checked too. `a / 0.0` throws `ArgumentError`.

**T2.4 Env functions and order.**
- Files: `julia/MiraDAC/src/env.jl` (new).
- Do: `init!`, `clear!`, `current_env`, `default_env`, `with_order(f, n)` (save the current
  order; `change_order`; `try f() finally` restore: `restore_order` when the saved order is the
  max order, else `change_order(saved)`), eps getter and setter.
- Acceptance: nested `with_order(3) do; with_order(2) do … end; end` restores 3 and then the
  original order; the order is restored also when `f` throws.
- *As built:* the `DAEnv` handle type (properties `order`, `max_order`, `nvars`,
  `full_length`, `poolsize`, `count`, `remain`, `retired`) is defined here; its constructor and
  `close` are Stage 7. `init!` returns the default env; eps is `get_eps()` / `set_eps!(x)`;
  `with_order` throws `ArgumentError` for an order above `max_order`.

**T2.5 Memory behavior.**
- Files: `julia/MiraDAC/test/test_memory.jl` (new).
- Acceptance tests:
  1. With a pool of 64, a loop of 100 000 iterations of `x = x * 0.5 + davar(1)` completes (the
     retry path runs the GC and recovers), and `MiraDAC.POOL_RETRIES[] > 0`.
  2. Keeping more live vectors than the pool holds throws `PoolExhaustedError`, not a crash.
  3. After `clear!()`, using an old vector throws `EnvError`; a `GC.gc()` afterwards does not
     crash.
  4. Finalizers on other threads: `Threads.@spawn` tasks that only create and drop DA vectors of
     **their own** env; `GC.gc()` from the main thread; no crash, and after draining, `count()`
     of each env is back to its base count. Run with `julia -t 4`.
- *As built:* test 4 makes its envs with `mdac_env_make` and selects one per task with
  `mdac_env_exchange` (`DAEnv(...)`/`with_env` are Stage 7). It first failed with 4 threads
  (`PoolExhaustedError`) for the reasons in the T2.1 note, which the queue and retry changes
  there fix.

### Stage 3 — NDA algorithms, bulk API, performance gate

**T3.1 C API lists and algorithms.** `mdac_ndalist_*` (new, free, length, get (a copy), set,
push, from an array of `mdac_nda*`), and `der`, `integ`, `substitute` (number, NDA, multiple
bases), `compose` (NDA map, float point, complex point), `inv_map`,
`evaluate_map(list, const double* pts, size_t npts, double* out)` (points contiguous, one after
another), `exponents(env, …)` (size-query convention). Catch2 tests against the C++ functions.
- *As built:* `capi/src/capi_list.cpp`. A list holds vectors of one env (`push`, `set`, `from`
  return `MDAC_ERR_ENV` otherwise), so a Julia drain frees it with that env;
  `mdac_ndalist_env` is null for an empty list, and `mdac_ndalist_length` returns its value.
  Map results are new lists (`mdac_ndalist_compose`, `_substitute`, `_inv_map` → `mdac_ndalist**`),
  single-vector results new vectors (`mdac_nda_der`, `_integ`, `_substitute_d`, `_substitute`,
  `_substitute_multi`); the point forms write `length(m)` values (`compose_z`: `(re, im)` pairs).
  `mdac_ndalist_evaluate_map` uses a copy of the Python `eval_points` kernel (bit-equal to
  `compose_d`); `mdac_nda_eval(v, pt, n, out)` was added for the T3.2 callable. The checks
  C++ only asserts are made as in the Python binding: base index → `MDAC_ERR_INDEX`, duplicate
  ids, wrong lengths, non-zero constant parts of an `inv_map` → `MDAC_ERR_VALUE`. The ABI
  version stays 1 (functions only added).

**T3.2 Julia wrappers.** `NDAList` type (the constructor from `Vector{NDA}` copies each element),
`der`, `integ`, `substitute`, `compose`, `inv_map`, `evaluate_map(m, pts)` (A.6: an `nvars × N`
matrix of points as columns → a `length(m) × N` matrix), a callable `(v::NDA)(pt::AbstractVector)`,
`exponents(env)::Matrix{Int32}` (`nvars × full_length`, one column per monomial). Tests ported
from `python/tests/test_algorithms.py`, with the reference files in `test/`.
- *As built:* `julia/MiraDAC/src/algorithms.jl`, tests in `test/test_algorithms.jl`.
  `NDAList <: AbstractVector{NDA}` (`l[i]` a copy, `l[i] = v` copies into the element,
  `push!`), freed through its own deferred-free queue (`LIST_QUEUE`, drained by the list
  functions and by every retry); every map argument also takes a `Vector{NDA}`. Signatures:
  `der(v, i)`, `integ(v, i)`, `substitute(v, i, x::Real|NDA)`, `substitute(v, ids, xs)` (`v` an
  `NDA` or a map; a map gives an `NDAList`), `compose(m, v)` (an `NDAList` for DA arguments,
  `Vector{Float64}`/`Vector{ComplexF64}` for a point), `inv_map(m, dim=length(m))`. A
  variable index outside `1:nvars` throws `BoundsError`; `evaluate_map` with a wrong number of
  rows `DimensionMismatch`.

**T3.3 Singular map.** `inv_map` of a singular map throws `ArgumentError` (C++ throws
`std::invalid_argument` since commit `3394e13`). Test both a zero row and a rank-deficient map.
- *As built:* in `test/test_algorithms.jl` ("singular map").

**T3.4 Bench script.** `julia/MiraDAC/bench/bench_ops.jl` with its own `bench/Project.toml`
(BenchmarkTools, and MiraDAC via `Pkg.develop(path=...)`). The cases of
`python/bench/bench_ops.py` (numeric) plus in-place `add!`, `mul!`, `exp!`; the C++ side from a
`bench_cpp --driven` child, interleaved; prints the table (case, C++ ns, Julia ns, overhead,
ratio, rule, pass/FAIL) and exits 1 on any failure.
- *As built:* the timing protocol is `bench_ops.py`'s (thread CPU time, minimum sample per
  case over interleaved rounds), which BenchmarkTools does not implement (wall time, its own
  sampling), so `bench/Project.toml` has MiraDAC (`[sources]` path), Printf and Random, not
  BenchmarkTools. `julia/dev_setup.jl` also writes the library preference into
  `bench/LocalPreferences.toml` (git-ignored), since a dependency reads its preferences from
  the active project. `iadd` is `add!(c, c, b)` (C++ `c += b`); `composition` is
  `compose(m, n)`, which returns a new list; `add!`, `mul!`, `exp!` are compared with the C++
  `add`, `mul`, `exp` cases. Default C++ binary: `build/python/bench/bench_cpp`.

**T3.5 Gate.** Run `taskset -c 5 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl`
on an otherwise idle machine (wait while the 1-minute load average is above 2). The gate of A.9
must pass. If it does not:
1. check `@code_warntype` on the failing operation and fix type instability first;
2. check the cost of the free-queue check on the fast path;
3. if it still fails, report with the full table.
Record the final table in `julia/MiraDAC/bench/REPORT.md`.
- *As built (first run, before the revised gate; superseded below):* the gate failed.
  Remedies 1 and 2 found nothing (`@code_warntype` clean; the empty free-queue check ~7 ns).
  (a) Allocating operations cost ~500 ns over C++ at every size: the pool of 400 fills every
  ~400 operations and the A.5 retry's `GC.gc(false)` (~100 µs) adds ~290 ns per operation;
  meeting 150 ns needs a design change to A.5 (for example a much larger default pool, or
  freeing without waiting for the GC), not a fix inside this plan. (b) `exp` and `compose`
  loops throw `PoolExhaustedError`: `src/engine.cpp` leaks its `pool.alloc()` temporaries when
  the pool runs out mid-function (one failed `da::exp` leaks 1 slot, `da_composition` 5), so the
  retry cannot recover; this needs a new **[C++]** task (exception-safe engine temporaries, with
  a Catch2 test that a failed `exp` leaves `count()` unchanged), recorded as `@test_broken` in
  `test/test_memory.jl`. (c) The in-place rule "overhead <= 60 ns at every size" is below the
  measurement noise at the large sizes. The bench warm-up needs two consecutive samples of
  >= 10 ms (one long collection otherwise fixed a tiny repetition count), and a case that
  throws is reported as a failed row.
- *Resolved 2026-09-30:* (b) is fixed in C++ by commit `a6be310` (engine temporaries are freed
  when a kernel throws; test `[numeric_leak]`); (a) and (c) are answered by the revised gate in
  A.9 (user decision). Remaining: turn the `@test_broken` in `test/test_memory.jl` into `@test`,
  add an `exp`/`compose` loop to the memory tests, update `bench_ops.jl` to the revised gate and
  pool sizes, and rerun the gate.
- *As built (gate passes, 2026-09-30):* the test and bench items above are done (`@test` for
  the retry in `exp`/`compose` loops with a pool of 64; `bench_ops.jl` has one rule per case by
  C++ time and the A.9 pools). The load rule is replaced by: pinned with `taskset -c 5`, started
  when core 5 is at most 25% busy over a 3 s sample of `/proc/stat` (the load average of this
  desktop rarely drops below 2). With the C++ fixes below in `build/` and `bench_cpp` (commit
  `57615a0`), the C++, ASan, C API, Python and Julia suites pass and the gate passes on all 27
  cases in two runs of two: in-place operations below 1 µs at most +17 ns over C++, allocating
  ones +145–219 ns, every case of 1 µs or more within ratio 1.10 (closest `n3o4/composition`
  1.086 and 1.078; `n6o10/mul_const` 0.868 and 0.910). Before, the gate failed for two
  systematic costs in the C++ library, each found from Julia and fixed by a **[C++]** commit:
  (1) `Layout::init_prod_index` allocated every row of `prdidx` separately, so in a fragmented
  heap (any Julia process; a C program that fragments its heap first showed the same) the rows
  scattered and the multiplication kernel paid in cache and TLB misses (`mul`, `exp`,
  `composition` 11–22% over C++; a C loop over `mdac_ndalist_compose` equalled C++, so not the
  binding); fixed by commit `3dec5a2` (`prdidx` one block owned by `Layout`). (2) `Pool::free`
  zeroed the whole slot, and a Julia result is freed only after a garbage collection, when its
  64 KB slot at `n6o10` has left the cache, so each allocating operation there paid ~2.5 µs of
  cold writes (`n6o10/mul_const` failed in 4 of 5 runs, ratios up to 1.379); fixed by commit
  `57615a0` (doubles zeroed in `assign`, not in `free`). Both final tables and the history are
  in `julia/MiraDAC/bench/REPORT.md`.

### Stage 4 — CNDA

**T4.1 C API.** `mdac_cnda_*`: new (from two NDAs; from `re, im` doubles in an env), free, copy,
`real`/`imag` (copies out), `set_real`/`set_imag`, `to_string`; arithmetic with CNDA, NDA
(`std::complex` templates), double, complex (`_z`), both sides, allocating and in-place;
functions `sqrt exp log asin acos atan asinh acosh atanh pow_i pow_d abs`; CNDA lists;
`cd_composition` (3 forms) and `compose` with complex points. Catch2 tests.
- *As built:* `capi/src/capi_cnda.cpp`, tests `[cnda]` in `capi/test/test_capi.cc`; the list
  handle helpers moved to `capi/src/common.h`. `mdac_cnda_new(re, im, out)` takes a null `im`
  for a zero imaginary part; `mdac_cnda_new_z(env, re, im, out)`. Operators: the suffix names the
  right operand and the prefix the left one, `mdac_cnda_<op>` (CNDA, CNDA), `_n`/`n<op>` (with an
  NDA), `_d`/`d<op>`, `_z`/`z<op>`, and `mdac_nda_<op>_z`/`mdac_nda_z<op>` (NDA with a complex,
  giving a CNDA), each also as `_into` (out first; `MDAC_CNDA_BINOP_DECL`). CNDA⊕NDA uses the
  `std::complex` templates (`da.h` has none); NDA − complex keeps `+imag(z)` as `da.h` defines it.
  Every `_into` form computes the C++ result in a temporary and copies it part by part into
  `out` (its slots are kept; `count()` is unchanged after the call), so it can also return
  `MDAC_ERR_POOL`. `mdac_cnda_neg` was added. The three `cd_composition` forms return new lists:
  `mdac_ndalist_compose_c` (NDA map, CNDA arguments), `mdac_cndalist_compose`,
  `mdac_cndalist_compose_n` (CNDA map, NDA arguments); C++ only asserts their checks, so
  arguments other than `nvars` vectors give `MDAC_ERR_VALUE` and vectors of other envs
  `MDAC_ERR_ENV`. "`compose` with complex points" is T3.1's `mdac_ndalist_compose_z` (C++ has no
  CNDA-map point form); a test checks it against `cd_composition` at constant CNDA arguments.
  The first line of `mdac_cnda_to_string` names the parts' slots, so the test compares the rest
  with the C++ text. The ABI version stays 1 (functions only added).

**T4.2 Julia.** `CNDA` type, `Base.real`/`Base.imag`, `Base.conj` only if C++ has it (check; skip
otherwise and note it), operators with `CNDA`, `NDA`, `Real`, `Complex` on both sides, functions,
in-place forms, `CNDAList`. Tests ported from `python/tests/test_cnda.py` and the reference files
`cd_calculation_*.txt`, `cd_composition_*.txt`, `da_composition_cd_*.txt`. Gate rerun for
`CNDA*CNDA` and `exp(CNDA)`.
- *As built:* `julia/MiraDAC/src/cnda.jl`, tests `test/test_cnda.jl` (with ports of
  `compare_cd_with_file`). `conj` is not defined: C++ has none. `CNDA(z::Number)` also takes a
  real. A.6 names no part setters, so the Julia API has none (the C functions are wrapped in
  `capi.jl`). `cd_composition` is `compose(m, v)` with `CNDA`s in `m` or `v` (`CNDAList` or
  `Vector`), returning a `CNDAList`. Functions are C++'s complex set (`sin(::CNDA)` is a
  `MethodError`); `abs(::CNDA)` is a `Float64`; `a / 0.0` throws `ArgumentError`. The in-place
  forms take temporary slots (T4.1), so they go through the `alloc_call` retry (`into_call`),
  and still allocate nothing in Julia (`@allocated == 0`). `CNDA` and `CNDAList` have their own
  free queues (`CNDA_QUEUE`, `CLIST_QUEUE`). The bench gained `cmul` (`c = a * b`) and `cexp`
  (`c = exp(a)`) on `CNDA(a, b)`, `CNDA(b, a)`, compared with `bench_cpp`'s cases of the same
  names under the allocating rule.
- *As built (gate rerun, 2026-10-01): FAILS on one case, `n3o4/cexp`.* Two runs of the protocol
  of T3.5 (`taskset -c 5`, core 5 at most 25% busy): every other case passes, including all
  `cmul` cases and `cexp` at `n6o6`/`n6o10` (ratios 1.000–1.045); `n3o4/cexp` (C++ ≈ 3.0 µs,
  ratio rule) is at 1.303 and 1.240 (earlier runs of the same code: 1.251, 1.199, and 1.265
  with a C++ pool of 10 000). Not the binding: a C loop over `mdac_cnda_exp_into` equals C++
  (3000 vs 2955 ns, pool 400). The cost is the C++ `Pool` free list: complex `exp` (`exp(re) *
  (cos(im) + i sin(im))`) takes ~20 temporary slots per call, and `Pool::alloc` hands out slots in
  free order; Julia frees results in GC order, which scatters the free list over the whole pool,
  so the temporaries land on unrelated pages. Reproduced from C with no Julia: pool 10 000 fresh
  3188 ns, after allocating 9 900 vectors and freeing them in random order 3743 ns (+17%; pool
  400: 3072 → 3140). The fix belongs in the C++ pool (for example `alloc` returning the lowest
  free slot, so temporaries stay packed), a new **[C++]** task, or the gate needs a user
  decision; both are outside this stage. Tables in `julia/MiraDAC/bench/REPORT.md`.
- *As built (gate rerun after `efc0f70` and the 2026-10-01 gate): passes in 2 runs of 4; the other
  two fail on `n6o6/mul_const` only.* `build/` and `bench_cpp` at commit `efc0f70`; the C++, ASan,
  C API (`[cnda]` 1426 assertions), Python (212) and Julia (504) suites pass. Same protocol (`taskset
  -c 5`, core 5 at most 25% busy over 3 s). `n3o4/cexp` now passes (ratio 1.056–1.095, overhead
  157–252 ns) as do all `cmul`/`cexp` cases. `n6o6/mul_const` (`a * 2.0`) sits on the 400 ns
  overhead limit: +319, +423, +404, +380 ns (ratio ≈ 2.3). The stack free list made C++ faster
  there (387–399 → 296–307 ns, it reuses the slot it has just freed), while Julia, whose results
  are freed only after a GC, got slower (605–617 → 627–727 ns): the cold-slot cost A.9's
  2026-10-01 revision describes, here on a 7.4 KB slot below 1 µs. Removing it needs the
  `dascope` work (A.9) or a gate decision, both outside this stage. Tables in
  `julia/MiraDAC/bench/REPORT.md`.
- *As built (final gate rerun, allocating cases informational per A.9 2026-10-01 later): the
  blocking gate PASSES; Stage 4 is done.* `bench_ops.jl` prints `gate (blocking): PASS/FAIL` over
  the in-place cases (strict rules; the only verdict that sets the exit status) and `allocating
  (informational): pass/fail` naming any allocating case over "overhead <= 400 ns or ratio <= 1.5".
  No case is non-allocating other than the in-place ones (`iadd`, `add!`, `mul!`, `exp!`), so
  those form the blocking set. One run (same protocol, core 5 17% busy; `build/` at `efc0f70`, no
  rebuild needed; C++, ASan, C API, `[cnda]` 1426 assertions, Python 212 and Julia 504 pass):
  exit 0, blocking PASS (in-place below 1 µs at most +12 ns, above at most ratio 1.018),
  allocating informational pass (`n6o6/mul_const` +380.6 ns, `n3o4/cexp` 1.083, `cmul`/`cexp`
  1.007–1.344). Table in `julia/MiraDAC/bench/REPORT.md`.

### Stage 5 — Symbolic: SymExpr and SDA

All symbolic C API functions are compiled only with `DA_WITH_SYMBOLIC`; otherwise they exist and
return `MDAC_ERR_UNSUPPORTED`. `mdac_has_symbolic()` reports the build; Julia's
`MiraDAC.HAS_SYMBOLIC` mirrors it, and the symbolic test files are skipped when it is false.

**T5.1 C API SymExpr.** `mdac_expr_*`: from double, int64, string (`SymEngine::parse`; bad input
→ `MDAC_ERR_VALUE`), symbol(name), free, copy, `to_string` (size-query), `to_double` (error if
free symbols remain), arithmetic with expr and double on both sides, `pow`, `neg`, `eq`, `hash`,
`subs(expr, n, const mdac_expr* const* keys, const mdac_expr* const* vals, out)`, `expand`,
`diff`, `free_symbols` (handles into a caller buffer, size-query), `simplify`.
- *As built:* `capi/src/capi_sda.cpp` (with T5.2), tests `[sym]` in `capi/test/test_capi.cc`.
  Every symbolic function body is wrapped in `MDAC_SYM(...)` (`common.h`), which in a numeric-only
  build becomes `return mdac::unsupported();` (`MDAC_ERR_UNSUPPORTED`, with a `mdac_last_error`
  message); the free functions then do nothing and the list length/env functions return 0/null
  (`MDAC_SYM_ELSE`). Names: `mdac_expr_new_d`, `_new_i` (`int64_t`, `<stdint.h>` now included),
  `_parse`, `_symbol`, `_copy`, `_free`, `_to_string`, `_to_double`, `mdac_expr_<op>` / `_<op>_d` /
  `_d<op>` for `add sub mul div pow`, `_neg`, `_eq(a, b, int*)`, `_hash(a, uint64_t*)`, `_subs`,
  `_expand`, `_diff` (`MDAC_ERR_VALUE` if the variable is not a symbol), `_simplify`
  (`da::simplified_expr`), `_free_symbols`; plus `mdac_expr_is_zero` (`da::is_zero`, for Julia
  `iszero`). `MDAC_CATCH` maps `SymEngine::SymEngineException` (parse errors, a symbol without a
  value in `evaluate`, ...) to `MDAC_ERR_VALUE`. The negative-exponent check moved from
  `capi_nda.cpp` to `common.h` (`mdac::exponents`). The ABI version stays 1 (functions only added).

**T5.2 C API SDA.** `mdac_sda_*` mirroring the NDA functions that compile for `T = Expression`
(try each; list the skipped ones in a comment, as the Python T5.1 did): lifecycle, `con`,
`coeff(exps) → expr`, `coeffs` (new expr handles, size-query), arithmetic with SDA, NDA
(`interop.h` mixed operators), expr and double on both sides, allocating and in-place; math
functions; `promote(nda)`, `svar`; `der integ substitute compose` for SDA;
`evaluate(sda, n, keys, vals, out_nda)` (the vector overload, `interop.h:175`); per-coefficient
`simplify`, `expand`, `subs`.
- *As built:* operand naming as for CNDA (T4.1): `mdac_sda_<op>` (SDA, SDA), `_n`/`n<op>` (NDA),
  `_e`/`e<op>` (Expr), `_d`/`d<op>` (double), `mdac_nda_<op>_e`/`mdac_nda_e<op>` (NDA with an Expr,
  giving an SDA), each also `_into` (`MDAC_SDA_BINOP_DECL`). Every `_into` form computes the C++
  result in a temporary and copy-assigns it into `out` (its slot is kept), so it can also return
  `MDAC_ERR_POOL`. Lifecycle: `mdac_sda_new(env, expr_or_null)`, `_new_d`, `_var` (`svar`),
  `_promote`, `_copy`, `_free`, `_env`; inspection `_con`, `_set_con`, `_length`, `_nterms`,
  `_coeffs`, `_coeff`, `_set_coeff`, `_index_term`, `_iszero`, `_clean`, `_reset`, `_to_string`
  (the first line names the slot, so the test compares the rest with C++). Not provided, listed in
  the file comment and the header: `from_coeffs`, `norm`, `abs`, eval at a point, the eps argument
  of `iszero`/`clean`, `asinh acosh atanh`, `inv_map`, `evaluate_map` (no SDA version in C++, or
  the kernel ignores symbolic coefficients). Lists `mdac_sdalist_*` as for NDA; algorithms
  `mdac_sda_der`, `_integ`, `_substitute_d`, `_substitute`, `_substitute_multi`,
  `mdac_sdalist_substitute`, `_compose`, `_compose_d` (length(m) new Expr handles). Checks C++ only
  asserts are made as for NDA. A numeric-only scratch build (`-DWITH_SYMBOLIC=OFF`) compiles and
  passes `capi_core`, with `mdac_has_symbolic() == 0` and the stubs returning `MDAC_ERR_UNSUPPORTED`.

**T5.3 Julia SymExpr.** `SymExpr` type: constructors from `Real` and `String`,
`dasymbols("a b")` → tuple, `Base.string`, `show`, `Float64(x)`, `==`, `hash`, arithmetic, `^`,
`subs(x, dict)`, `expand`, `diff(x, s)`, `free_symbols`, `simplify`. Tests ported from
`python/tests/test_expr.py`.
- *As built:* `julia/MiraDAC/src/symbolic.jl` (with T5.4), tests `test/test_expr.jl` (the sympy
  tests are not ported). `MiraDAC.HAS_SYMBOLIC` is a `Bool` global set in `__init__` from
  `mdac_has_symbolic()`; `runtests.jl` includes the symbolic files only when it is true (checked
  against a numeric-only scratch build: 504 tests pass, symbolic files skipped). SymEngine.jl 0.13
  exports `subs`, `expand`, `free_symbols` (and `coeff`, which MiraDAC exports since Stage 2), so,
  per A.6, MiraDAC does not export `subs`, `expand`, `free_symbols` (use `MiraDAC.subs`, ...);
  `simplify` and `dasymbols` are exported, `diff` extends `Base.diff`, `iszero` `Base.iszero`.
  `SymExpr(::Integer)` is an exact integer (beyond `Int64` through its decimal text), another
  `Real` a double; `Float64(x)` with free symbols throws `ArgumentError`; `show` prints
  `SymExpr("a + b")`, `print`/`string` the text. Arithmetic with an `Integer` operand is exact (as
  in Python). `SymExpr`s are freed through their own deferred-free queue (`EXPR_QUEUE`, no env:
  every drain frees them), since SymEngine's reference counts are not thread-safe
  (`WITH_SYMENGINE_THREAD_SAFE=OFF`); symbolic objects are therefore for one thread at a time.
  Julia 1.13 turns a literal `x^-1` into `inv(x)`, so `Base.literal_pow` is defined for `NDA`,
  `CNDA`, `SymExpr` and `SDA` to call `^`.

**T5.4 Julia SDA.** `SDA` type, `sdavar(i)`, `promote_sda(v)`, operators (with SDA, NDA,
SymExpr, Real), functions, in-place forms, `SDAList`, algorithms,
`evaluate(s, dict::AbstractDict{SymExpr,<:Real})::NDA`. Tests ported from
`python/tests/test_sda.py` and `python/tests/test_symbolic.py` (stored-reference tolerance 1e-12,
as the Python suite uses since commit `7ebbf09`).
- *As built:* tests `test/test_sda.jl`, `test/test_symbolic.jl`. `SDA()` (zero), `SDA(x::SymExpr)`,
  `SDA(::Integer)` (exact), `SDA(::Real)`; an `Integer` operand of an SDA operator or in-place form
  is an exact `SymExpr`. As for CNDA, the Julia API has no setters (`set_coeff`, `con =`), so the
  Python tests of them are not ported. `coeffs`/`coeff`/`con` return `SymExpr`s; `iszero(::SDA)`;
  `show` as for NDA. In-place forms `add! sub! mul! div!` for every operand pair and `exp!` ... for
  the SDA functions go through the `alloc_call` retry (they take temporary slots). `SDA` and
  `SDAList` have their own free queues (`SDA_QUEUE`, `SLIST_QUEUE`); a 500-iteration loop with a
  pool of 32 exercises the retry. Algorithms as for NDA: `der`, `integ`, `substitute` (number, SDA,
  several), `compose(m, v)` (an `SDAList`; at a `Real` point a `Vector{SymExpr}`; unlike the NDA
  kernel, the SDA one at a point takes nvars + length(m) temporary slots, so it too goes through
  the retry, tested with the pool filled by garbage).
  `simplify`, `MiraDAC.expand`, `MiraDAC.subs` act per coefficient and return a new SDA;
  `evaluate(s, d::AbstractDict{<:Any,<:Real})` converts each key with `SymExpr`, so T5.5's
  `SymEngine.Basic` keys work; `ArgumentError` if a symbol has no value.

**T5.5 SymEngine.jl extension.** `julia/MiraDAC/ext/MiraDACSymEngineExt.jl`, with
`[weakdeps] SymEngine` and `[extensions] MiraDACSymEngineExt = "SymEngine"` in `Project.toml`:
`SymExpr(b::SymEngine.Basic)` and `SymEngine.Basic(x::SymExpr)` via strings; the `SDA`
constructors and `evaluate` keys accept `SymEngine.Basic`. Tests in
`test/test_symengine_ext.jl` (SymEngine added to the test target's extras).
- *As built:* `[compat] SymEngine = "0.13"` (tested with 0.13.2, which loads `SymEngine_jll`'s
  libsymengine 0.11.2). SymEngine.jl prints with its Julia printer (`x^2`, `exp(1)`) and parses
  `**` too; MiraDAC's parser accepts `^`, so both directions round-trip ordinary expressions
  (tested: powers, rationals, `sin`, `exp`, `pi`, `E`, floats). The imaginary unit does not
  (SymEngine.jl prints `im`); complex symbolic work belongs in CSDA. The test uses
  `import SymEngine`: `using` would make the unqualified `coeff` ambiguous.

**T5.6 Coexistence of the two SymEngine libraries.**
- Test: in one Julia process, `using SymEngine` then `using MiraDAC`, run SDA arithmetic and
  SymEngine.jl arithmetic, and compare the results with known values; then the same with the
  opposite load order. Each order runs in its own subprocess (`run(`$(Base.julia_cmd()) …`)`).
- If both orders pass: record "no clash" in an *As built* note.
- If either fails (wrong result, crash, or undefined-symbol error), apply the **[build]
  fallback**: add `setup_symengine.sh --static` (builds the pinned SymEngine with
  `BUILD_SHARED_LIBS=OFF` and `-fPIC` into a separate prefix), link it into `libmiradac_c` with
  `-Wl,--exclude-libs,ALL` so none of its symbols are exported, and rerun the test. Report if it
  still fails.
- *As built: no clash; the fallback was not applied.* `test/test_coexistence.jl` runs, in a
  child Julia process per order (`import SymEngine` then `using MiraDAC`, and the reverse), SDA
  arithmetic (`exp(a + sdavar(1)) * (1 + davar(2))`, its coefficients and `evaluate`) and
  SymEngine.jl arithmetic (`expand((x + y)^3)`, `subs`, conversion both ways); both orders pass.
  Run standalone, both processes have both libraries loaded (SymEngine_jll's `libsymengine.so`
  0.11.2 and the pinned `libsymengine.so.0.14`).

### Stage 6 — CSDA

**T6.1 C API.** `mdac_csda_*`: lifecycle, `real`/`imag`, `to_string` (the C++ `operator<<` for
CSDA exists since commit `3394e13`), arithmetic with CSDA, SDA, expr, double and complex on both
sides, allocating and in-place; functions (the CNDA set); `promote(cnda)`; `cd_composition` for
`T = Expression`; `evaluate(csda, …) → cnda`.
- *As built:* `capi/src/capi_csda.cpp`, tests `[csda]` in `capi/test/test_capi.cc`; the Expr, SDA
  and CSDA handle helpers moved from `capi_sda.cpp` to `common.h`. Every function is wrapped in
  `MDAC_SYM` as in T5.1 (a numeric-only scratch build compiles, and `capi_core` checks that the
  CSDA stubs return `MDAC_ERR_UNSUPPORTED`). Lifecycle: `mdac_csda_new(re, im_or_null)`,
  `_new_z(env, re, im)`, `_promote(cnda)`, `_copy`, `_free`, `_env`, `_real`/`_imag` (copies),
  `_set_real`/`_set_imag`, `_to_string`. Operators named as for CNDA (T4.1): `mdac_csda_<op>`
  (CSDA, CSDA), `_s`/`s<op>` (SDA), `_e`/`e<op>` (Expr), `_d`/`d<op>`, `_z`/`z<op>`, and
  `mdac_sda_<op>_z`/`mdac_sda_z<op>` (SDA with a complex, giving a CSDA), each also `_into`
  (`MDAC_CSDA_BINOP_DECL`; temporary then part-by-part copy, as for CNDA, so it can return
  `MDAC_ERR_POOL`). No CSDA operator with an NDA or a CNDA (C++ has none; the Python binding has
  none either). `mdac_csda_neg`, `_pow_i`, `_pow_d`, the CNDA function set (`sqrt exp log asin
  acos atan asinh acosh atanh`, both forms) and `mdac_csda_abs` (an SDA, C++'s modulus).
  `mdac_csda_evaluate(v, n, keys, vals, out_cnda)` uses the vector overload, as for SDA. Lists
  `mdac_csdalist_*` as for CNDA; `cd_composition` forms, as new lists: `mdac_sdalist_compose_c`
  (SDA map, CSDA arguments), `mdac_csdalist_compose`, `mdac_csdalist_compose_s` (CSDA map, SDA
  arguments), with the checks of the CNDA forms. The ABI version stays 1 (functions only added).
  The tests keep symbolic constant parts out of the function and `pow` inputs: `pow(z, 3)` of a CSDA
  with symbolic constants took 12 s at order 4, and the complex functions ran out of memory
  (16 GB) on such input; with numeric constants and symbolic higher coefficients (the usage of
  `test_symbolic_cd.cc`) the whole `[csda]` set takes about 1 s.

**T6.2 Julia.** `CSDA` type and methods as for CNDA, plus `SymExpr` operands; `CSDAList`. Tests
ported from `python/tests/test_csda.py`, including the checks of
`examples/example_complex_symbolic.cc`.
- *As built:* `julia/MiraDAC/src/csda.jl`, tests `test/test_csda.jl` (plus a CSDA `evaluate`
  with `SymEngine.Basic` keys in `test/test_symengine_ext.jl`). `CSDA(re, im)`, `CSDA(re)`,
  `CSDA()` (exact zero parts), `CSDA(z::Number)` (doubles); `promote_sda(::CNDA)` gives a CSDA;
  `real`/`imag` copies, no part setters (A.6, as for CNDA). Operators and in-place forms (`add!
  sub! mul! div!`) for every pair of T6.1, with an `Integer` operand an exact `SymExpr` (as for
  SDA); `-`, `^`, `literal_pow`; `sqrt! exp!` ... for the CNDA function set; `abs(::CSDA)` is an
  `SDA`; `evaluate(z, d)` a `CNDA`; `compose(m, v)` for the three forms (a `CSDAList`). `sin`,
  `MiraDAC.erf`, ... of a CSDA and CSDA ⊕ NDA/CNDA are `MethodError`s. `CSDA` and `CSDAList`
  have their own free queues (`CSDA_QUEUE`, `CSLIST_QUEUE`); a 300-iteration loop with a pool of
  40 exercises the retry. Not ported: the list form of `evaluate` and CSDA operators with a
  `SymEngine.Basic` (Julia has neither for SDA). The CNDA check `@allocated == 0` for in-place
  forms is not made for CSDA (nor was it for SDA): `@allocated` counts SymEngine's GMP
  allocations, since SymEngine shares the `libgmp.so.10` Julia loads and Julia routes GMP's
  memory functions through its counted allocator (`add!` 960 bytes, `exp!` ≈ 20 KB per call;
  SDA `add!` 240 bytes, the same for every call, so no Julia object is left behind).

### Stage 7 — Multi-env API

**T7.1 C API.** `mdac_*_import(env, v, out)` (via `da::import_to`) and
`mdac_nda_promote_to(env, v, out)` (via `da::promote_to`); `mdac_*_env(v)` for every type.
- *As built:* `mdac_nda_import`, `mdac_cnda_import` (part by part), `mdac_sda_import`,
  `mdac_csda_import` and `mdac_nda_promote_to` (an `mdac_sda*`), next to each type's `_env`
  function (the `_env` functions already existed for every type and list). Each takes an
  `EnvGuard` on `v`'s env (`MDAC_ERR_ENV` if it was closed) and checks the target with `live`
  (moved from `capi_core.cpp` to `common.h`): null or closed → `MDAC_ERR_ENV`; different layouts
  → `MDAC_ERR_VALUE` (C++ `std::invalid_argument`). The symbolic three are `MDAC_SYM` stubs in a
  numeric-only build (checked in `capi: mdac_has_symbolic`). No `mdac_cnda_promote_to`: the
  plan names none, and Julia's `promote_sda(env, ::CNDA)` promotes the two parts. ABI version
  stays 1 (functions only added). Tests `[multienv]` in `capi/test/test_capi.cc`.

**T7.2 Julia `DAEnv`.** `DAEnv(order, nvars, poolsize; table=false)` (does not change the current
env), `==` and `hash` by pointer, properties (`order`, `max_order`, `nvars`, `full_length`,
`poolsize`, `count`, `remain`, `retired`), `close(env)` (idempotent), `with_env(f, env)`
(exchange, then `try f() finally` restore), `with_order(f, env, n)`, an `env` keyword on every
constructor and on `davar`/`sdavar`, `import_vec(env, v)`, `promote_sda(env, v)`.
- *As built:* in `src/env.jl`; `import_vec` and `promote_sda(env, v)` methods next to each type.
  `==` and `hash` are Julia's defaults for the immutable `DAEnv` struct, which compare and hash
  its one field, the pointer (tested). `with_env` throws `EnvError` for a closed env, exchanges
  with `mdac_env_exchange` and restores the previous env in `finally`. `with_order(f, n)` is now
  `with_order(f, current_env(), n)`. The `env` keyword (default `nothing`: the current env) is on
  `NDA(x)`, `NDA(coeffs)`, `davar`, `CNDA(z::Number)`, `SDA()`, `SDA(x)` (every form), `sdavar`,
  `CSDA()`, `CSDA(z::Number)`; constructors from DA parts (`CNDA(re, im)`, `CSDA(re, im)`) take
  the env of their parts, and lists have none. `promote_sda(env, ::CNDA)` is
  `CSDA(promote_sda(env, real(v)), promote_sda(env, imag(v)))` (see T7.1). Exports `with_env`,
  `import_vec`; `close` extends `Base.close`.

**T7.3 Tests** (port of `python/tests/test_multienv.py`): two envs with different orders;
operations on env-B vectors while env A is current; mixing envs throws `EnvError`; `close` with
live vectors, then use throws `EnvError` and a GC afterwards is safe; `with_order` affects one env
only; `inv_map` in a non-default env. Gate rerun.
- *As built:* `test/test_multienv.jl` (129 tests, in `runtests.jl` after `test_env.jl`): env
  creation and properties, `==`/`hash`/`show`, `with_env` nesting and restore on throw, the
  `env` keyword on every constructor, `import_vec`/`promote_sda(env, v)` for every type (layout
  mismatch → `ArgumentError`), mixing envs → `EnvError`, `close` with live NDA/CNDA vectors (use,
  `import_vec`, `with_env`, `with_order` on the closed env → `EnvError`; then `GC.gc()` and
  `MiraDAC.drain!()` are safe), closing the default env, `with_order(f, env, n)` on one env only,
  `inv_map` in a non-default env (also after `clear!()`), a default env without user vectors
  retired (not freed) by `clear!()`, a failed `DAEnv(4, 3, 2^31)` (`ErrorException`, the C++
  `std::bad_alloc` as `MDAC_ERR_RUNTIME` per A.3, where Python has `MemoryError`) keeping the
  current env, `erf` of an SDA in a non-default `table=true` env, and NDA, CNDA, SDA and CSDA operations
  (arithmetic, functions, in-place, `der`/`integ`/`substitute`/`compose`/`inv_map`/
  `evaluate_map`, `evaluate`) on env-B vectors giving the same result with env A current as with
  B current. Not ported: Python's `TypeError` checks (dispatch rejects those arguments with
  `MethodError` before any env logic) and the subprocess clean-exit test (the Julia suite itself
  exits after closing envs with live vectors). Suites: C++ and ASan 1600
  assertions, C API 3 ctest targets (`[multienv]` 657 assertions; numeric-only scratch build 55
  cases pass), Python 212, Julia 1156. Gate rerun (2026-10-01, two runs): both exit 0, blocking PASS
  (in-place below 1 µs at most +23.2 / +15.9 ns, above at most ratio 1.024 / 1.056); allocating
  informational pass in run 1, and in run 2 (loaded machine) only `n6o6/mul_const` over, +426.8 ns. Table in `julia/MiraDAC/bench/REPORT.md`.

### Stage 8 — Docs, examples, CI

**T8.1 Examples.** `julia/MiraDAC/examples/*.jl`: ports of `examples/examples.cc`,
`example_interop.cc`, `example_complex_da.cc`, `example_1_symbolic.cc`,
`example_complex_symbolic.cc`; each run by the test suite in a subprocess.
- *As built:* same names with `.jl` (from the Python ports, `python/examples/`). Each script
  calls `init!`, runs `main()` and `clear!()`; vectors print with `display` (the C++ table).
  `examples.jl` uses the allocating `substitute`/`compose` (A.6 has no output-argument forms) and
  skips C++'s `weighted_norm` (not in A.6); `example_complex_symbolic.jl` checks the evaluation
  with a local port of `compare_cd_vectors` (relative, per coefficient) and exits 1 on a
  mismatch, as does `example_interop.jl`. Output equals the C++ examples' up to float formatting
  and slot numbers (checked for `examples.jl`). `test/test_examples.jl` runs each in its own
  Julia process in a temporary directory (exit code 0, no "FAIL"), skipping the three symbolic
  ones when `HAS_SYMBOLIC` is false.

**T8.2 Docstrings and Documenter.** A docstring on every exported name;
`julia/MiraDAC/docs/` (Documenter.jl, its own `docs/Project.toml`): Home, Getting started, NDA,
CNDA, Symbolic, Multiple envs, Performance (the in-place API and the pool/GC behavior of A.5),
API reference (`@autodocs`). `makedocs` runs without warnings (`warnonly = false`).
- *As built:* only the module itself lacked a docstring; `runtests.jl` now checks
  `isempty(Docs.undocumented_names(MiraDAC))` (Julia >= 1.11 only, the function is newer than
  1.10). `docs/Project.toml` has Documenter 1 and MiraDAC by `[sources]` path; `julia/dev_setup.jl`
  instantiates it and writes its `LocalPreferences.toml` as for `bench` (and runs
  `Pkg.develop(path=...)` there on Julia 1.10, which ignores `[sources]`). The topical pages are
  prose with plain `julia` blocks (not run by Documenter; every block was run by hand), so each
  docstring appears once, in the API reference. `HTML(repolink = nothing, edit_link = nothing)`,
  since there is no public repository (otherwise Documenter warns). Output in `docs/build`
  (git-ignored). Build: `julia --project=julia/MiraDAC/docs julia/MiraDAC/docs/make.jl`, no
  warnings.

**T8.3 README.** A "Julia" section in the root `README.md`: build steps (C API and
`dev_setup.jl`), the API example, and the thread and pool notes.
- *As built:* section "Julia" after "Python", with the plan's commands, an example ported from
  the Python README example (run as a script), and notes on pool/GC, threads, envs, 1-based
  indices and SymEngine.jl.

**T8.4 CI.** Extend the T0.5 job: run the examples and `docs/make.jl` (no deploy).
- *As built:* steps "Examples" (each script in `$RUNNER_TEMP`, failing on a non-zero exit or
  "FAIL"; symbolic ones skipped in the numeric-only cell) and "Docs (Documenter, no deploy)". The
  YAML parses (`yaml.safe_load`); the job's commands were run locally in order with Julia 1.13.1
  for both matrix values of `symbolic` (the OFF cell in a scratch copy of `julia/` with
  `build/` linked to a numeric-only scratch build). The real CI run is pending a push; Julia 1.10
  was not run locally (not installed).

**T8.5 Numeric-only build.** Build with `-DWITH_SYMBOLIC=OFF -DDA_BUILD_CAPI=ON`; the Julia suite
passes with the symbolic files skipped. Add it to the CI job matrix.
- *As built:* matrix axis `symbolic: [ON, OFF]` on the `julia` job (so 4 cells); the SymEngine
  cache and setup steps run only for ON, and the configure passes
  `-DWITH_SYMBOLIC=${{ matrix.symbolic }}`. Locally: numeric-only scratch build, `ctest -R capi`
  3/3, Julia suite 612 tests pass (symbolic test files and symbolic examples skipped), examples and docs build.
  Symbolic build after Stage 8: C++ and ASan 1600 assertions, C API 3/3, Python 212, Julia 1163.

### Stage 9 — Distribution (needs the user's decisions first)

Stop and ask before starting this stage: the General registry and Yggdrasil require public
source, and the repository is private (as of 2026-09-30).

**T9.1 BinaryBuilder recipe.** `julia/binarybuilder/build_tarballs.jl` building `libmiradac_c`
with the pinned SymEngine linked statically with hidden symbols (the T5.6 fallback, used here in
any case, because a JLL must not clash with `SymEngine_jll`) and GMP from `GMP_jll`. Platforms:
`x86_64-linux-gnu` first, then those of `supported_platforms()` that build. Run locally with
BinaryBuilder.jl (Docker) and test the tarball with the Julia suite.
- *As built (first step, x86_64-linux-gnu, 2026-10-02):* `julia/binarybuilder/build_tarballs.jl`
  (`MiraDAC` 1.1.0) and `julia/binarybuilder/README.md` (the commands). Sources: the v1.1.0 commit
  `f32dc07` and the pinned SymEngine commit as a `GitSource`. An `ArchiveSource` of the pin's tarball
  is refused by BinaryBuilderBase: it is GitHub's automatic archive, whose checksum GitHub does not
  guarantee. The script checks the GitSource commit against the pin file of the checked-out MiraDAC.
  SymEngine is built with the pin's options but `BUILD_SHARED_LIBS=OFF`, and PIC, into
  `srcdir/symengine-static` (not `${prefix}`). MiraDAC is configured with
  `-DDA_BUILD_CAPI=ON -DDA_BUILD_TESTS=OFF -DDA_BUILD_EXAMPLES=OFF -DWITH_SYMBOLIC=ON`, and only
  target `miradac_c` is built. The stamp is written truthfully (`BUILD_SHARED_LIBS=OFF`), so the
  pin check cannot pass: `DA_IGNORE_SYMENGINE_PIN=ON` turns that into a warning, and the warning
  names only `stamp does not have BUILD_SHARED_LIBS=ON`. A stamp claiming `ON` would be false. This
  library never meets symengine.py, so the shared-library part of the pin does not apply.
  `-DCMAKE_SKIP_BUILD_RPATH=ON`: the build-tree binary is shipped, and without this it carried
  `RUNPATH /usr/local/lib`; the audit sets `$ORIGIN`. The script fails if `nm -D --defined-only`
  shows anything but `mdac_*` and toolchain symbols. Result: 561 `mdac_*` exports and nothing else,
  no SymEngine symbol in the dynamic table (it is linked in), `NEEDED` libgmp.so.10 plus the C/C++
  runtime, and `RUNPATH $ORIGIN`. `GMP_jll` dependency: build against 6.2.1, compat `6.2.1` (Julia
  1.10's stdlib GMP_jll is 6.2.1; 1.11+ ship 6.3.0, the same `libgmp.so.10`). `julia_compat="1.10"`,
  GCC 10.
  - Tools: Julia 1.13.1 (BinaryBuilder runs on it, no other Julia needed), BinaryBuilder 0.6.6,
    BinaryBuilderBase 1.44.0. Runner: `BINARYBUILDER_RUNNER=userns`. The Docker runner built too,
    but it runs `sudo chown` on the build directory after each step, which fails without a terminal
    or passwordless sudo. `expand_cxxstring_abis` now yields only `x86_64-linux-gnu-cxx11` (cxx03
    needs `old_abis=true`), so one tarball (tree hash `9834a1bc…`); the build takes about 5 min.
    `--deploy=local` writes the JLL to the first depot's `dev/MiraDAC_jll` (`JULIA_PKG_DEVDIR` is
    ignored), with a `file://` artifact URL.
  - Julia package: `MiraDAC_jll` is in `[deps]` (compat `1.1`). `libmiradac` is the preference
    when set, else `MiraDAC_jll.libmiradac_c` (the soname constant, which the JLL's `__init__`
    has dlopened by full path), else `""` (the platform has no JLL build), which `__init__` reports.
    The ABI check is unchanged. `runtests.jl` logs the loaded library's path. Until T9.2,
    `MiraDAC_jll` is unregistered, so every environment needs `Pkg.develop` of the local JLL:
    `dev_setup.jl` does it for `julia/MiraDAC`, `bench` and `docs` when `MIRADAC_JLL` is set. On
    Julia >= 1.11, `develop` writes a machine-local `[sources]` path into the tracked
    `Project.toml`, so `dev_setup.jl` restores that file afterwards and only the untracked
    manifest keeps the path; later runs then work without `MIRADAC_JLL`. **The
    CI `julia` job (T0.5) cannot instantiate until T9.2, unless it builds or fetches the JLL.**
  - Tests: (1) a fresh environment with the local JLL, MiraDAC developed, and no `libmiradac`
    preference (the package's `LocalPreferences.toml` moved away, because `Pkg.test` copies it).
    `Pkg.test("MiraDAC")` loads `depot/artifacts/9834a1bc…/lib/libmiradac_c.so` and passes 1228/1228,
    including T5.6 coexistence with SymEngine_jll 0.12.0. (2) `MIRADAC_JLL=… julia julia/dev_setup.jl`,
    then `Pkg.test()` loads `build/capi/libmiradac_c.so` and passes 1228/1228. (3)
    `ctest --test-dir build -R capi` passes 3/3. Julia 1.10 was not run (not installed).
  - For the other platforms: each must build GMP-linked static SymEngine with the toolchain file.
    The export check uses `nm -D`, which needs a Mach-O variant on macOS and an export-table check
    on Windows. `-Wl,--exclude-libs,ALL` and the version script `capi/src/exports.map` are GNU ld
    options, so `capi/CMakeLists.txt` may need a branch for Apple ld (`-exported_symbols_list`) and
    for MinGW. That is a CMake change, to be planned there.
- *As built (second step, `aarch64-apple-darwin` + `x86_64-apple-darwin`, 2026-10-02, branch
  `zhanghe9704/jll-macos`):* the CMake change above is in (`capi/CMakeLists.txt` compiles the C API
  sources as an OBJECT library and links with `-Wl,-exported_symbols_list,<generated list>`; the
  list comes from `nm -U` over the objects, since most of the API is declared through macros in
  `miradac.h` — the header cannot be parsed). The numeric-only cross-check compiled both targets
  with Clang 18/libc++ with zero warnings and no source change; an include-hygiene commit added
  direct `<algorithm>`/`<utility>`/`<type_traits>` includes (libc++ prunes transitive ones). The
  recipe's audit is per format now: `nm -gU` and the `_mdac_` pattern on darwin (the sandbox nm is
  cctools nm: no `-D`, no `--defined-only`), `nm -D --defined-only` and `mdac_` elsewhere. Both
  darwin legs build the full symbolic recipe (static SymEngine + `GMP_jll`, ~5 min each) and pass
  the audit: 561 `mdac_*` exports, nothing else, install name `@rpath/libmiradac_c.dylib`,
  `@rpath/libgmp.10.dylib` the only non-runtime dependency. The darwin builds need
  `BINARYBUILDER_AUTOMATIC_APPLE=true` (Apple SDK terms) and the MiraDAC GitSource at a commit
  with the per-linker capi CMake (after v1.1.0); they were verified with the recipe driven against
  a `DirectorySource` of the branch. Runtime tests are for the machines that can run them.
- *As built (third step, `x86_64-w64-mingw32` + `i686-w64-mingw32`, 2026-10-02, branch
  `zhanghe9704/jll-windows`):* the MinGW CMake branch is in: `capi/CMakeLists.txt` routes WIN32 to
  a `.def` file generated from the compiled objects (`gen_export_list.cmake`, `DEF_FILE` mode,
  plain undecorated names — MinGW ld matches them with or without the i686 leading underscore),
  and passing a `.def` switches MinGW ld from export-everything to exactly the listed names.
  Passing the ELF options through had been worse than a link failure: MinGW ld accepts
  `--version-script` silently on PE, exports nothing, and the DLL ships with an empty export
  table. The DLL links `-static`: libgcc, libstdc++ and winpthread (pulled in by the emutls
  behind `thread_local`) are baked in, so the imports stay KERNEL32/msvcrt plus the GMP DLL from
  `GMP_jll` — the same only-GMP story as the other platforms. The recipe's audit grew a dll
  branch: `nm -D` does not read PE export tables (it reports "no symbols"), so the exported names
  come from `objdump -p` (`[Ordinal/Name Pointer] Table`; works for both bitnesses, the 32-bit
  objdump omits the `+base[...]` column). `-DSymEngine_DIR` points at `<prefix>/CMake` on MinGW,
  where SymEngine installs its CMake package on WIN32 (`lib/cmake/symengine` elsewhere). The
  numeric-only cross-check compiled both targets with MinGW GCC 10.2, zero warnings, 561
  `mdac_*` exports each; both symbolic legs build the full recipe and pass the audit the same
  way. A latent Win64 bug fixed in passing: the double→exact-integer promotion used `long`
  (32-bit on Win64, undefined past 2^31); now `int64_t` (`DAEnv::promote`, the interop
  copy-with-promotion, and the cast in `mdac_expr_new_i`). Tree hashes `882c8c2c…` (x86_64) and
  `81df14c7…` (i686); the Windows tarball holds `bin/libmiradac_c.dll` (the PE RUNTIME
  destination), which `LibraryProduct` finds. Verified with the recipe driven against a
  `DirectorySource` of the branch (the GitSource still pins v1.1.0). Runtime tests are for the
  machines that can run them (`jll-ci-tests`).

**T9.2 `MiraDAC_jll`.** Submit the recipe to Yggdrasil (from the user's GitHub account); once
merged, make `MiraDAC_jll` the default library in `src/MiraDAC.jl` (the preference still
overrides it for development).

**T9.3 Registration.** Register `MiraDAC.jl` in General from the `julia/MiraDAC` subdirectory
(Registrator with `subdir`); package versions follow the C++ project version.

### Deferred (not in this plan)

- CxxWrap.jl bindings (A.2).
- Zero-copy interop with SymEngine.jl (needs the same SymEngine version and one shared library;
  A.1).
- Multi-threaded use of one env (needs locks in the C++ pools).
- Automatic differentiation / ChainRules integration.
