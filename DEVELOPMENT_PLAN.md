# MiraDAC — Development Plan

**Goal:** Combine the numerical DA library (`ref/tpsa`, "NDA") and the symbolic DA library
(`ref/tpsa_sym`, "SDA") into a single C++ library that supports both numerical and symbolic
Differential Algebra, and computation between them, from one unified, templated codebase.

This document has two parts:

1. **Part A — Design** (the agreed architecture).
2. **Part B — Staged implementation plan** (explicit, ordered steps with per-stage acceptance
   criteria, written to be executable by a small model with no further design decisions required).

> **Ground rules for any implementer**
> - The two reference repos live in `ref/tpsa` and `ref/tpsa_sym`. **They are read-only.** Do not
>   modify anything under `ref/`. You may freely copy and adapt code *out* of them.
> - `ref/tpsa` and `ref/tpsa_sym` are the **source of truth for algorithms and for expected
>   results**. When in doubt about behavior, the new code must reproduce the reference output.
> - Every stage must compile and pass its tests before the next stage begins. Do not skip ahead.
> - All design decisions are already made and recorded in Part A. Do not invent new ones. If a
>   genuine ambiguity appears, stop and ask rather than guessing.

---

## Part A — Design

### A.1 Background: why this merge is clean

The two reference libraries descend from the same engine (Lingyun Yang's `tpsa.cpp`, extended by
He Zhang in `tpsa_extend.cc`). They are **logically identical**; they differ in exactly two ways:

1. **Coefficient type.** NDA stores `double` coefficients in its memory pool
   (`double** advecpool`). SDA stores `SymEngine::Expression` (`SymEngine::Expression** advecpool`).
2. **Namespacing.** SDA is wrapped in `namespace SymbTPSA` / `namespace SymbDA`; NDA is in the
   global namespace.

The monomial-indexing math (`order_index`, `base`, the product-index table `prdidx`, the hash `H`,
`FULL_VEC_LEN`) and every arithmetic kernel are byte-for-byte the same in logic. Therefore the
correct unification is to **template the engine on the coefficient type `T`** and instantiate it
for `double` and `SymEngine::Expression`.

### A.2 The central architectural split

The engine has two separable halves. Keeping them separate is what makes the whole design work:

- **Layout layer (type-independent).** Everything that depends only on `(number of variables, order
  scheme)`: the monomial table `base`, `order_index`, the product-index table `prdidx`, the hash
  `H`, and `FULL_VEC_LEN`. This does **not** depend on `T`.
- **Coefficient layer (templated on `T`).** The memory pool that stores coefficients, and the
  arithmetic kernels that read/write them.

Two `DAVector`s of different coefficient type in the *same* environment **share one layout**. This
is the foundation of cheap interop (see A.7).

### A.3 Naming: `da` is the umbrella for both NDA and SDA

- **`da` is the general name** for differential algebra in this library, covering both the
  numerical flavor (NDA) and the symbolic flavor (SDA). Everything lives in **namespace `da`**.
  (The repository is `MiraDAC`; the library, namespace, and headers are all named `da`.)
- General templated type: `da::DAVector<T>`.
- The two flavors are aliases of that one template:
  - `using NDA = da::DAVector<double>;` (numerical DA)
  - `using SDA = da::DAVector<SymEngine::Expression>;` (symbolic DA, when symbolic is enabled)
- Free functions and macros keep the `da_` / `DA_` convention from the reference: `da::da_init`,
  `da::da_clear`, `da::da_composition`, `DA_CHECK_ENV`, `DA_WITH_SYMBOLIC`, etc.
- **Bases accessor.** The reference exposes the bases as a global object named `da` (`da[i]`).
  Because the namespace is now also named `da`, the bases object is renamed to **`base`**, accessed
  as `da::base[i]`. A variable cannot share the name of its enclosing namespace under a
  `using namespace da;` directive, so this single rename versus the reference is deliberate.
- The legacy global `DAVector` / `SymbDA::DAVector` names are **not** preserved; this is a new
  library with a new API. (A thin compatibility shim may be added later if needed, but it is not in
  scope for v1.)

### A.4 Environments (`DAEnv`) and the multiple-environment policy

- All formerly-global engine state is encapsulated in a **`DAEnv`** object. A `DAEnv` owns:
  - one **Layout** (shared by both coefficient types), and
  - one `Pool<double>`, and (when `WITH_SYMBOLIC`) one `Pool<SymEngine::Expression>`.
- **Multiple environments are supported but expected to be rare.** The single-environment path must
  feel and perform exactly like the reference libraries.
- There is always an implicit **default / current environment**, configured by the classic call
  `da::da_init(order, num_vars, pool_size)`. A single-env user never sees a `DAEnv` object and
  never passes an env to anything.
- A `DAVector<T>` stores a `DAEnv*` (the env it belongs to) plus its pool slot index.
  - The `DAEnv*` is set at **construction time** from a thread-local **current env**.
  - Copy/move constructors **inherit the source vector's env** (no ambiguity).
  - Operations **never** consult the current env; they use the env stored in their operands.
- **Cross-environment computation is forbidden.** Operating on two vectors from different envs is an
  error. To combine them, the user must explicitly transfer one into the other's env via
  `env.import(vec)` / `da::convert(vec, env)`. Same-layout transfer is a cheap slot copy; different
  layout is an error in v1.
- Multi-env API surface is intentionally minimal: `da_make_env(...)`, `da_select_env(env)`,
  `da_current_env()`, plus `import`/`convert`. No registry, no scheduling.

### A.5 The environment guard

- Every binary and mixed operation begins with a pointer-equality check
  `lhs.env_ == rhs.env_`; on mismatch it throws `std::logic_error`.
- The check is **compile-time switchable** via `DA_CHECK_ENV` (default **on**). Building with
  `-DDA_CHECK_ENV=0` strips it for users who have profiled a hot numerical loop and want it gone.
- A public `bool da::same_env(const DAVector<T>&, const DAVector<U>&)` is always available so users
  can also assert explicitly.
- **Cost rationale (already decided):** the check is one branch-predicted pointer compare, hoisted
  out of all inner loops; negligible for `double`, parts-per-million for `Expression`. It is kept on
  by default because the failure it prevents (a slot index used against the wrong pool) is silent
  memory corruption.

### A.6 The memory pool: templated RAII `Pool<T>` with an O(1) free-list

- `Pool<T>` is an RAII class owned by `DAEnv`. Its **constructor** allocates one contiguous block of
  `full_len * n` elements and builds the free-list; its **destructor** releases it. This replaces
  the manual `ad_reserve` / `ad_clear` global management and gives deterministic teardown (this is
  expected to fix the SDA "segfault after all tests pass", which is a static-destruction-order
  problem between the global pool and SymEngine globals).
- **Slot model (unchanged from the reference):** creating a `DAVector` pops a free slot
  (O(1) free-list pop, *not* a heap allocation); destroying it pushes the slot back (O(1)). The pool
  is pre-reserved and reused — that is the entire point of the pool.
- **POD fast-path vs non-POD path** — selected with `if constexpr (std::is_trivially_copyable_v<T>)`:
  - `double` (trivially copyable): zero a slot with `memset`, copy with `memcpy`.
  - `SymEngine::Expression` (holds a reference-counted pointer): zero a slot by assigning
    `Expression(0)` per element (releases the RCP); copy with per-element `operator=`. `memset`/
    `memcpy` are **forbidden** for this type — they corrupt refcounts.
- The free-list logic (`adlist`, `ad_flag`/head, `ad_end`/tail) only manipulates slot **indices** and
  is therefore identical for both types — written once, not specialized.

### A.7 Interop: NDA ↔ SDA

- **Same-type operations never convert.** `NDA·NDA → NDA`, `SDA·SDA → SDA`.
- **Mixed `NDA`-with-`SDA` operators are overloaded and return `SDA`.** A mixed op promotes the NDA
  operand to a temporary SDA by a **same-layout slot copy** with `expr[i] = Expression(dbl[i])`
  (no re-indexing, no algebra — linear in length), then calls the existing SDA kernel. This is the
  only automatic conversion, and only in the unavoidable direction (a result with symbolic
  coefficients cannot stay numeric).
- **SDA → NDA is not a conversion; it is evaluation.** An SDA becomes numeric only by substituting
  values for its symbols. Provide `da::evaluate(sda, symbol_values) -> NDA`, built on the existing
  `eval` / `eval_funs` machinery from `ref/tpsa_sym/src/sda.cc`.
- Promotion requires identical layouts (same env → pure copy; different env, same layout → copy +
  retag; different layout → error).

### A.8 Performance model (zero-cost abstraction)

- Templates are **monomorphized**: `Pool<double>` and the `double` kernels compile to essentially the
  same machine code as the current `tpsa.cpp`. No virtual dispatch anywhere in the kernels.
- `if constexpr` resolves the POD/non-POD choice at compile time — no runtime type test.
- The `DAEnv*` indirection is resolved once at the top of each operation into local raw pointers;
  inner loops run on raw pointers exactly as today.
- The `MonomialScheme` abstraction (A.9) is consulted **only at table-build time** (`da_init`); the
  kernels read the precomputed flat `prdidx` table and never call into the scheme at runtime.
- Compile-time/binary-size cost (two instantiations) is contained by **explicit instantiation** in
  `.cpp` files and by the `WITH_SYMBOLIC` flag (numeric-only build instantiates only `double` and
  needs no SymEngine).
- C++ standard: **C++17** (for `if constexpr`). C++20 is acceptable if a toolchain benefit appears,
  but C++17 is the target.

### A.9 Forward-compatibility seams (designed-in, not built in v1)

- **Per-variable order.** The layout layer is built around a `MonomialScheme` abstraction with a
  `UniformOrder` implementation (the only one in v1). A future `PerVariableOrder` scheme (variable
  `i` capped at order `o_i`, optional total-order cap) changes only how `base` / `order_index` /
  `prdidx` are generated and how products are truncated; the kernels are untouched because the
  truncation is baked into `prdidx`. **Not implemented in v1**, but the seam must exist so it is a
  localized addition later.
- **Complex DA.** Templating on `T` means `T = std::complex<double>` could later subsume much of the
  hand-rolled `std::complex<DAVector>` / `cd_composition` code as a third instantiation. Not in v1;
  v1 keeps the existing complex support as-is on the numerical side (see Stage 5 note).

### A.10 Build configuration summary

- `WITH_SYMBOLIC` (CMake option, default `ON`): when `OFF`, SymEngine is not required and only the
  numerical library is built.
- `DA_CHECK_ENV` (compile definition, default `1`): the env guard.
- Library targets: `daShared` (`.so`/`.dll`) and `daStatic` (`.a`). Public umbrella header
  `include/da/da.h`. Tests via Catch2 (v2.13.x header, copied from a reference repo's
  `catch.hpp`).

---

## Part B — Staged implementation plan

### Conventions used in every stage

- **Directory of the new library:** the repository root (this file's directory). Reference code is
  under `ref/`. Create new code under `include/da/`, `src/`, `test/`, `examples/`.
- Each stage lists: **Goal**, **Files**, **Steps**, **Source mapping** (what to port from `ref/`),
  and **Acceptance** (must pass before continuing).
- "Port" means: copy the reference function, adapt it to the new structure (template parameter,
  `DAEnv`/`Layout`/`Pool` access instead of globals, `da` namespace), and keep its algorithm
  unchanged.
- Run tests from the build directory unless a test's comment says otherwise. (The reference numeric
  tests must be run from the `test/` folder; preserve that requirement.)

---

### Stage 0 — Project scaffolding

**Goal:** an empty but correct build system; both libraries configure and link; Catch2 runs a
trivial test. No DA logic yet.

**Files**
- `CMakeLists.txt` (root)
- `include/da/da.h` (umbrella header, empty for now)
- `src/placeholder.cpp` (temporary, compiles to nothing meaningful)
- `test/CMakeLists.txt`, `test/catch.hpp` (copy from `ref/tpsa/test/catch.hpp`),
  `test/catch_main.cc` (copy from `ref/tpsa/test/catch_main.cc`), `test/test_smoke.cc`
- `examples/CMakeLists.txt` (empty target for now)
- `LICENSE` (copy `ref/tpsa/LICENSE.md`)

**Steps**
1. Root `CMakeLists.txt`: `cmake_minimum_required(VERSION 3.16)`, `project(da CXX)`,
   `set(CMAKE_CXX_STANDARD 17)`, `set(CMAKE_CXX_STANDARD_REQUIRED ON)`.
2. Add `option(WITH_SYMBOLIC "Build symbolic (SDA) support" ON)` and
   `set(DA_CHECK_ENV 1 CACHE STRING "Enable env guard")`; add
   `add_compile_definitions(DA_CHECK_ENV=${DA_CHECK_ENV})`.
3. If `WITH_SYMBOLIC`, `find_package(SymEngine REQUIRED CONFIG)` and capture
   `${SYMENGINE_INCLUDE_DIRS}` / `${SYMENGINE_LIBRARIES}` (mirror `ref/tpsa_sym/CMakeLists.txt`).
   Define `DA_WITH_SYMBOLIC` for the compiler when on.
4. Define `daShared` and `daStatic` from the (currently placeholder) sources; set
   `PUBLIC` include dir `include/`. Link SymEngine to both when `WITH_SYMBOLIC`.
5. `test/`: build `run_tests` from `catch_main.cc` + `test_smoke.cc`, link `daStatic`. Add a
   trivial `TEST_CASE("smoke"){ REQUIRE(1+1==2); }`.

**Acceptance**
- `cmake -DWITH_SYMBOLIC=ON .. && make` succeeds.
- `cmake -DWITH_SYMBOLIC=OFF .. && make` succeeds **without SymEngine installed**.
- `run_tests` passes the smoke test.

---

### Stage 1 — Layout layer (type-independent)

**Goal:** the monomial-indexing machinery as a standalone, fully-tested `Layout` class. No
coefficients, no `T`.

**Files**
- `include/da/monomial_scheme.h` — `MonomialScheme` interface + `UniformOrder`.
- `include/da/layout.h`, `src/layout.cpp`.

**Steps**
1. `MonomialScheme` (abstract): describes the admissible exponent set. `UniformOrder{nv, nd}` is the
   only implementation in v1 (total order ≤ `nd`). Provide: `num_vars()`, `max_total_order()`,
   `full_len()` (= `C(nv+nd, nd)` for uniform).
2. `Layout` holds and computes: `gnv`, `gnd`, `FULL_VEC_LEN`, `order_index[]`, `base[]`, `prdidx[][]`,
   `H`. Build them in the `Layout` constructor from a `MonomialScheme`.
3. Port the generation code from `ref/tpsa/src/tpsa.cpp`: `gcd`, `comb_num`, `init_order_index`,
   `init_prod_index`, `init_base`, `choose`, `within_limit`. Convert file-static globals into
   `Layout` members. Keep the algorithms identical.
4. Port `ADOrderTable` (from `ref/tpsa/include/tpsa_extend.h` / `tpsa_extend.cc`) into the layout:
   `order_table` (index → exponent vector) and `order_index`-map (exponent vector → index),
   `generate_order_table()`, `orders(i)`, `find_index(orders)`.
5. Add `change_order(new_order)` / `restore_order()` (port `ad_change_order`/`ad_restore_order` from
   `tpsa_extend.cc`) operating on `Layout` members.

**Source mapping:** `ref/tpsa/src/tpsa.cpp` (lines defining `comb_num`, `init_order_index`,
`init_prod_index`, `init_base`); `ref/tpsa/src/tpsa_extend.cc` (`ADOrderTable`, `ad_change_order`,
`ad_restore_order`, `ad_generate_order_table`).

**Acceptance (new unit tests in `test/test_layout.cc`)**
- For `(nv=3, nd=4)`: `full_len()==35`; for `(nv=6, nd=6)`: `full_len()==924`;
  for `(nv=6, nd=10)`: `full_len()==8008`. (Values cross-checked against `ref/tpsa` README Table 1.)
- `order_index` is monotonically increasing and `order_index[nd+1]==full_len()`.
- Round-trip: for every index `i` in `[0, full_len())`, `find_index(orders(i)) == i`.
- Spot-check `prdidx`: pick simple monomials (e.g. `x*y`) and verify the product index equals the
  index of the summed exponent vector via `find_index`.

---

### Stage 2 — Templated pool `Pool<T>`

**Goal:** the RAII pool with the O(1) free-list and POD/non-POD slot primitives, tested for
`double` (and a non-trivial proxy type).

**Files**
- `include/da/pool.h` (template definition; small enough to be header-only, or paired with explicit
  instantiation in `src/pool.cpp`).

**Steps**
1. `template<class T> class Pool` with members: `T* block_`, `std::vector<T*> slot_`,
   `std::vector<unsigned> len_` (= `adveclen`), free-list `std::vector<unsigned> free_`,
   `unsigned head_, tail_, full_len_, size_`.
2. Rule-of-five: non-copyable, movable; destructor frees `block_`.
3. `reserve(full_len, n)`: allocate one block of `full_len*n` `T`, wire `slot_[i]`, build free-list,
   zero lengths. Port from `ref/tpsa/src/tpsa_extend.cc::ad_reserve`.
4. `assign() -> unsigned`: free-list pop (port `ad_assign`). `alloc() -> unsigned`: `assign()` then
   set const 0, length 1 (port `ad_alloc` semantics). `free(i)`: push slot to tail (port `ad_free`).
5. `reset(i)`, `clean(i, eps)`, `pool_clean(idx)`: port from `tpsa_extend.cc`.
6. Slot primitives with `if constexpr (std::is_trivially_copyable_v<T>)`:
   `zero_slot(T*)`, `copy_slot(const T*, T*, len)` — `memset`/`memcpy` vs typed loops (A.6).
7. Optionally over-align `block_` (e.g. 32 bytes) for the trivially-copyable path to aid
   auto-vectorization.

**Acceptance (`test/test_pool.cc`, instantiate `Pool<double>` and `Pool<std::string>` as a non-POD
proxy so this stage does not require SymEngine)**
- `reserve(35, 100)` then `assign` 100 slots; the 101st `assign` reports pool exhaustion (matches
  reference behavior: prints "Run out of vectors" — here throw or return sentinel; pick throw and
  document it).
- `free` then `assign` returns a slot (round-trip); lengths reset to 0 on assign.
- `zero_slot`/`copy_slot` correct for both `double` and `std::string` (non-POD path must not corrupt).
- Move-construct a pool; original is emptied, moved-to works; no double free under ASan.

---

### Stage 3 — `DAEnv` and environment management

**Goal:** assemble Layout + pools into an environment; implement the default/current-env model.

**Files**
- `include/da/env.h`, `src/env.cpp`.

**Steps**
1. `class DAEnv` owns: `Layout layout_`, `Pool<double> pool_d_`, and (under `DA_WITH_SYMBOLIC`)
   `Pool<SymEngine::Expression> pool_e_`.
2. `template<class T> Pool<T>& pool();` with specializations returning `pool_d_` / `pool_e_`.
3. Construction: `DAEnv(order, num_vars, pool_size, bool table)` builds `layout_` from a
   `UniformOrder` scheme and reserves both pools to `pool_size`.
4. Thread-local current-env pointer + `da_current_env()`, `da_select_env(DAEnv&)`,
   `DAEnv& da_make_env(order, nv, pool, table)`.
5. Classic entry points: `da_init(order, nv, pool, table=...)` creates/configures the **default**
   env and selects it; `da_clear()` tears it down. Match the reference defaults: numeric
   `table=false`, symbolic `table=true` (see `ref/tpsa/src/da.cc` and `ref/tpsa_sym/src/sda.cc`
   `da_init`).
6. `env.import<T_from, T_to>(slot)` / free function `convert` per A.7 (same-layout copy; promotion
   double→Expression when crossing types). Cross-layout → throw.
7. `same_env(...)` and the `DA_CHECK_ENV`-gated guard helper (a small inline function/macro used by
   all binary ops).

**Acceptance (`test/test_env.cc`)**
- `da_init(4,3,100)`; `da_current_env()` returns it; layout `full_len()==35`.
- `da_make_env` + `da_select_env` switch the construction target; switching back restores.
- Tear down and re-init in the same process works (no leaks under ASan).

---

### Stage 4 — Numerical kernels (`T = double`)

**Goal:** all arithmetic kernels, templated, instantiated and validated for `double`. Reproduce the
numerical engine exactly.

**Files**
- `include/da/engine.h` (templated kernel declarations in `namespace da::detail`).
- `src/engine.cpp` (definitions + `template ... <double>` explicit instantiation).

**Steps**
1. Declare kernels as `template<class T>` free functions taking `(const Layout&, Pool<T>&, slots…)`:
   `add`, `sub`, `mult`, `mult_const`, `add_const`, `div`, `div_c`, `c_div`, `copy`, `reset`,
   `set_const`, `pok`, `pek`, `elem`, `var`, `subst`/`composition`, `derivative`, `integrate`,
   `truncate`, `clean`, `n_element`, `zero_check`, `norm`, `weighted_norm`.
2. Port each body from `ref/tpsa/src/tpsa.cpp` and `ref/tpsa/src/tpsa_extend.cc`, replacing global
   `advec[iv]`, `adveclen[iv]`, `base`, `prdidx`, `gnv`, `gnd`, `FULL_VEC_LEN` with `pool.slot(iv)`,
   `pool.len(iv)`, `layout.base()`, `layout.prdidx()`, `layout.num_vars()`, etc. **Resolve these into
   local raw pointers at the top of each kernel; keep inner loops on raw pointers.**
3. Keep the composition / group-composition routines (`ad_composition` overloads in
   `tpsa_extend.cc`) — these are the performance-critical paths from the README tables.
4. Explicitly instantiate the `double` versions at the bottom of `engine.cpp`.

**Source mapping:** `ref/tpsa/src/tpsa.cpp` (all `ad_*` kernels: `ad_add`, `ad_sub`, `ad_mult`,
`ad_mult_const`, `ad_add_const`, `ad_div`, `ad_div_c`, `ad_c_div`, `ad_copy`, `ad_reset`, `ad_const`,
`ad_pok`, `ad_pek`, `ad_elem`, `ad_var`, `ad_subst`, `ad_sqrt`, `ad_exp`, `ad_log`, `ad_sin`,
`ad_cos`, `ad_inverse`) and `ref/tpsa/src/tpsa_extend.cc` (composition, substitute, derivative,
integrate, norm, clean, n_element, zero_check).

**Acceptance**
- Port the numerical reference tests from `ref/tpsa/test/tests.cc` (3 `TEST_CASE`s) into
  `test/test_numeric_engine.cc`, rewritten against the new API, and they pass.
- Add a direct equivalence check: for a handful of random low-order vectors, the new `mult`/`add`/
  `composition` results equal results computed by linking the reference `ref/tpsa` library (or
  match values hard-coded from running the reference). Cross-check at least one composition against
  README Table 1 dimensions (order 4, 210 terms; order 6, 924 terms).

---

### Stage 5 — `DAVector<double>` wrapper, operators, math functions (NDA complete)

**Goal:** a complete, shippable **numerical** library (`WITH_SYMBOLIC=OFF` builds and works fully).

**Files**
- `include/da/davector.h` (`template<class T> struct DAVector`), `src/davector.cpp`.
- `include/da/functions.h`, `src/functions.cpp` (math functions for `double`).
- `include/da/base.h` (the `Base`/`da[i]` bases), in the current env.

**Steps**
1. Implement `DAVector<T>` members from `ref/tpsa/include/da.h` / `src/da.cc`: constructors
   (default, copy, move, from `double`, from `int`, from `std::vector`), `=`, `+=`, `-=`, `*=`, `/=`,
   `print`, `con`, `length`, `element`/`set_element`, `derivative`, `reset`, `clean`, `iszero`,
   `norm`, `to_vector`, static `dim/order/full_length`, destructor (returns the slot via
   `pool.free`).
   - Each `DAVector` stores `DAEnv* env_` (set from current env at construction) and `unsigned slot_`.
   - Copy/move inherit `env_` from the source.
2. Free operators (`+ - * /` in all combinations with `double`) and unary `+/-`, `==`. Each binary
   op runs the `DA_CHECK_ENV` guard first.
3. Math functions (`sqrt, exp, log, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh,
   atanh, pow, abs, erf, atan2`) — port from `ref/tpsa/src/da.cc`.
4. `Base` + the convenience `da::base[i]` accessor bound to the current env (port `Base` from
   `ref/tpsa/include/da.h`, renaming the global object from `da` to `base` per A.3).
5. **Complex DA (numerical):** port the existing `std::complex<DAVector>` operators and
   `cd_composition` from `ref/tpsa` as-is (kept, not refactored — see A.9). Keep them behind the
   numerical build.
6. Top-level `da_*` free functions mirroring `ref/tpsa/include/da.h` (`da_der`, `da_int`,
   `da_substitute`, `da_composition`, `da_count`, `da_remain`, `da_poolsize`, `da_full_length`,
   `da_set_eps`, etc.) delegating to the current env.

**Acceptance**
- `WITH_SYMBOLIC=OFF` build is fully functional and links **without SymEngine**.
- Port and pass all numerical tests (`ref/tpsa/test/tests.cc`).
- Build and run the numerical examples (`ref/tpsa/examples/examples.cc`,
  `example_complex_da.cc`) adapted to the new API; output matches the reference.

---

### Stage 6 — Symbolic enablement (`T = SymEngine::Expression`)

**Goal:** instantiate the same engine for `Expression`; reach parity with `ref/tpsa_sym`.

**Files**
- `include/da/symbolic_ops.h` (the `Expression` scalar operator helpers).
- Extend `engine.cpp`, `davector.cpp`, `functions.cpp` with `Expression` explicit instantiation.
- `src/symbolic.cpp` (`is_zero`, `simplified_expr`).

**Steps**
1. Guard all symbolic code with `#ifdef DA_WITH_SYMBOLIC`.
2. Port the scalar operator overloads for `SymEngine::Expression` from
   `ref/tpsa_sym/include/symbolic.h` into `symbolic_ops.h`.
3. Add the `Expression`-specific `DAVector` API from `ref/tpsa_sym/include/sda.h`/`src/sda.cc`:
   `DAVector(SymEngine::Expression)`, `con()` returning `Expression`, `set_element(Expression)`,
   `reset_const(Expression)`, `=/+=/-=/*=//=` with `Expression`, `simplify()`, `eval`, `eval_funs`.
4. Confirm the kernels (Stage 4) instantiate for `Expression` unchanged — the only differences are
   the `Pool` non-POD path (already handled in Stage 2) and the `Expression` scalar ops. The
   reference SDA fork only changed coefficient type and a few signatures
   (`ad_pok`, `ad_mult_const`, `ad_add_const`, `ad_elem`, `ad_pek`); verify these map to the typed
   kernel paths. Source: diff of `ref/tpsa/src/tpsa.cpp` vs `ref/tpsa_sym/src/tpsa.cpp`.
5. Explicitly instantiate `Pool<Expression>`, the kernels, `DAVector<Expression>`, and the math
   functions for `Expression` (the subset present in `ref/tpsa_sym`: `sqrt, exp, log, sin, cos, tan,
   asin, acos, atan, sinh, cosh, tanh, pow, erf` — see `sda.h`).
6. `SDA` alias and `da::base[i]` bases working in symbolic mode.

**Acceptance**
- `WITH_SYMBOLIC=ON` build links against SymEngine.
- Port and pass the symbolic tests (`ref/tpsa_sym/tests/tests.cc`).
- Build/run symbolic examples (`example_1_fundamental.cc`, `example_2_eval.cc`,
  `example_3_eval_funs.cc`) adapted to the new API; output matches the reference (e.g. the
  `sqrt(da1)` expansion in the `ref/tpsa_sym` README).
- **Teardown check:** run the full symbolic test binary under ASan; confirm there is **no
  segfault on exit** (the RAII pool should resolve the reference's known issue). If a fault remains,
  investigate static-destruction order before proceeding.

---

### Stage 7 — Interop (mixed NDA ↔ SDA)

**Goal:** natural mixed-type expressions and evaluation, per A.7.

**Files**
- `include/da/interop.h`, `src/interop.cpp` (guarded by `DA_WITH_SYMBOLIC`).

**Steps**
1. Implement `promote(const NDA&) -> SDA`: allocate an `Expression` slot in the NDA's env, copy with
   `expr[i] = Expression(dbl[i])`, copy length, retain env. (Same-layout slot copy.)
2. Overload mixed operators: `operator+(-,*,/)` for `(NDA, SDA)` and `(SDA, NDA)`, each running the
   env guard, promoting the NDA side, and returning `SDA`.
3. `evaluate(const SDA&, const std::map<symbol,double>&) -> NDA` using `eval`/`eval_funs` from
   Stage 6; allocate the result in a numerical pool of the same env.
4. `import`/`convert` (Stage 3) finalized for the double↔Expression and cross-env same-layout cases.

**Acceptance (`test/test_interop.cc`)**
- For random low-order vectors: `nda * sda` equals `promote(nda) * sda` coefficient-by-coefficient.
- `evaluate(promote(nda), {})` reproduces `nda` (numeric round-trip within `eps`).
- `evaluate(sda, values)` matches direct numeric computation for a known closed-form case.
- Mixed op between vectors of two **different** envs throws when `DA_CHECK_ENV=1`.

---

### Stage 8 — Multiple environments

**Goal:** validate the multi-env policy and the guard end-to-end.

**Files**
- `test/test_multienv.cc`; possibly small additions to `env.cpp`.

**Steps**
1. Confirm two `DAEnv`s with different `(order, nv)` coexist; vectors built in each carry the right
   `env_`.
2. Confirm `import`/`convert` moves a vector between two **same-layout** envs (copy + retag) and
   that the result computes correctly in the target env.
3. Confirm cross-env operations throw (guard on) and that a different-layout `import` throws.

**Acceptance**
- All multi-env tests pass with `DA_CHECK_ENV=1`.
- Re-run the **entire** numeric + symbolic + interop suite with `-DDA_CHECK_ENV=0`; everything except
  the guard-specific tests passes (those are skipped/compiled out cleanly).

---

### Stage 9 — Packaging, docs, benchmarks

**Goal:** make it installable, documented, and verify performance parity.

**Steps**
1. CMake `install()` rules for `daShared`, `daStatic`, and `include/da/` headers; export a
   `daConfig.cmake`. Provide `WITH_SYMBOLIC` install variants.
2. `README.md`: overview, build (numeric-only and full), quick-start for NDA, SDA, and mixed usage;
   carry over the relevant guidance from both reference READMEs.
3. Doxygen config (adapt from the reference repos' `.github/workflows/doxygen.yml`).
4. New interop example: a single program that builds an NDA and an SDA, combines them, and
   `evaluate`s the result back to numbers.
5. **Benchmark:** reproduce `ref/tpsa` README Table 1 (single composition, orders 2–10) with the new
   `double` path; confirm timings are within noise of the reference (no regression from templating).

**Acceptance**
- `make install` works for both build variants.
- Benchmark numbers are on par with `ref/tpsa` (no measurable regression).
- Docs build; interop example runs and prints expected output.

---

## Appendix — API mapping cheat-sheet (reference → new)

| Reference (NDA, global) | Reference (SDA, `SymbDA::`) | New (`da`) |
|---|---|---|
| `DAVector` | `SymbDA::DAVector` | `da::DAVector<double>` (`da::NDA`) / `da::DAVector<Expression>` (`da::SDA`) |
| `da_init(o,nv,np)` | `SymbDA::da_init(...)` | `da::da_init(o,nv,np)` (configures default env) |
| `da_clear()` | `SymbDA::da_clear()` | `da::da_clear()` |
| `da[i]` (global `Base da`) | `SymbDA::da[i]` | `da::base[i]` (current env) |
| `ad_*` engine fns (globals) | same, `Expression` | `da::detail::*<T>(Layout&, Pool<T>&, …)` |
| `advecpool/advec/adveclen` | same, `Expression` | `da::Pool<T>` members |
| `gnv/gnd/FULL_VEC_LEN/base/prdidx` | same | `da::Layout` members |
| (n/a) | `eval` / `eval_funs` | `da::evaluate(SDA, values) -> NDA` |
| (n/a) | (n/a) | `da::convert` / `DAEnv::import`, mixed operators |

**Do-not-break checklist for implementers**
- Inner kernel loops must operate on raw local pointers (no `env_`/`vector` indirection inside loops).
- Never `memset`/`memcpy` an `Expression` slot.
- Never operate across envs without an explicit `import`/`convert`.
- Keep `ref/` untouched; always diff your behavior against it.
