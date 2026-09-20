# Development Plan — Symbolic Complex DA (Approach A: `std::complex<SDA>`)

**Goal.** Add complex Differential Algebra over *symbolic* coefficients, structured as a
**real part + imaginary part pair of symbolic DA vectors** — i.e. `std::complex<da::SDA>`,
mirroring the existing numerical `std::complex<da::NDA>`. This lets users separate and manipulate
the real and imaginary DA maps independently, exactly as with the numerical complex DA.

**Chosen approach (A).** Generalize the existing, already-tested numerical complex-DA layer from the
concrete type `NDA` (= `DAVector<double>`) to the template `DAVector<T>`, and instantiate it for
`T = SymEngine::Expression` (the `SDA` type). We do **not** use SymEngine's built-in complex numbers
(that was "Approach B", rejected because it cannot separate real/imag DA maps).

> **For the implementer (read this first)**
> - Namespace is `da`. `NDA = DAVector<double>`, `SDA = DAVector<SymEngine::Expression>`.
> - The numerical complex layer is the **template to copy from**; the reference project
>   `ref/tpsa` (READ-ONLY) is the ultimate source of truth for the algorithms.
> - **Never break the numerical build or numerical complex DA.** After every stage, the
>   `WITH_SYMBOLIC=OFF` suite must still pass (currently **696 assertions**), and the numerical
>   complex tests (`[numeric_cd]`) must still pass.
> - All new symbolic code must be inside `#ifdef DA_WITH_SYMBOLIC`. The OFF build must not reference
>   SymEngine.
> - Work stage by stage; each stage must build and pass its tests before the next.
> - `ref/` is read-only. Do not modify it.

---

## Where the numerical complex-DA code lives (what you will generalize)

| Piece | File(s) | Notes |
|---|---|---|
| `get_real`/`get_imag` (reinterpret trick), `cd_copy` | `include/da/da.h` | Type-agnostic already; just template. |
| Complex operators `+ - * /` (with `double`, `complex<double>`, `complex<NDA>`) | `include/da/da.h` (inline) | ~40 inline functions concrete on `NDA`. |
| Complex transcendental fns `exp,sqrt,log,asin,acos,atan,asinh,acosh,atanh,pow` | `include/da/functions.h` (decl) + `src/functions.cpp` (def) | Implemented via real/imag DA arithmetic using the scalar DA functions. |
| `abs(complex<NDA>)` | `include/da/functions.{h,cpp}` | Returns `double` for numerical — **type-specific**, see Stage C2. |
| `cd_composition` (3 overloads) + `complex_da_pow_int_pos` | `include/da/da.h` (decl) + `src/base.cpp` (def) | Composition of complex-DA maps. |
| `read_cd_from_file` / `compare_cd_with_file` / `compare_cd_vectors` | `include/da/da.h` + `src/base.cpp` | Numerical only (compares to numeric files). Do NOT generalize; symbolic uses a different check (Stage C5). |
| Reference algorithms | `ref/tpsa/include/da.h`, `ref/tpsa/src/da.cc` | Authoritative source; the numerical files above were ported from here. |

**Key fact that makes A cheap:** the complex transcendental functions are built entirely from
real/imag **DA-vector** arithmetic and the scalar DA functions (`sqrt,exp,log,sin,cos,tan,asin,
acos,atan,sinh,cosh,tanh,pow`). `SDA` already provides all of those. So generalizing over `T`
and instantiating for `Expression` requires no new math — only re-typing and instantiation.

---

## Build & test commands (use these in every stage)

Portable (on any machine with SymEngine, e.g. your Linux/WSL):
```
./build.sh --symbolic --symengine-dir <PREFIX>/lib/cmake/symengine --gmp-dir <PREFIX> --build-dir build-on
# numerical regression:
./build.sh --build-dir build-off
```

Exact form inside the current dev sandbox (SymEngine built at /tmp/deps, GMP at /tmp/root):
```
export PATH=/sessions/ecstatic-busy-fermat/.local/bin:$PATH
# ON build:
cmake -S <repo> -B /tmp/da_on -DWITH_SYMBOLIC=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DSymEngine_DIR=/tmp/deps/lib/cmake/symengine \
  -DCMAKE_INCLUDE_PATH=/tmp/root/usr/include/x86_64-linux-gnu \
  -DCMAKE_LIBRARY_PATH=/tmp/root/usr/lib/x86_64-linux-gnu
cmake --build /tmp/da_on -j4      # resume if a 45s call is killed (build dir persists)
cd <repo>/test && LD_LIBRARY_PATH=/tmp/deps/lib:/tmp/root/usr/lib/x86_64-linux-gnu /tmp/da_on/test/run_tests
# OFF regression:
cmake -S <repo> -B /tmp/da_off -DWITH_SYMBOLIC=OFF && cmake --build /tmp/da_off -j4
cd <repo>/test && /tmp/da_off/test/run_tests
```

> **Mount gotcha (sandbox only):** editing a file may not bump its mtime, so `cmake` can skip
> recompiling it. After editing a `.cpp`/`.h`, force a rebuild of the affected translation unit:
> `find <build> -name '<file>.cpp.o' -delete` then rebuild. Run tests from the repo `test/`
> directory (reference data files live there).

---

## Type-policy decisions (fixed — do not re-litigate)

1. **Complex value type:** `std::complex<DAVector<T>>`. Aliases to add:
   `using CNDA = std::complex<NDA>;` and (guarded) `using CSDA = std::complex<SDA>;`.
2. **Complex scalar constants** stay `std::complex<double>` (numeric constants promoted into the
   DA). Symbolic scalar constants enter through the real/imag `SDA` parts (which may hold
   `SymEngine::Expression` symbols). Do **not** introduce SymEngine's `I`.
3. **`abs`:** keep the numerical `double abs(const std::complex<NDA>&)` unchanged. For symbolic,
   add `SDA abs(const std::complex<SDA>&)` returning `sqrt(re*re + im*im)` (a magnitude DA). Do NOT
   templatize `abs` into one function; provide the `SDA` overload separately (guarded).
4. **`std::complex<T>` with non-float `T`** is technically outside the C++ standard but works under
   libstdc++ for the arithmetic we use. The numerical side already relies on this; we are extending
   the same pattern, not adding new risk. Keep using our own `da::` complex operators/functions
   rather than `std::exp`/`std::sqrt` (those are only defined for float/double/long double).
5. **File-based comparison** (`compare_cd_with_file`) is numeric-only. Symbolic correctness is
   checked by *evaluation* (Stage C5), not file compares.

---

## Stage C0 — Preparation & scoping (no code changes)

**Goal:** confirm prerequisites and record the exact inventory to generalize.

**Steps**
1. Build and run BOTH suites (commands above). Record baselines: OFF = 696 assertions; ON = current
   symbolic count. These must never regress.
2. In `include/da/functions.h`, confirm `SDA` (Expression) already has scalar
   `sqrt,exp,log,sin,cos,tan,asin,acos,atan,sinh,cosh,tanh,pow` (it does — instantiated in
   `src/functions.cpp`). List any missing; if a complex function needs a scalar function not present
   for `Expression`, note it (none expected).
3. List every `std::complex<NDA>` symbol in `include/da/da.h`, `include/da/functions.h`,
   `src/functions.cpp`, `src/base.cpp` (grep `complex<NDA>`). This is your work checklist.

**Acceptance:** both baselines reproduced; checklist written into the PR/commit message.

---

## Stage C1 — Templatize the complex operator layer (header-only)

**Goal:** make the complex `get_real/get_imag/cd_copy` and all `+ - * /` operators work for
`std::complex<DAVector<T>>` for both `T=double` and `T=Expression`, with no behavior change for
numerical.

**Files:** `include/da/da.h`.

**Steps**
1. Convert `get_real`, `get_imag` (both const and non-const) to
   `template<class T> DAVector<T>& get_real(std::complex<DAVector<T>>&)` etc. (The reinterpret-cast
   body is unchanged; it is type-agnostic.)
2. Convert `cd_copy` overloads to templates over `DAVector<T>` (keep the `complex<double>` and
   `double` source overloads).
3. Convert every complex operator (`+ - * /`) currently written on `NDA` to
   `template<class T>` over `std::complex<DAVector<T>>` and `DAVector<T>`. Keep the three groups:
   (DAVector ⊗ complex<double>), (complex<DAVector> ⊗ double), (complex<DAVector> ⊗ complex<double>),
   (complex<DAVector> ⊗ complex<DAVector>), and DAVector⊗complex<double>. The bodies are unchanged —
   only `NDA` → `DAVector<T>`.
4. Because these are header-inline templates, no explicit instantiation is needed; they compile on
   demand for whichever `T` is used.
5. Watch for ambiguity/overload clashes with `SDA`'s own scalar operators. If any arise, constrain
   the complex templates (e.g. SFINAE `std::enable_if` that `T` is `double` or `Expression`, or a
   small trait `is_da_coeff<T>`), and document it. Prefer minimal constraints.

**Acceptance**
- OFF build: 696 assertions still pass (numerical complex operators unchanged in behavior).
- ON build: still passes.
- Add a tiny compile-only check in `test/test_symbolic_cd.cc` (created here, guarded) that
  constructs `CSDA z(re, im);` and computes `z+z, z-z, z*z, z/CSDA(...)` and reads back
  `get_real(w).con()` — just to prove the operator templates instantiate for `Expression`.

---

## Stage C2 — Templatize complex transcendental functions

**Goal:** `exp, sqrt, log, asin, acos, atan, asinh, acosh, atanh, pow` for `std::complex<DAVector<T>>`,
instantiated for `double` and `Expression`; plus symbolic `abs`.

**Files:** `include/da/functions.h`, `src/functions.cpp`.

**Steps**
1. In `functions.h`, change the declarations from `std::complex<NDA>` to
   `template<class T> std::complex<DAVector<T>> exp(const std::complex<DAVector<T>>&);` etc. (Keep
   them declared in the header; define in the `.cpp` and explicitly instantiate — do not make them
   header-inline, to preserve compile times and match the existing pattern.)
2. In `functions.cpp`, change each definition body `NDA` → `DAVector<T>` and wrap as
   `template<class T> ...`. The algorithms (real/imag decomposition using scalar DA `sin/cos/exp/
   sqrt/log/atan/...`) are unchanged.
3. At the bottom of `functions.cpp`, add explicit instantiations:
   - Always: `template std::complex<DAVector<double>> exp(...);` … for all listed functions.
   - Guarded `#ifdef DA_WITH_SYMBOLIC`: the same list for `SymEngine::Expression`.
4. `abs`: keep `double abs(const std::complex<NDA>&)` exactly as-is. Add (guarded) a separate
   `SDA abs(const std::complex<SDA>&)` returning `sqrt(get_real(v)*get_real(v) +
   get_imag(v)*get_imag(v))`. Declare in `functions.h` (guarded), define+nothing-to-instantiate
   (concrete function) in `functions.cpp` (guarded).
5. `pow(complex, int)` and `pow(complex, double)`: templatize likewise; the `double`-exponent
   version uses `exp(order*log(v))`, which is fine for symbolic.

**Acceptance**
- OFF: 696 pass (numerical complex functions unchanged).
- ON: builds; extend `test/test_symbolic_cd.cc` to call `exp/sqrt/log` on a `CSDA` and confirm it
  compiles and runs (value checks come in Stage C5).

---

## Stage C3 — Templatize `cd_composition` and `complex_da_pow_int_pos`

**Goal:** composition of complex symbolic DA maps.

**Files:** `include/da/da.h` (declarations), `src/base.cpp` (definitions).

**Steps**
1. Generalize `complex_da_pow_int_pos` and the `cd_composition` overloads that take
   `std::complex<DAVector>` and `DAVector` from `NDA` to `DAVector<T>` (`template<class T>`).
   The bodies use the engine composition (`da::detail`/`ad_composition`) and complex-DA arithmetic,
   all of which are already generic/available for `Expression`.
2. Explicitly instantiate the composition functions for `double` (always) and `Expression`
   (guarded), at the bottom of `base.cpp`.
3. The overload that composes with **numeric** `complex<double>` values maps to *evaluation* for
   symbolic and is out of scope here — do NOT instantiate a `complex<double>`-number composition for
   `Expression`. Only the map-into-map composition (`complex<DAVector<T>>` into `complex<DAVector<T>>`,
   and `DAVector<T>` into `complex<DAVector<T>>`) is generalized.
4. Do NOT generalize `read_cd_from_file`/`compare_cd_*` — leave them `NDA`-only.

**Acceptance**
- OFF: 696 pass (numerical `cd_composition` unchanged).
- ON: builds; add a compile+run smoke of `cd_composition` on small `CSDA` maps in the test file.

---

## Stage C4 — Public API surface & aliases

**Goal:** make symbolic complex DA ergonomic and discoverable.

**Files:** `include/da/da.h` (aliases + guarded re-exports), optionally `include/da/interop.h`.

**Steps**
1. Add `using CNDA = std::complex<NDA>;` (always) and, guarded, `using CSDA = std::complex<SDA>;`.
2. Ensure the umbrella header `include/da/da.h` exposes the complex operators/functions for `SDA`
   (they are templates/instantiated, so just make sure nothing is accidentally `#ifdef`-excluded).
3. (Optional, recommended) In `interop.h` (guarded): add
   `CSDA promote(const CNDA&)` — promote a numerical complex DA to symbolic by promoting each part
   (reuse the existing `promote(NDA)->SDA`), and
   `CNDA evaluate(const CSDA&, const std::map<...>& values)` — evaluate each part via the existing
   `evaluate(SDA,...)->NDA`. These make numerical↔symbolic complex interop natural and are the basis
   of the Stage C5 cross-check.

**Acceptance:** OFF 696 pass; ON builds; `CSDA` usable in a user program.

---

## Stage C5 — Tests (authoritative validation)

**Goal:** prove symbolic complex DA is correct, using the numerical complex DA as ground truth via
evaluation.

**File:** `test/test_symbolic_cd.cc` (entire body guarded by `#ifdef DA_WITH_SYMBOLIC`; contributes
nothing to the OFF build). Wire into the test target (it already compiles to nothing when OFF).

**Test design**
1. **Self-initialize** the DA environment at the top of each TEST_CASE (`da::da_init(order,nv,pool)`)
   so cases run in isolation (follow the fix applied to the numerical CD test).
2. **Structural checks:** build `CSDA z(re, im)` with symbolic `re,im` (using `SymEngine::Expression`
   symbols and `da::base[i]`); verify `get_real`/`get_imag` return the expected parts; check a few
   coefficients via `.con()` / `element(...)` and `simplify()`.
3. **Evaluation cross-check (the important one):** for each operation under test
   (`+,-,*,/,exp,sqrt,log,asin,acos,atan,pow`, and `cd_composition`):
   a. Build symbolic inputs `CSDA` whose coefficients are symbols `a,b,c,...`.
   b. Compute the symbolic result `Rs = op(...)`.
   c. Pick concrete numeric values for the symbols; `evaluate` `get_real(Rs)` and `get_imag(Rs)` to
      `NDA` (via `evaluate(SDA,map)`), forming a numerical `CNDA` `R_from_sym`.
   d. Independently compute the SAME operation numerically: substitute the numeric values into the
      inputs to get numerical `CNDA` inputs, run the numerical `op`, giving `R_num`.
   e. `REQUIRE` that `R_from_sym` and `R_num` agree coefficient-by-coefficient within `1e-12`
      (use `compare_cd_vectors` on the numerical results, or compare real/imag DA coefficients).
   This validates the symbolic implementation against the already-trusted numerical one for arbitrary
   symbol values.
4. Cover: all operators; each transcendental function; `abs` (symbolic magnitude vs numerical
   magnitude after evaluation); at least one `cd_composition` case.

**Acceptance**
- ON build: all `test_symbolic_cd` assertions pass; total ON count increases accordingly.
- OFF build: unchanged (696), because the whole file is guarded out.
- Run under ASan (ON) if feasible on the target machine: no leaks/UB across the complex-DA
  destructors and the reinterpret-based `get_real/get_imag`.

---

## Stage C6 — Example, docs, packaging

**Goal:** ship it usefully.

**Steps**
1. `examples/example_complex_symbolic.cc` (guarded): construct a `CSDA` map, do complex arithmetic +
   a function (e.g. `sqrt`), separate and print `get_real(...)`/`get_imag(...)` with `.print()`,
   then `evaluate` to numbers. Wire into `examples/CMakeLists.txt` under `WITH_SYMBOLIC=ON`.
2. `README.md`: add a "Symbolic complex DA" subsection showing `CSDA`, real/imag separation, and the
   evaluate-to-numeric workflow. State clearly that symbolic complex is the `std::complex<SDA>`
   (real+imag pair) model.
3. Doxygen: the new templates/instantiations are picked up automatically; ensure the guarded symbols
   are documented (brief `@brief` on the `CSDA` alias and `abs`/`evaluate` complex overloads).
4. No CMake target changes beyond adding the example and the test file (both already handled by the
   existing `WITH_SYMBOLIC` wiring).

**Acceptance**
- `./build.sh --symbolic ...` builds libs + tests + `example_complex_symbolic`, and the example runs
  and prints sensible real/imag parts.
- `./build.sh` (OFF) unaffected: 696 pass.

---

## Risks & mitigations (tell the implementer)

- **Overload ambiguity** when templatizing operators (Stage C1): if the compiler complains about
  ambiguous `operator+`/`*` between the complex templates and `SDA` scalar operators, constrain the
  complex templates with a small `is_da_coeff<T>` trait (`double` or `Expression` only). Add the
  constraint only where needed.
- **`std::complex<T>` UB pedantry:** compilers may warn under `-std=c++17 -pedantic`. This mirrors the
  existing numerical code; keep our own `da::` complex functions (never `std::exp` on
  `complex<DAVector>`). Do not "fix" by switching to `std::` math.
- **Symbolic blow-up / slow simplify:** complex functions multiply expression size. Keep test cases
  low order (order 3–4, 2–3 vars). Use `simplify()` only where needed for assertions.
- **Regression discipline:** re-run the OFF suite (696) after EVERY stage. If it drops, you changed
  numerical behavior — revert and re-do the generalization without touching the `double` path.
- **Do not** generalize `read_cd_from_file`/`compare_cd_with_file` to symbolic; symbolic validation
  is by evaluation (Stage C5), not file compare.

## Definition of done

- `WITH_SYMBOLIC=ON`: `std::complex<SDA>` supports `+ - * /` (with double/complex<double>/complex<SDA>),
  `exp sqrt log asin acos atan asinh acosh atanh pow abs`, and `cd_composition`; all validated by the
  evaluation cross-check against numerical complex DA.
- `WITH_SYMBOLIC=OFF`: byte-for-byte the same behavior as today (696 assertions), no SymEngine
  dependency, numerical complex DA unchanged.
- Example + README + Doxygen updated. `ref/` untouched.
