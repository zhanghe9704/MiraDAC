# MiraDAC.jl benchmark

The Julia benchmark (plan A.9, T3.x) is measured against the same C++ baseline as the Python
one: `python/bench/bench_cpp.cpp`, driven round by round.

## Build the C++ baseline

```sh
eval "$(scripts/setup_symengine.sh --print-env)"
cmake -S . -B build -G Ninja -DDA_BUILD_CAPI=ON -DDA_BUILD_BENCH=ON -DSymEngine_DIR="$SymEngine_DIR"
cmake --build build --target bench_cpp
```

The binary is `build/python/bench/bench_cpp`. Then `julia julia/dev_setup.jl` points both the
package and this environment (`julia/MiraDAC/bench`) at `build/capi/libmiradac_c.so`.

## Run

```sh
taskset -c 3 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl [--rounds N] [--cpp PATH] [--scoped] [--batch N]
```

Julia cases, per size: `add`, `mul`, `iadd` (`add!(c, c, b)`), `mul_const`, `exp`,
`composition` (`compose(m, n)`), and the in-place API `add!`, `mul!`, `exp!` (compared with the
C++ `add`, `mul`, `exp`). With `--scoped`, the allocating cases also run inside
`dascope`, one scope per batch of `--batch N` operations (default 16), as `scoped_*` (plan A.5b). `cmul` (`c = a * b`) and `cexp`
(`c = exp(a)`) on `CNDA(a, b)`, `CNDA(b, a)`; the symbolic cases come with Stage 5.

## Protocol (as in `python/bench/bench_ops.py`)

- Every case first doubles its repetition count until two samples in a row take >= 10 ms
  (warm-up).
  Each round then gives every case 0.2 s of samples (at least one); a case reports the minimum
  sample over all rounds, in thread CPU time. The data is the same in every round
  (`mt19937(12345)`, coefficients uniform in [-1, 1]). The C++ side uses a pool of 400; the
  Julia side a realistic one (A.9): 10 000 slots at `n3o4` and `n6o6`, 2 000 at `n6o10`.
- Cases, for sizes `n3o4`, `n6o6`, `n6o10` (nvars, order): `add`, `mul`, `iadd`, `mul_const`,
  `exp`, `cmul`, `cexp`, `composition`; plus the symbolic cases `sym_n3o3/*`.
- The driver starts a fresh `bench_cpp --driven` child, which inherits the CPU pinning. For
  each of 15 rounds (default) it writes one line to the child's stdin and waits for the answer
  `done` (the child has run one C++ round), then runs its own round. After the last round it
  closes stdin; the child prints `{"case": ns_per_op, ...}` as JSON on stdout and exits.
- Transparent huge pages are disabled for the driver process (`prctl(PR_SET_THP_DISABLE)`), as
  the C++ side does for itself.
- Run pinned to one core: `taskset -c 3 julia --project=julia/MiraDAC/bench <script>`.

## Gate (plan A.9, revised 2026-09-30)

Each case gets one rule, chosen by its C++ time:

- C++ time < 1 µs: overhead <= 60 ns in place (`iadd`, `add!`, `mul!`, `exp!`), <= 400 ns for
  an operation returning a new object.
- C++ time >= 1 µs: ratio Julia / C++ <= 1.10.

Allocating cases outside a scope: overhead <= 400 ns or ratio <= 1.5 (revised 2026-10-01).
Scoped cases (`scoped_*`, plan A.5b): overhead <= 150 ns below 1 µs, ratio <= 1.10 from 1 µs.
With `--scoped` the script prints `gate A.5b` (the `scoped_*` cases) and `gate A.9` (the
others) apart, then `gate:` for both; it exits 1 if either fails. The unscoped `n6o6`
allocating cases have little A.9 margin and can fail when the machine is loaded.
