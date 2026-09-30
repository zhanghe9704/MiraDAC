# MiraDAC Python binding — performance report

The cost of driving MiraDAC from Python, per operation, against the same operation in C++
(plan T2.5, T8.4). Run of 2026-09-30, at the end of Stage 8.

## Method

`python/bench/bench_ops.py` times each case in Python and, interleaved round by round, in a
`bench_cpp --driven` child process (the C++ baseline, `python/bench/bench_cpp.cpp`):

- Dense random DA vectors (every coefficient set); symbolic cases give coefficient `i` its own
  symbol.
- Warm-up: the repetition count doubles until one sample takes at least 10 ms.
- 15 rounds; in each round every case gets 0.2 s of samples. A case reports its **minimum**
  sample, in **thread CPU time** (`CLOCK_THREAD_CPUTIME_ID` / `time.thread_time_ns`).
- Both processes pinned to one core (`taskset -c 5`), transparent huge pages off (`prctl`).

Gate: C++ time below 1 µs, overhead (Python − C++) at most 150 ns for an operation that returns a
new object and at most 100 ns for an in-place one (`iadd`); C++ time 1 µs or more, ratio
(Python / C++) at most 1.10; `sym_n3o3/exp`, ratio at most 1.05.

Reproduce:

```bash
eval "$(scripts/setup_symengine.sh --print-env)"
cmake -S . -B build-bench -G Ninja -DCMAKE_BUILD_TYPE=Release -DDA_BUILD_BENCH=ON \
      -DDA_BUILD_TESTS=OFF -DDA_BUILD_EXAMPLES=OFF
cmake --build build-bench
taskset -c 5 .venv/bin/python python/bench/bench_ops.py
```

## Machine and toolchain

| Item | Value |
|---|---|
| CPU | Intel Core i7-6700K, 4 cores / 8 threads, 4.0 GHz (4.2 GHz max), `powersave` governor |
| Memory | 48 GB |
| OS | Debian 13.7, Linux 6.12.107+deb13-amd64 |
| Compiler | GCC 14.2.0 (Debian 14.2.0-19), `-O3 -DNDEBUG` (CMake `Release`) for both the C++ baseline and the extension |
| Python | CPython 3.13.5 |
| Binding | nanobind 3.1.0 (`NB_STATIC`, `NOMINSIZE`, no LTO), scikit-build-core 1.1.0 |
| SymEngine | commit `153b7e98` (0.14.0 + 14 commits), GMP integers, the pinned build |
| numpy | 2.5.3 |
| Load | load average 1.31 at start, 1.52 at end (desktop browser on other cores) |

## Results

Case names are `n<nvars>o<order>/<op>`: `add` (`c = a + b`), `mul` (`c = a * b`), `iadd`
(`c += b`), `mul_const` (`c = a * 2.0`), `exp` (`c = exp(a)`), `cmul`/`cexp` (the same on
CNDA), `composition` (an nvars-map composed with an nvars-map, output list given). `sym_n3o3`
are SDA at 3 variables, order 3.

| case | C++ ns | Python ns | overhead ns | ratio | gate |
|---|---:|---:|---:|---:|---|
| n3o4/add | 40.9 | 125.1 | 84.2 | 3.060 | overhead ≤ 150: pass |
| n3o4/mul | 171.7 | 269.4 | 97.7 | 1.569 | overhead ≤ 150: pass |
| n3o4/iadd | 12.2 | 60.9 | 48.7 | 4.991 | overhead ≤ 100: pass |
| n3o4/mul_const | 31.4 | 126.7 | 95.3 | 4.035 | overhead ≤ 150: pass |
| n3o4/exp | 744.2 | 832.5 | 88.3 | 1.119 | overhead ≤ 150: pass |
| n3o4/cmul | 814.6 | 932.6 | 118.0 | 1.145 | overhead ≤ 150: pass |
| n3o4/cexp | 3 142.4 | 3 419.4 | 277.0 | 1.088 | ratio ≤ 1.10: pass |
| n3o4/composition | 9 659.4 | 9 717.0 | 57.6 | 1.006 | ratio ≤ 1.10: pass |
| n6o6/add | 594.5 | 679.3 | 84.8 | 1.143 | overhead ≤ 150: pass |
| n6o6/mul | 11 510.8 | 11 748.1 | 237.3 | 1.021 | ratio ≤ 1.10: pass |
| n6o6/iadd | 188.3 | 239.7 | 51.4 | 1.273 | overhead ≤ 100: pass |
| n6o6/mul_const | 449.5 | 530.1 | 80.6 | 1.179 | overhead ≤ 150: pass |
| n6o6/exp | 69 457.7 | 69 705.9 | 248.2 | 1.004 | ratio ≤ 1.10: pass |
| n6o6/cmul | 48 491.1 | 48 527.5 | 36.4 | 1.001 | ratio ≤ 1.10: pass |
| n6o6/cexp | 183 504.0 | 185 905.9 | 2 401.9 | 1.013 | ratio ≤ 1.10: pass |
| n6o6/composition | 23 560 754 | 23 892 186 | 331 432 | 1.014 | ratio ≤ 1.10: pass |
| n6o10/add | 8 976.3 | 7 601.9 | −1 374.4 | 0.847 | ratio ≤ 1.10: pass |
| n6o10/mul | 428 585.6 | 426 311.6 | −2 274.0 | 0.995 | ratio ≤ 1.10: pass |
| n6o10/iadd | 2 237.6 | 2 269.2 | 31.6 | 1.014 | ratio ≤ 1.10: pass |
| n6o10/mul_const | 6 796.5 | 5 383.9 | −1 412.6 | 0.792 | ratio ≤ 1.10: pass |
| n6o10/exp | 4 258 436.5 | 4 266 826.0 | 8 389.5 | 1.002 | ratio ≤ 1.10: pass |
| n6o10/cmul | 1 750 679.6 | 1 742 919.8 | −7 759.9 | 0.996 | ratio ≤ 1.10: pass |
| n6o10/cexp | 9 700 277.0 | 9 663 996.5 | −36 280.5 | 0.996 | ratio ≤ 1.10: pass |
| n6o10/composition | 9 889 671 392 | 9 941 714 806 | 52 043 414 | 1.005 | ratio ≤ 1.10: pass |
| sym_n3o3/mul | 34 277.8 | 34 148.9 | −128.9 | 0.996 | ratio ≤ 1.10: pass |
| sym_n3o3/exp | 70 449.9 | 70 457.8 | 7.9 | 1.000 | ratio ≤ 1.05: pass |

**Gate: PASS** on all 26 cases.

## Reading the numbers

- A bound call costs a fixed 50–120 ns (loop iteration, nanobind dispatch, the env guard,
  allocating the result object): 49–51 ns in place, 81–118 ns for a new object. That dominates
  only below about 1 µs of C++ work, i.e. at low order.
- From order 6 up, and for every symbolic case, Python is within the noise of C++ (ratios
  0.99–1.02, apart from the cases below).
- `n6o10/add` and `n6o10/mul_const` are 15–21% faster from Python. The C++ `c = a + b` returns by
  value and moves into `c`, and the move constructor takes and frees an extra pool slot, zeroing
  it at full length (8008 doubles at (6,10)). The binding constructs the result in the Python
  object directly and skips that (plan T2.6 step 3).
- `n3o4/cexp` (1.088) is the case closest to its limit; under heavier machine load it has gone
  above 1.10 in earlier runs (plan T5.6, T7.4).
- LTO on the extension and the library was tried and dropped: 0.993× median over the cases.
