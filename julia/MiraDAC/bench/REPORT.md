# MiraDAC.jl benchmark report (plan T3.5)

## Stage 3 result: the gate of plan A.9 PASSES (2 runs of 2, all 27 cases)

Run: `taskset -c 5 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl`
(15 rounds, Julia 1.13.1, one thread; pools of A.9: 10 000 slots at `n3o4` and `n6o6`, 2 000 at
`n6o10`; C++ side `build/python/bench/bench_cpp --driven`, pool 400), each run started when
core 5 was at most 25% busy over a 3 s sample of `/proc/stat` (i7-6700K; core 5 shares a
physical core with core 1, not isolated), 2026-09-30. Library and `bench_cpp`: `build/` at
commit 57615a0 plus the working tree (`capi/`, `julia/`).

In-place operations below 1 µs are at most +17 ns over C++ (limit 60 ns), allocating ones
below 1 µs +145–219 ns (limit 400 ns); every case of 1 µs or more is within ratio 1.10
(`n3o4/composition` 1.086 and 1.078 are the closest).

### Run 1 (started 22:40:37, core 5 busy 24%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        40.6           186.3         145.7   4.589  overhead <= 400: pass
n3o4/mul                       169.1           325.0         155.9   1.922  overhead <= 400: pass
n3o4/iadd                       12.4            17.6           5.2   1.420  overhead <= 60: pass
n3o4/mul_const                  28.1           175.1         147.0   6.231  overhead <= 400: pass
n3o4/add!                       40.6            24.6         -16.0   0.607  overhead <= 60: pass
n3o4/mul!                      169.1           153.8         -15.3   0.910  overhead <= 60: pass
n3o4/exp!                      733.8           751.1          17.3   1.024  overhead <= 60: pass
n3o4/exp                       733.8           904.1         170.3   1.232  overhead <= 400: pass
n3o4/composition              9973.1         10828.2         855.1   1.086  ratio <= 1.10: pass
n6o6/add                       550.2           650.8         100.6   1.183  overhead <= 400: pass
n6o6/mul                     11459.0         11643.2         184.2   1.016  ratio <= 1.10: pass
n6o6/iadd                      189.1           196.6           7.5   1.040  overhead <= 60: pass
n6o6/mul_const                 399.2           608.9         209.7   1.525  overhead <= 400: pass
n6o6/add!                      550.2           264.4        -285.8   0.481  overhead <= 60: pass
n6o6/mul!                    11459.0         11187.9        -271.1   0.976  ratio <= 1.10: pass
n6o6/exp!                    69314.7         69365.8          51.1   1.001  ratio <= 1.10: pass
n6o6/exp                     69314.7         69416.4         101.7   1.001  ratio <= 1.10: pass
n6o6/composition          23950024.0      23680272.0     -269752.0   0.989  ratio <= 1.10: pass
n6o10/add                     7552.6          6771.7        -780.9   0.897  ratio <= 1.10: pass
n6o10/mul                   437576.5        439748.8        2172.2   1.005  ratio <= 1.10: pass
n6o10/iadd                    2290.3          2248.0         -42.3   0.982  ratio <= 1.10: pass
n6o10/mul_const               5367.6          4656.4        -711.2   0.868  ratio <= 1.10: pass
n6o10/add!                    7552.6          3584.1       -3968.5   0.475  ratio <= 1.10: pass
n6o10/mul!                  437576.5        424715.6      -12860.9   0.971  ratio <= 1.10: pass
n6o10/exp!                 4292508.5       4292550.0          41.5   1.000  ratio <= 1.10: pass
n6o10/exp                  4292508.5       4264638.5      -27870.0   0.994  ratio <= 1.10: pass
n6o10/composition      10130588801.0   10001935566.0  -128653235.0   0.987  ratio <= 1.10: pass
gate: PASS
```

### Run 2 (started 22:50:49, core 5 busy 7%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        39.8           184.5         144.7   4.636  overhead <= 400: pass
n3o4/mul                       167.3           322.4         155.1   1.927  overhead <= 400: pass
n3o4/iadd                       12.2            17.5           5.3   1.435  overhead <= 60: pass
n3o4/mul_const                  27.6           172.1         144.5   6.237  overhead <= 400: pass
n3o4/add!                       39.8            23.7         -16.1   0.595  overhead <= 60: pass
n3o4/mul!                      167.3           148.5         -18.8   0.888  overhead <= 60: pass
n3o4/exp!                      726.0           742.5          16.5   1.023  overhead <= 60: pass
n3o4/exp                       726.0           898.3         172.3   1.237  overhead <= 400: pass
n3o4/composition              9820.5         10585.1         764.6   1.078  ratio <= 1.10: pass
n6o6/add                       532.7           732.7         200.0   1.375  overhead <= 400: pass
n6o6/mul                     11399.7         11507.7         108.0   1.009  ratio <= 1.10: pass
n6o6/iadd                      187.8           195.7           7.9   1.042  overhead <= 60: pass
n6o6/mul_const                 386.9           605.7         218.8   1.565  overhead <= 400: pass
n6o6/add!                      532.7           256.7        -276.0   0.482  overhead <= 60: pass
n6o6/mul!                    11399.7         11018.1        -381.6   0.967  ratio <= 1.10: pass
n6o6/exp!                    68813.1         69003.6         190.5   1.003  ratio <= 1.10: pass
n6o6/exp                     68813.1         68790.6         -22.5   1.000  ratio <= 1.10: pass
n6o6/composition          23489984.0      23482858.0       -7126.0   1.000  ratio <= 1.10: pass
n6o10/add                     7385.2          6722.3        -662.9   0.910  ratio <= 1.10: pass
n6o10/mul                   426720.0        426247.1        -472.9   0.999  ratio <= 1.10: pass
n6o10/iadd                    2296.1          2240.5         -55.6   0.976  ratio <= 1.10: pass
n6o10/mul_const               5155.1          4688.6        -466.5   0.910  ratio <= 1.10: pass
n6o10/add!                    7385.2          3553.9       -3831.3   0.481  ratio <= 1.10: pass
n6o10/mul!                  426720.0        423228.5       -3491.5   0.992  ratio <= 1.10: pass
n6o10/exp!                 4257764.0       4264823.2        7059.2   1.002  ratio <= 1.10: pass
n6o10/exp                  4257764.0       4272307.0       14543.0   1.003  ratio <= 1.10: pass
n6o10/composition       9964273232.0    9911352058.0   -52921174.0   0.995  ratio <= 1.10: pass
gate: PASS
```

## Stage 4 rerun (CNDA cases `cmul`, `cexp`): FAILS on `n3o4/cexp` only

Same command and protocol as above, 2026-10-01, `build/` at commit 57615a0 plus the working tree.
Every case passes except `n3o4/cexp` (C++ ≈ 3.0 µs, ratio rule): 1.303 and 1.240. `cmul` passes
at every size; `cexp` at `n6o6`/`n6o10` is within 1.045.

Cause (plan T4.2, *As built*): not the binding (a C loop over `mdac_cnda_exp_into` takes 3000 ns,
C++ 2955 ns), but the C++ `Pool` free list. Complex `exp` uses ~20 temporary slots per call, and
`Pool::alloc` hands slots out in free order. Julia frees results in GC order, so the free list
ends up scattered over the pool of 10 000 and the temporaries land on unrelated pages. The same
C loop after 9 900 vectors were freed in random order: 3188 → 3743 ns. The fix needs a [C++]
change to the pool or a decision on the gate.

### Run 1 (started 01:53:03, core 5 busy 15%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        39.1           189.0         149.9   4.834  overhead <= 400: pass
n3o4/mul                       167.5           322.9         155.4   1.928  overhead <= 400: pass
n3o4/iadd                       12.2            17.5           5.3   1.433  overhead <= 60: pass
n3o4/mul_const                  27.6           174.9         147.3   6.337  overhead <= 400: pass
n3o4/add!                       39.1            23.8         -15.3   0.608  overhead <= 60: pass
n3o4/mul!                      167.5           148.1         -19.4   0.884  overhead <= 60: pass
n3o4/exp!                      729.0           739.8          10.8   1.015  overhead <= 60: pass
n3o4/exp                       729.0           909.4         180.4   1.247  overhead <= 400: pass
n3o4/cmul                      771.2          1080.9         309.7   1.402  overhead <= 400: pass
n3o4/cexp                     2954.6          3851.3         896.7   1.303  ratio <= 1.10: FAIL
n3o4/composition              9884.9         10631.3         746.4   1.076  ratio <= 1.10: pass
n6o6/add                       532.0           738.8         206.8   1.389  overhead <= 400: pass
n6o6/mul                     11405.6         11647.3         241.7   1.021  ratio <= 1.10: pass
n6o6/iadd                      187.8           195.7           7.9   1.042  overhead <= 60: pass
n6o6/mul_const                 387.2           614.0         226.8   1.586  overhead <= 400: pass
n6o6/add!                      532.0           260.8        -271.2   0.490  overhead <= 60: pass
n6o6/mul!                    11405.6         10927.5        -478.1   0.958  ratio <= 1.10: pass
n6o6/exp!                    68392.7         69345.1         952.4   1.014  ratio <= 1.10: pass
n6o6/exp                     68392.7         70046.9        1654.2   1.024  ratio <= 1.10: pass
n6o6/cmul                    47191.7         48786.8        1595.1   1.034  ratio <= 1.10: pass
n6o6/cexp                   177360.1        183386.0        6025.9   1.034  ratio <= 1.10: pass
n6o6/composition          23463020.0      23768607.0      305587.0   1.013  ratio <= 1.10: pass
n6o10/add                     7424.5          6907.3        -517.2   0.930  ratio <= 1.10: pass
n6o10/mul                   431561.9        435087.2        3525.3   1.008  ratio <= 1.10: pass
n6o10/iadd                    2245.9          2261.0          15.1   1.007  ratio <= 1.10: pass
n6o10/mul_const               5133.5          4803.2        -330.3   0.936  ratio <= 1.10: pass
n6o10/add!                    7424.5          3615.1       -3809.4   0.487  ratio <= 1.10: pass
n6o10/mul!                  431561.9        424589.6       -6972.3   0.984  ratio <= 1.10: pass
n6o10/exp!                 4308086.8       4300821.5       -7265.3   0.998  ratio <= 1.10: pass
n6o10/exp                  4308086.8       4277402.0      -30684.8   0.993  ratio <= 1.10: pass
n6o10/cmul                 1752471.0       1774668.6       22197.6   1.013  ratio <= 1.10: pass
n6o10/cexp                 9668906.5       9763848.5       94942.0   1.010  ratio <= 1.10: pass
n6o10/composition      10013428742.0   10067892308.0    54463566.0   1.005  ratio <= 1.10: pass
gate: FAIL
```

### Run 2 (started 02:02:49, core 5 busy 24%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        39.8           188.8         149.0   4.744  overhead <= 400: pass
n3o4/mul                       167.1           323.3         156.2   1.935  overhead <= 400: pass
n3o4/iadd                       12.2            17.8           5.6   1.461  overhead <= 60: pass
n3o4/mul_const                  27.8           175.5         147.7   6.312  overhead <= 400: pass
n3o4/add!                       39.8            24.4         -15.4   0.614  overhead <= 60: pass
n3o4/mul!                      167.1           147.9         -19.2   0.885  overhead <= 60: pass
n3o4/exp!                      728.3           745.0          16.7   1.023  overhead <= 60: pass
n3o4/exp                       728.3           943.3         215.0   1.295  overhead <= 400: pass
n3o4/cmul                      774.2          1130.6         356.4   1.460  overhead <= 400: pass
n3o4/cexp                     3018.0          3740.8         722.8   1.240  ratio <= 1.10: FAIL
n3o4/composition              9991.2         10774.9         783.7   1.078  ratio <= 1.10: pass
n6o6/add                       537.3           736.7         199.4   1.371  overhead <= 400: pass
n6o6/mul                     11447.6         11825.7         378.1   1.033  ratio <= 1.10: pass
n6o6/iadd                      188.1           195.8           7.7   1.041  overhead <= 60: pass
n6o6/mul_const                 387.2           616.8         229.6   1.593  overhead <= 400: pass
n6o6/add!                      537.3           257.6        -279.7   0.480  overhead <= 60: pass
n6o6/mul!                    11447.6         10980.9        -466.7   0.959  ratio <= 1.10: pass
n6o6/exp!                    68867.9         69481.2         613.3   1.009  ratio <= 1.10: pass
n6o6/exp                     68867.9         70963.6        2095.7   1.030  ratio <= 1.10: pass
n6o6/cmul                    47261.4         49407.2        2145.8   1.045  ratio <= 1.10: pass
n6o6/cexp                   178416.6        185661.7        7245.1   1.041  ratio <= 1.10: pass
n6o6/composition          23669656.0      23785319.0      115663.0   1.005  ratio <= 1.10: pass
n6o10/add                     7410.4          6732.1        -678.3   0.908  ratio <= 1.10: pass
n6o10/mul                   432115.4        435518.3        3402.9   1.008  ratio <= 1.10: pass
n6o10/iadd                    2229.1          2328.4          99.3   1.045  ratio <= 1.10: pass
n6o10/mul_const               5193.4          4686.7        -506.7   0.902  ratio <= 1.10: pass
n6o10/add!                    7410.4          3650.2       -3760.2   0.493  ratio <= 1.10: pass
n6o10/mul!                  432115.4        424908.6       -7206.8   0.983  ratio <= 1.10: pass
n6o10/exp!                 4325521.0       4296303.5      -29217.5   0.993  ratio <= 1.10: pass
n6o10/exp                  4325521.0       4287361.5      -38159.5   0.991  ratio <= 1.10: pass
n6o10/cmul                 1772787.1       1770555.5       -2231.6   0.999  ratio <= 1.10: pass
n6o10/cexp                 9702159.0       9744949.0       42790.0   1.004  ratio <= 1.10: pass
n6o10/composition      10080597661.0   10104250488.0    23652827.0   1.002  ratio <= 1.10: pass
gate: FAIL
```

## History: why the gate failed before

Earlier runs of this gate failed for three reasons, each a systematic cost in the C++ library,
each fixed there:

1. **Engine temporaries leaked when the pool ran out** (commit a6be310). `exp` and `compose`
   loops threw `PoolExhaustedError` even with the A.5 retry, because a kernel that threw
   mid-function did not free its `pool.alloc()` temporaries (one failed `da::exp` leaked 1
   slot, `da_composition` 5). Now freed when a kernel throws (Catch2 `[numeric_leak]`);
   `test/test_memory.jl` checks the retry in `exp`/`compose` loops with a pool of 64.
2. **Product-index rows scattered in a fragmented heap** (commit 3dec5a2). The cases dominated
   by the multiplication kernel (`mul`, `exp`, `composition` at `n6o6`/`n6o10`) were 11–14%
   over C++, `n3o4/composition` 1.12–1.22. Not the binding: a C loop over
   `mdac_ndalist_compose` + `mdac_ndalist_free` took the same time as C++ `da_composition`
   (`n3o4` 10 720 vs 10 671 ns, `n6o6` 25.75 vs 25.72 ms). The same C code ran 10–14% slower in
   any process with a fragmented heap (a Julia process, or a C program that first allocates
   and frees every other of 100 000 random-size blocks: 27.4 → 31.9 ms at `n6o6`), because
   `Layout::init_prod_index` allocated each row of `prdidx` separately and the kernel, which
   walks them, paid in cache and TLB misses. A fresh C++ process (`bench_cpp`) got the rows
   nearly contiguous. `prdidx` is now one block owned by `Layout` (Catch2 "prdidx rows are
   contiguous").
3. **Cold free-time zeroing of 64 KB slots** (commit 57615a0). With 3dec5a2 in, the gate still
   failed on `n6o10/mul_const` (`a * 2.0`, C++ ≈ 6.6 µs) in 4 of 5 runs (ratios 1.379, 1.200,
   1.141, 1.122; 0.720 once). `Pool::free` zeroed the whole slot (and `alloc` zeroed it again).
   C++ frees `c`'s old slot right after writing it, so the zeroing hit cache; Julia frees a
   result only after a garbage collection (A.5), when its 64 KB slot has left the cache (the
   pool cycles 128 MB between collections), so each allocating operation at `n6o10` paid
   ~2.5 µs of cold writes (`a * 2.0` from Julia 7.6–8.2 µs; `mul!(c, a, 2.0)` 2.8 µs). Slots of
   doubles are now zeroed in `assign` (where they are about to be written), not in `free`;
   `n6o10/mul_const` is at 0.868 and 0.910 above, and C++ itself gains (≈ 6.6 → 5.2–5.4 µs).

The gate itself was revised once (A.9, approved by the user 2026-09-30): allocating operations
cannot reach the Python limit of 150 ns, since every new Julia object pays a finalizer, the
deferred free and its share of the GC that the pool-exhaustion retry triggers; the in-place API
is the performance path with the 60 ns limit, and cases of 1 µs or more have a ratio rule.

Remedies 1 and 2 of T3.5 found nothing in the binding: `@code_warntype` is clean on every
benchmarked operation, and the empty free-queue check costs ~7 ns.
