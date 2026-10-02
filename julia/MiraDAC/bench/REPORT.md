# MiraDAC.jl benchmark report (plan T3.5)

## Stage 3 result: the gate of plan A.9 PASSED (2 runs of 2, all 27 cases; first revision of the gate, before efc0f70)

## A.5b result (branch `julia-alloc-scope`): all 21 `dascope` cases meet the strict rule (2 runs of 2)

Run: `taskset -c 3 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl --scoped`
(15 rounds, Julia 1.13.1, one thread; pools of A.9; `scoped_*` cases run one `dascope` per batch
of 16 operations). Each run started when core 3 was at most 25% busy over a 3 s sample of
`/proc/stat`, 2026-10-02. Library and `bench_cpp`: worktree `build/` (Release) at commit 43f61f5
plus the working tree. Core 3 (not 5, where another benchmark ran), so absolute numbers differ a
little from the sections below.

There are 21 scoped cases: `scoped_add`, `scoped_mul`, `scoped_mul_const`, `scoped_exp`,
`scoped_cmul`, `scoped_cexp` and `scoped_composition` at each of the three sizes. The strict rule
for them (plan A.5b) is overhead <= 150 ns when C++ < 1 µs and ratio <= 1.10 when C++ >= 1 µs.
**All 21 pass in both runs (`gate A.5b (scoped_* cases): PASS`).** The closest are
`n3o4/scoped_exp` (+106 and +109 ns), `n3o4/scoped_cmul` (+97 and +93 ns) and `n3o4/scoped_cexp`
(1.075 and 1.090). The script prints the two gates on separate lines. Run 3 exits 0
(`gate: PASS`). Run 4 exits 1 because one **unscoped** case fails A.9:
`n6o6/mul_const` at +424 ns against a limit of 400 (see below). The same cases in and out of a
scope:

| case | unscoped overhead or ratio (run 3 / run 4) | scoped (run 3 / run 4) |
|---|---|---|
| n3o4/add | +49 / +166 ns | +64 / +59 ns |
| n3o4/mul_const | +170 / +169 ns | +67 / +69 ns |
| n3o4/cmul | +235 / +238 ns | +97 / +93 ns |
| n6o6/add | +294 / +337 ns | +77 / +91 ns |
| n6o6/mul_const | +377 / **+424** ns | +58 / +81 ns |
| n6o10/add | 1.257 / 1.291 | 0.954 / 0.980 |
| n6o10/mul_const | 1.460 / 1.423 | 0.945 / 0.965 |

**Unscoped `n6o6` cases and A.9.** These cases sit at the A.9 limit with or without this change.
An independent verifier saw `n6o6/add` fail twice (+410 and +439 ns, load average 6.0 and 3.8),
and run 4 here saw `n6o6/mul_const` fail at a load average of 1.8. An A/B against HEAD, by the
verifier (unscoped, n6o6, pool 10 000, core 3, three alternating pairs, best of 60 × 40k), shows
no slowdown: `a + b` took 923.8 / 898.4 / 888.5 ns at HEAD and 906.7 / 891.0 / 901.7 ns in the
worktree; `a * 2.0` took 789.2 / 786.6 / 789.4 ns at HEAD and 818.4 / 793.8 / 795.3 ns in the
worktree. The overhead does not come from the Julia layer. It comes from the cold slot. Raw C
calls with no Julia objects show it (`mdac_nda_add`, `mdac_nda_mul_d` through `ccall`, n6o6,
pool 10 000, core 3, load 2.5, 6 repetitions):

| | each result freed at once (warm slot, as C++ and `dascope`) | 9 000 results kept, then freed (cold slots, as unscoped Julia between GCs) | unscoped Julia `a + b` / `a * 2.0` |
|---|---|---|---|
| add | 396–446 ns | 706–750 ns | 795–849 ns |
| mul_const | 160–174 ns | 410–435 ns | 732–755 ns |

So about 250–300 ns of the unscoped overhead is the cold 7.4 KB slot that every result gets
when results are freed only after a GC. That cost depends on memory traffic from other processes,
not only on core 3. The finalizer and deferred free add about 100 ns. Plan A.5b item 7 keeps the
unscoped behaviour unchanged, so this cost stays outside a scope. A scope removes both parts:
`n6o6/scoped_add` is +77 / +91 ns.

### Run 3 (started 01:36:12, core 3 busy 17%, load average 4.29)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        34.0            82.8          48.8   2.435  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       165.5           339.0         173.5   2.048  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.1           5.5   1.435  overhead <= 60: pass
n3o4/mul_const                  21.5           191.1         169.6   8.888  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       34.0            24.9          -9.1   0.733  overhead <= 60: pass
n3o4/mul!                      165.5           153.1         -12.4   0.925  overhead <= 60: pass
n3o4/exp!                      696.9           710.0          13.1   1.019  overhead <= 60: pass
n3o4/exp                       696.9           867.9         171.0   1.245  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      739.1           974.1         235.0   1.318  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2753.9          2980.2         226.3   1.082  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition             10042.4         10324.0         281.6   1.028  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       411.2           704.9         293.7   1.714  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11455.1         11919.7         464.6   1.041  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      193.1           201.8           8.7   1.045  overhead <= 60: pass
n6o6/mul_const                 304.1           681.5         377.4   2.241  overhead <= 400 or ratio <= 1.5: pass
n6o6/add!                      411.2           267.7        -143.5   0.651  overhead <= 60: pass
n6o6/mul!                    11455.1         11301.8        -153.3   0.987  ratio <= 1.10: pass
n6o6/exp!                    69833.8         69673.0        -160.8   0.998  ratio <= 1.10: pass
n6o6/exp                     69833.8         70696.8         863.0   1.012  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    47600.6         47830.9         230.3   1.005  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   175447.2        177547.4        2100.2   1.012  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24280524.0      24429802.0      149278.0   1.006  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5839.6          7343.0        1503.4   1.257  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   438742.4        444690.9        5948.5   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2279.9          2269.8         -10.1   0.996  ratio <= 1.10: pass
n6o10/mul_const               3545.5          5176.8        1631.3   1.460  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5839.6          3673.6       -2166.0   0.629  ratio <= 1.10: pass
n6o10/mul!                  438742.4        440044.2        1301.8   1.003  ratio <= 1.10: pass
n6o10/exp!                 4437890.8       4430725.0       -7165.8   0.998  ratio <= 1.10: pass
n6o10/exp                  4437890.8       4469936.5       32045.7   1.007  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1773474.8       1798933.4       25458.6   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9711598.0       9663136.0      -48462.0   0.995  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      10810513801.0   11003655350.0   193141549.0   1.018  overhead <= 400 or ratio <= 1.5: pass
n3o4/scoped_add                 34.0            98.4          64.4   2.894  overhead <= 150: pass
n3o4/scoped_mul                165.5           230.0          64.5   1.390  overhead <= 150: pass
n3o4/scoped_mul_const            21.5            88.8          67.3   4.132  overhead <= 150: pass
n3o4/scoped_exp                696.9           802.6         105.7   1.152  overhead <= 150: pass
n3o4/scoped_cmul               739.1           836.3          97.2   1.132  overhead <= 150: pass
n3o4/scoped_cexp              2753.9          2961.8         207.9   1.075  ratio <= 1.10: pass
n3o4/scoped_composition         10042.4         10422.1         379.7   1.038  ratio <= 1.10: pass
n6o6/scoped_add                411.2           487.8          76.6   1.186  overhead <= 150: pass
n6o6/scoped_mul              11455.1         11731.7         276.6   1.024  ratio <= 1.10: pass
n6o6/scoped_mul_const           304.1           362.3          58.2   1.192  overhead <= 150: pass
n6o6/scoped_exp              69833.8         70358.8         525.0   1.008  ratio <= 1.10: pass
n6o6/scoped_cmul             47600.6         47625.2          24.6   1.001  ratio <= 1.10: pass
n6o6/scoped_cexp            175447.2        176709.5        1262.3   1.007  ratio <= 1.10: pass
n6o6/scoped_composition      24280524.0      24363830.0       83306.0   1.003  ratio <= 1.10: pass
n6o10/scoped_add              5839.6          5572.3        -267.3   0.954  ratio <= 1.10: pass
n6o10/scoped_mul            438742.4        443610.5        4868.1   1.011  ratio <= 1.10: pass
n6o10/scoped_mul_const          3545.5          3351.3        -194.2   0.945  ratio <= 1.10: pass
n6o10/scoped_exp           4437890.8       4447704.2        9813.5   1.002  ratio <= 1.10: pass
n6o10/scoped_cmul          1773474.8       1789462.4       15987.6   1.009  ratio <= 1.10: pass
n6o10/scoped_cexp          9711598.0       9649526.0      -62072.0   0.994  ratio <= 1.10: pass
n6o10/scoped_composition   10810513801.0   11093600145.0   283086344.0   1.026  ratio <= 1.10: pass
gate A.5b (scoped_* cases): PASS
gate A.9 (unscoped cases): PASS
gate: PASS
exit=0
```

### Run 4 (started 01:52:31, core 3 busy 20%, load average 1.79)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        34.0           199.7         165.7   5.872  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       162.2           294.9         132.7   1.818  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.5            18.0           5.5   1.437  overhead <= 60: pass
n3o4/mul_const                  21.2           190.6         169.4   8.993  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       34.0            24.4          -9.6   0.719  overhead <= 60: pass
n3o4/mul!                      162.2           152.4          -9.8   0.940  overhead <= 60: pass
n3o4/exp!                      699.6           716.8          17.2   1.025  overhead <= 60: pass
n3o4/exp                       699.6           876.1         176.5   1.252  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      732.6           970.1         237.5   1.324  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2725.3          2927.8         202.5   1.074  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition              9955.6         10239.4         283.8   1.029  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       406.2           743.6         337.4   1.831  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11505.3         11968.5         463.2   1.040  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      193.1           202.0           8.9   1.046  overhead <= 60: pass
n6o6/mul_const                 303.9           727.6         423.7   2.394  overhead <= 400 or ratio <= 1.5: FAIL
n6o6/add!                      406.2           267.8        -138.4   0.659  overhead <= 60: pass
n6o6/mul!                    11505.3         11363.2        -142.1   0.988  ratio <= 1.10: pass
n6o6/exp!                    69845.2         69903.8          58.6   1.001  ratio <= 1.10: pass
n6o6/exp                     69845.2         70806.1         960.9   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    47838.7         48153.9         315.2   1.007  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   175756.4        180526.7        4770.3   1.027  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24262892.0      24331791.0       68899.0   1.003  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5888.0          7599.5        1711.5   1.291  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   440485.9        444231.1        3745.2   1.009  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2280.1          2326.2          46.1   1.020  ratio <= 1.10: pass
n6o10/mul_const               3719.6          5294.6        1575.0   1.423  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5888.0          3740.4       -2147.6   0.635  ratio <= 1.10: pass
n6o10/mul!                  440485.9        436230.4       -4255.5   0.990  ratio <= 1.10: pass
n6o10/exp!                 4424821.7       4381513.8      -43308.0   0.990  ratio <= 1.10: pass
n6o10/exp                  4424821.7       4420244.0       -4577.7   0.999  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1780441.8       1787869.4        7427.6   1.004  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9633570.0       9612635.0      -20935.0   0.998  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      10977853102.0   10942671186.0   -35181916.0   0.997  overhead <= 400 or ratio <= 1.5: pass
n3o4/scoped_add                 34.0            93.3          59.3   2.745  overhead <= 150: pass
n3o4/scoped_mul                162.2           232.4          70.2   1.433  overhead <= 150: pass
n3o4/scoped_mul_const            21.2            90.4          69.2   4.263  overhead <= 150: pass
n3o4/scoped_exp                699.6           808.6         109.0   1.156  overhead <= 150: pass
n3o4/scoped_cmul               732.6           825.8          93.2   1.127  overhead <= 150: pass
n3o4/scoped_cexp              2725.3          2971.2         245.9   1.090  ratio <= 1.10: pass
n3o4/scoped_composition          9955.6         10383.5         427.9   1.043  ratio <= 1.10: pass
n6o6/scoped_add                406.2           497.1          90.9   1.224  overhead <= 150: pass
n6o6/scoped_mul              11505.3         11691.4         186.1   1.016  ratio <= 1.10: pass
n6o6/scoped_mul_const           303.9           385.1          81.2   1.267  overhead <= 150: pass
n6o6/scoped_exp              69845.2         70496.4         651.2   1.009  ratio <= 1.10: pass
n6o6/scoped_cmul             47838.7         47525.2        -313.5   0.993  ratio <= 1.10: pass
n6o6/scoped_cexp            175756.4        177874.4        2118.0   1.012  ratio <= 1.10: pass
n6o6/scoped_composition      24262892.0      24269701.0        6809.0   1.000  ratio <= 1.10: pass
n6o10/scoped_add              5888.0          5772.4        -115.6   0.980  ratio <= 1.10: pass
n6o10/scoped_mul            440485.9        442230.3        1744.4   1.004  ratio <= 1.10: pass
n6o10/scoped_mul_const          3719.6          3591.1        -128.5   0.965  ratio <= 1.10: pass
n6o10/scoped_exp           4424821.7       4429662.5        4840.8   1.001  ratio <= 1.10: pass
n6o10/scoped_cmul          1780441.8       1784595.2        4153.4   1.002  ratio <= 1.10: pass
n6o10/scoped_cexp          9633570.0       9728517.0       94947.0   1.010  ratio <= 1.10: pass
n6o10/scoped_composition   10977853102.0   10922372024.0   -55481078.0   0.995  ratio <= 1.10: pass
gate A.5b (scoped_* cases): PASS
gate A.9 (unscoped cases): FAIL
gate: FAIL
exit=1
```


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

## Stage 7 rerun (multi-env API): blocking gate PASSED (allocating cases informational, A.9 as of 2026-10-01 later)

Same command and protocol as the Stage 4 result below, 2026-10-01, `build/` at commit 43f61f5 plus
the working tree (C API rebuilt with the T7.1 functions). One run, started 18:08:54 with core 5 4%
busy and a 1-minute load average of 1.45; another benchmark (a different checkout) ran pinned to
core 3 during the run. Exit status 0: blocking PASS (in-place cases below 1 µs at most +23.2 ns,
of 1 µs or more at most ratio 1.024); every allocating case passes too (`n3o4/composition`
+390.8 ns and `n6o6/mul_const` +376.5 ns closest to 400 ns; `n6o6/mul` +413.5 ns passes by its
ratio 1.036).

### Run (core 5 busy 4%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        34.1            77.3          43.2   2.266  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       160.5           330.2         169.7   2.057  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.3            18.1           5.8   1.471  overhead <= 60: pass
n3o4/mul_const                  21.7           182.1         160.4   8.391  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       34.1            25.3          -8.8   0.741  overhead <= 60: pass
n3o4/mul!                      160.5           149.8         -10.7   0.933  overhead <= 60: pass
n3o4/exp!                      682.4           705.6          23.2   1.034  overhead <= 60: pass
n3o4/exp                       682.4           870.7         188.3   1.276  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      728.4           950.6         222.2   1.305  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2726.3          2955.1         228.8   1.084  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition              9832.1         10222.9         390.8   1.040  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       394.8           692.9         298.1   1.755  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11344.2         11757.7         413.5   1.036  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      189.0           195.9           6.9   1.036  overhead <= 60: pass
n6o6/mul_const                 303.0           679.5         376.5   2.243  overhead <= 400 or ratio <= 1.5: pass
n6o6/add!                      394.8           262.1        -132.7   0.664  overhead <= 60: pass
n6o6/mul!                    11344.2         10999.1        -345.1   0.970  ratio <= 1.10: pass
n6o6/exp!                    68047.5         67961.3         -86.2   0.999  ratio <= 1.10: pass
n6o6/exp                     68047.5         69236.5        1189.0   1.017  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    46415.9         47142.0         726.1   1.016  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   172767.4        178683.9        5916.5   1.034  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          23839917.0      23876974.0       37057.0   1.002  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5725.8          7054.8        1329.0   1.232  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   442304.3        442449.6         145.3   1.000  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2245.7          2298.5          52.8   1.024  ratio <= 1.10: pass
n6o10/mul_const               3577.3          4873.2        1295.9   1.362  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5725.8          3779.3       -1946.5   0.660  ratio <= 1.10: pass
n6o10/mul!                  442304.3        436542.9       -5761.4   0.987  ratio <= 1.10: pass
n6o10/exp!                 4357907.7       4354708.2       -3199.5   0.999  ratio <= 1.10: pass
n6o10/exp                  4357907.7       4335909.0      -21998.7   0.995  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1797798.5       1774046.4      -23752.1   0.987  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9677705.0       9655920.5      -21784.5   0.998  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      10383334977.0   10339470430.0   -43864547.0   0.996  overhead <= 400 or ratio <= 1.5: pass
gate (blocking): PASS
allocating (informational): pass
```

### Run 2 (2026-10-01, after the T7.3 test additions; machine loaded)

Same command, started 20:43:48 with a 1-minute load average of 4.87 (another checkout's
`bench_ops.jl --scoped` and `bench_cpp --driven` pinned to core 3). Exit status 0: blocking PASS
(in-place cases below 1 µs at most +15.9 ns, of 1 µs or more at most ratio 1.056). Allocating
(informational): `n6o6/mul_const` +426.8 ns over 400 ns, the case already at the limit in the
Stage 4 rerun after efc0f70; `n3o4/composition` +415.2 ns passes by its ratio 1.042.

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        35.1           202.4         167.3   5.767  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       167.7           337.3         169.6   2.012  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.2           5.6   1.446  overhead <= 60: pass
n3o4/mul_const                  21.4           184.6         163.2   8.626  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       35.1            25.2          -9.9   0.718  overhead <= 60: pass
n3o4/mul!                      167.7           153.4         -14.3   0.915  overhead <= 60: pass
n3o4/exp!                      704.8           720.7          15.9   1.023  overhead <= 60: pass
n3o4/exp                       704.8           887.2         182.4   1.259  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      742.7           990.5         247.8   1.334  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2912.8          3119.3         206.5   1.071  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition              9954.2         10369.4         415.2   1.042  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       412.7           752.6         339.9   1.824  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11708.1         12031.5         323.4   1.028  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      192.9           202.4           9.5   1.049  overhead <= 60: pass
n6o6/mul_const                 306.1           732.9         426.8   2.394  overhead <= 400 or ratio <= 1.5: FAIL
n6o6/add!                      412.7           268.6        -144.1   0.651  overhead <= 60: pass
n6o6/mul!                    11708.1         11562.0        -146.1   0.988  ratio <= 1.10: pass
n6o6/exp!                    69554.6         70088.2         533.6   1.008  ratio <= 1.10: pass
n6o6/exp                     69554.6         71918.3        2363.7   1.034  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    47793.1         48964.4        1171.3   1.025  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   178643.3        181020.1        2376.8   1.013  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24553202.0      24948577.0      395375.0   1.016  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5702.4          7206.4        1504.0   1.264  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   450895.3        446802.3       -4093.0   0.991  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2237.7          2288.2          50.5   1.023  ratio <= 1.10: pass
n6o10/mul_const               3535.2          5255.0        1719.8   1.486  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5702.4          3724.7       -1977.7   0.653  ratio <= 1.10: pass
n6o10/mul!                  450895.3        447900.9       -2994.4   0.993  ratio <= 1.10: pass
n6o10/exp!                 4408862.0       4656036.0      247174.0   1.056  ratio <= 1.10: pass
n6o10/exp                  4408862.0       4517933.8      109071.8   1.025  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1791932.8       1780008.5      -11924.3   0.993  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9796234.0       9744784.0      -51450.0   0.995  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      12189446324.0   12153511011.0   -35935313.0   0.997  overhead <= 400 or ratio <= 1.5: pass
gate (blocking): PASS
allocating (informational): fail n6o6/mul_const
```

## Stage 4 result: blocking gate PASSED (allocating cases informational, A.9 as of 2026-10-01 later)

Same command and protocol as above, 2026-10-01, `build/` at commit efc0f70 plus the working tree
(rebuilt, no work to do). `bench_ops.jl` now prints two verdicts: `gate (blocking)` over the
in-place cases under their strict rules (it alone sets the exit status), and `allocating
(informational)` with the allocating cases under "overhead <= 400 ns or ratio <= 1.5", which do
not block until `dascope` lands (A.9). One run, exit status 0: blocking PASS; the allocating
cases all pass too, `n6o6/mul_const` again closest to its limit (+380.6 ns). In-place cases below
1 µs are at most +12 ns over C++; in-place cases of 1 µs or more at most ratio 1.018.

### Run (core 5 busy 17%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        34.6           202.7         168.1   5.860  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       165.7           337.0         171.3   2.034  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.2           5.6   1.445  overhead <= 60: pass
n3o4/mul_const                  21.2           186.5         165.3   8.795  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       34.6            24.9          -9.7   0.719  overhead <= 60: pass
n3o4/mul!                      165.7           152.3         -13.4   0.919  overhead <= 60: pass
n3o4/exp!                      703.8           715.9          12.1   1.017  overhead <= 60: pass
n3o4/exp                       703.8           875.0         171.2   1.243  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      737.2           991.0         253.8   1.344  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2753.5          2981.0         227.5   1.083  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition             10078.1         10370.0         291.9   1.029  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       409.9           714.7         304.8   1.744  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11656.0         11993.1         337.1   1.029  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      193.7           202.1           8.4   1.043  overhead <= 60: pass
n6o6/mul_const                 305.5           686.1         380.6   2.246  overhead <= 400 or ratio <= 1.5: pass
n6o6/add!                      409.9           269.9        -140.0   0.658  overhead <= 60: pass
n6o6/mul!                    11656.0         11307.0        -349.0   0.970  ratio <= 1.10: pass
n6o6/exp!                    70842.8         70069.3        -773.5   0.989  ratio <= 1.10: pass
n6o6/exp                     70842.8         71208.9         366.1   1.005  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    48185.0         48538.2         353.2   1.007  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   177610.9        179182.0        1571.1   1.009  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24825923.0      24606175.0     -219748.0   0.991  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5659.9          7360.3        1700.4   1.300  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   440862.2        449699.4        8837.2   1.020  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2298.9          2289.0          -9.9   0.996  ratio <= 1.10: pass
n6o10/mul_const               3547.6          5179.5        1631.9   1.460  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5659.9          3578.1       -2081.8   0.632  ratio <= 1.10: pass
n6o10/mul!                  440862.2        445778.5        4916.3   1.011  ratio <= 1.10: pass
n6o10/exp!                 4414785.0       4493529.0       78744.0   1.018  ratio <= 1.10: pass
n6o10/exp                  4414785.0       4549097.8      134312.8   1.030  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1784727.2       1837396.9       52669.7   1.030  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9672664.0       9784561.0      111897.0   1.012  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      11788502815.0   11588281780.0  -200221035.0   0.983  overhead <= 400 or ratio <= 1.5: pass
gate (blocking): PASS
allocating (informational): pass
```

## Stage 4 rerun after commit efc0f70 and the 2026-10-01 gate: `n6o6/mul_const` at the 400 ns limit (2 runs of 4 pass)

Same command and protocol as above, 2026-10-01, `build/` rebuilt at commit efc0f70 (pool free list is a
stack) plus the working tree; gate of A.9 as revised 2026-10-01 (allocating operations: overhead
<= 400 ns or ratio <= 1.5 at every size). Four runs: 1 and 4 pass, 2 and 3 fail, each on the single case
`n6o6/mul_const` (`a * 2.0`, C++ ≈ 300 ns, so a ratio of 1.5 means 150 ns and the 400 ns overhead is the
binding limit): +319, +423, +404, +380 ns. Every other case passes in every run; `n3o4/cexp`, which
failed before, is at ratio 1.056, 1.083, 1.075, 1.095 (overhead 157–252 ns); `cmul` and `cexp` pass at
every size. In-place cases below 1 µs are at most +22 ns over C++.

Why `n6o6/mul_const` moved to the limit: with the stack free list C++ reuses the slot it has just freed,
so its `a * 2.0` got faster (387–399 ns in the Stage 3 tables, 296–307 ns now), while the Julia time
went from 605–617 ns to 627–727 ns. A Julia result is freed only after a GC (A.5), so the slot it gets
is not warm (7.4 KB at `n6o6`) — the cost that A.9's 2026-10-01 revision describes for slots of 1 µs
and more, here on a case below 1 µs. The fix is the separate `dascope` work (A.9) or a gate decision;
both are outside Stage 4.

### Run 1 (started 07:57:47, core 5 busy 23%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        36.2           195.5         159.3   5.400  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       164.6           338.4         173.8   2.056  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.4           5.8   1.460  overhead <= 60: pass
n3o4/mul_const                  21.9           188.1         166.2   8.587  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       36.2            24.4         -11.8   0.673  overhead <= 60: pass
n3o4/mul!                      164.6           153.0         -11.6   0.930  overhead <= 60: pass
n3o4/exp!                      711.5           716.9           5.4   1.008  overhead <= 60: pass
n3o4/exp                       711.5           885.2         173.7   1.244  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      744.8           981.3         236.5   1.317  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2788.5          2945.1         156.6   1.056  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition              9964.7         10326.0         361.3   1.036  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       410.9           739.6         328.7   1.800  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11532.4         11987.0         454.6   1.039  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      194.5           200.6           6.1   1.031  overhead <= 60: pass
n6o6/mul_const                 307.3           626.6         319.3   2.039  overhead <= 400 or ratio <= 1.5: pass
n6o6/add!                      410.9           271.8        -139.1   0.661  overhead <= 60: pass
n6o6/mul!                    11532.4         11429.9        -102.5   0.991  ratio <= 1.10: pass
n6o6/exp!                    68981.8         70884.6        1902.8   1.028  ratio <= 1.10: pass
n6o6/exp                     68981.8         71625.4        2643.6   1.038  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    48002.0         48955.2         953.2   1.020  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   177262.7        182346.8        5084.1   1.029  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24733721.0      24761990.0       28269.0   1.001  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5897.8          7245.6        1347.8   1.229  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   445564.8        453344.8        7780.0   1.017  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2354.1          2331.3         -22.8   0.990  ratio <= 1.10: pass
n6o10/mul_const               3551.2          5169.7        1618.5   1.456  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5897.8          3595.4       -2302.4   0.610  ratio <= 1.10: pass
n6o10/mul!                  445564.8        445693.1         128.3   1.000  ratio <= 1.10: pass
n6o10/exp!                 4482737.5       4570790.5       88053.0   1.020  ratio <= 1.10: pass
n6o10/exp                  4482737.5       4535911.5       53174.0   1.012  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1836516.0       1842864.6        6348.6   1.003  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9902260.0       9706133.0     -196127.0   0.980  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      11749810080.0   11715206380.0   -34603700.0   0.997  overhead <= 400 or ratio <= 1.5: pass
gate: PASS
```

### Run 2 (started 08:09:06, core 5 busy 8%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        35.2           200.6         165.4   5.699  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       164.4           333.3         168.9   2.028  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.1           5.5   1.436  overhead <= 60: pass
n3o4/mul_const                  21.6           181.2         159.6   8.390  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       35.2            24.3         -10.9   0.691  overhead <= 60: pass
n3o4/mul!                      164.4           152.1         -12.3   0.925  overhead <= 60: pass
n3o4/exp!                      702.9           712.0           9.1   1.013  overhead <= 60: pass
n3o4/exp                       702.9           878.4         175.5   1.250  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      744.7           976.6         231.9   1.311  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2775.4          3005.0         229.6   1.083  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition             10055.7         10206.4         150.7   1.015  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       406.8           749.9         343.1   1.843  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11536.0         12068.7         532.7   1.046  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      194.5           201.2           6.7   1.034  overhead <= 60: pass
n6o6/mul_const                 303.5           726.8         423.3   2.395  overhead <= 400 or ratio <= 1.5: FAIL
n6o6/add!                      406.8           268.6        -138.2   0.660  overhead <= 60: pass
n6o6/mul!                    11536.0         11441.9         -94.1   0.992  ratio <= 1.10: pass
n6o6/exp!                    69829.8         71009.4        1179.6   1.017  ratio <= 1.10: pass
n6o6/exp                     69829.8         70831.9        1002.1   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    48463.8         48694.5         230.7   1.005  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   176827.0        182103.3        5276.3   1.030  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24654577.0      24824903.0      170326.0   1.007  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5975.9          7479.1        1503.2   1.252  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   452128.0        460119.3        7991.3   1.018  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2318.2          2307.9         -10.3   0.996  ratio <= 1.10: pass
n6o10/mul_const               3653.2          5281.4        1628.2   1.446  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5975.9          3622.9       -2353.0   0.606  ratio <= 1.10: pass
n6o10/mul!                  452128.0        450490.5       -1637.5   0.996  ratio <= 1.10: pass
n6o10/exp!                 4604988.5       4580125.0      -24863.5   0.995  ratio <= 1.10: pass
n6o10/exp                  4604988.5       4609683.0        4694.5   1.001  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1819209.4       1863496.5       44287.1   1.024  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9849140.0       9890634.0       41494.0   1.004  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      11775328717.0   11679257555.0   -96071162.0   0.992  overhead <= 400 or ratio <= 1.5: pass
gate: FAIL
```

### Run 3 (started 08:21:03, core 5 busy 11%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        36.0           198.2         162.2   5.507  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       165.1           335.8         170.7   2.034  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.6            18.1           5.5   1.440  overhead <= 60: pass
n3o4/mul_const                  22.7           184.3         161.6   8.117  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       36.0            24.7         -11.3   0.687  overhead <= 60: pass
n3o4/mul!                      165.1           152.1         -13.0   0.921  overhead <= 60: pass
n3o4/exp!                      709.1           724.5          15.4   1.022  overhead <= 60: pass
n3o4/exp                       709.1           867.3         158.2   1.223  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      745.0           958.2         213.2   1.286  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2739.2          2945.1         205.9   1.075  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition             10043.7         10329.1         285.4   1.028  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       409.9           714.8         304.9   1.744  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11580.3         11954.7         374.4   1.032  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      193.7           203.1           9.4   1.048  overhead <= 60: pass
n6o6/mul_const                 302.4           706.5         404.1   2.336  overhead <= 400 or ratio <= 1.5: FAIL
n6o6/add!                      409.9           266.4        -143.5   0.650  overhead <= 60: pass
n6o6/mul!                    11580.3         11385.6        -194.7   0.983  ratio <= 1.10: pass
n6o6/exp!                    70788.1         70204.5        -583.6   0.992  ratio <= 1.10: pass
n6o6/exp                     70788.1         71504.8         716.7   1.010  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    47894.2         48307.9         413.7   1.009  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   177173.6        180992.2        3818.6   1.022  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          24406592.0      24696319.0      289727.0   1.012  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5765.6          7588.8        1823.2   1.316  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   445676.6        450457.2        4780.7   1.011  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2334.1          2320.5         -13.6   0.994  ratio <= 1.10: pass
n6o10/mul_const               3598.2          5286.9        1688.7   1.469  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5765.6          3934.5       -1831.1   0.682  ratio <= 1.10: pass
n6o10/mul!                  445676.6        448060.5        2383.9   1.005  ratio <= 1.10: pass
n6o10/exp!                 4442847.0       4394285.0      -48562.0   0.989  ratio <= 1.10: pass
n6o10/exp                  4442847.0       4578305.5      135458.5   1.030  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1838679.1       1841182.4        2503.3   1.001  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9696359.0       9629453.0      -66906.0   0.993  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      12346790689.0   12153145482.0  -193645207.0   0.984  overhead <= 400 or ratio <= 1.5: pass
gate: FAIL
```

### Run 4 (started 08:32:31, core 5 busy 22%)

```
case                          C++ ns        Julia ns   overhead ns   ratio  gate
n3o4/add                        34.2           189.9         155.7   5.553  overhead <= 400 or ratio <= 1.5: pass
n3o4/mul                       158.4           327.0         168.6   2.064  overhead <= 400 or ratio <= 1.5: pass
n3o4/iadd                       12.3            17.9           5.6   1.452  overhead <= 60: pass
n3o4/mul_const                  20.9           177.1         156.2   8.472  overhead <= 400 or ratio <= 1.5: pass
n3o4/add!                       34.2            24.7          -9.5   0.721  overhead <= 60: pass
n3o4/mul!                      158.4           148.9          -9.5   0.940  overhead <= 60: pass
n3o4/exp!                      681.8           703.3          21.5   1.032  overhead <= 60: pass
n3o4/exp                       681.8           875.6         193.8   1.284  overhead <= 400 or ratio <= 1.5: pass
n3o4/cmul                      714.8           963.9         249.1   1.348  overhead <= 400 or ratio <= 1.5: pass
n3o4/cexp                     2653.0          2905.2         252.2   1.095  overhead <= 400 or ratio <= 1.5: pass
n3o4/composition              9883.7         10068.3         184.6   1.019  overhead <= 400 or ratio <= 1.5: pass
n6o6/add                       399.0           693.9         294.9   1.739  overhead <= 400 or ratio <= 1.5: pass
n6o6/mul                     11196.1         11641.3         445.2   1.040  overhead <= 400 or ratio <= 1.5: pass
n6o6/iadd                      188.5           197.4           8.9   1.047  overhead <= 60: pass
n6o6/mul_const                 295.8           675.9         380.1   2.285  overhead <= 400 or ratio <= 1.5: pass
n6o6/add!                      399.0           261.0        -138.0   0.654  overhead <= 60: pass
n6o6/mul!                    11196.1         11126.7         -69.4   0.994  ratio <= 1.10: pass
n6o6/exp!                    69035.3         68714.8        -320.5   0.995  ratio <= 1.10: pass
n6o6/exp                     69035.3         68798.4        -236.9   0.997  overhead <= 400 or ratio <= 1.5: pass
n6o6/cmul                    46303.7         46948.0         644.3   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o6/cexp                   170434.1        175858.8        5424.7   1.032  overhead <= 400 or ratio <= 1.5: pass
n6o6/composition          23932567.0      23755553.0     -177014.0   0.993  overhead <= 400 or ratio <= 1.5: pass
n6o10/add                     5791.6          6868.4        1076.8   1.186  overhead <= 400 or ratio <= 1.5: pass
n6o10/mul                   425972.9        434008.5        8035.6   1.019  overhead <= 400 or ratio <= 1.5: pass
n6o10/iadd                    2253.6          2215.3         -38.3   0.983  ratio <= 1.10: pass
n6o10/mul_const               3454.8          4843.4        1388.6   1.402  overhead <= 400 or ratio <= 1.5: pass
n6o10/add!                    5791.6          3473.5       -2318.1   0.600  ratio <= 1.10: pass
n6o10/mul!                  425972.9        423026.0       -2946.9   0.993  ratio <= 1.10: pass
n6o10/exp!                 4419638.0       4253693.8     -165944.2   0.962  ratio <= 1.10: pass
n6o10/exp                  4419638.0       4262661.0     -156977.0   0.964  overhead <= 400 or ratio <= 1.5: pass
n6o10/cmul                 1750160.9       1722658.4      -27502.5   0.984  overhead <= 400 or ratio <= 1.5: pass
n6o10/cexp                 9445844.0       9574914.0      129070.0   1.014  overhead <= 400 or ratio <= 1.5: pass
n6o10/composition      10294772452.0   10231531070.0   -63241382.0   0.994  overhead <= 400 or ratio <= 1.5: pass
gate: PASS
```

## Stage 4 rerun (CNDA cases `cmul`, `cexp`): superseded: FAILED on `n3o4/cexp` only

Superseded by the rerun above (C++ pool change efc0f70 and the 2026-10-01 gate).

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
