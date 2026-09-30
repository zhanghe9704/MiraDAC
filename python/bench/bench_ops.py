"""Python counterpart of bench_cpp (plan T2.5), measured the same way.

Each case first doubles its repetition count until one sample takes >= 10 ms
(the warm-up). Then ROUNDS rounds each give every case 0.2 s of samples (at
least one), and a case reports its minimum sample, in thread CPU time. See
bench_cpp.cpp for why. The C++ numbers come from a fresh `bench_cpp --driven`
child (same CPU pinning, inherited), which runs one round before each Python
round. Prints a table and the Stage 2 gate (Stage 3 applies it to the CNDA
cases `cmul` and `cexp` too, Stage 5 to the symbolic cases, with ratio <= 1.05
for `sym_n3o3/exp`), and exits 1 if the gate fails.

Run: taskset -c 5 .venv/bin/python python/bench/bench_ops.py [--rounds N] [--cpp PATH]
"""

import argparse
import ctypes
import json
import math
import subprocess
import sys
from itertools import repeat
from pathlib import Path
from time import thread_time_ns

import numpy as np

import miradac as da

SAMPLE_NS = 10_000_000
PER_ROUND_NS = 200_000_000
SIZES = [(3, 4), (6, 6), (6, 10)]


def t_add(a, b, c, reps):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        c = a + b
    return thread_time_ns() - t0


def t_mul(a, b, c, reps):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        c = a * b
    return thread_time_ns() - t0


def t_iadd(a, b, c, reps):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        c += b
    return thread_time_ns() - t0


def t_mul_const(a, b, c, reps):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        c = a * 2.0
    return thread_time_ns() - t0


def t_exp(a, b, c, reps, exp=da.exp):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        c = exp(a)
    return thread_time_ns() - t0


def t_composition(m, n, out, reps, compose=da.compose):
    t0 = thread_time_ns()
    for _ in repeat(None, reps):
        compose(m, n, out)
    return thread_time_ns() - t0


reps_of = {}
best_ns = {}


def bench(name, fn, *args):
    reps = reps_of.get(name)
    if reps is None:
        reps = 1
        while fn(*args, reps) < SAMPLE_NS:
            reps *= 2
        reps_of[name] = reps
        best_ns[name] = math.inf
    spent = 0
    while spent < PER_ROUND_NS:
        t = fn(*args, reps)
        best_ns[name] = min(best_ns[name], t / reps)
        spent += t


def numeric(nvars, order):
    da.init(order, nvars, 400, True)
    rng = np.random.default_rng(12345)       # same data every round

    def rand():
        return da.NDA.from_coeffs(rng.uniform(-1.0, 1.0, da.full_length()))

    p = f"n{nvars}o{order}/"
    a, b, c = rand(), rand(), da.NDA()
    bench(p + "add", t_add, a, b, c)
    bench(p + "mul", t_mul, a, b, c)
    bench(p + "iadd", t_iadd, a, b, c)
    bench(p + "mul_const", t_mul_const, a, b, c)
    bench(p + "exp", t_exp, a, b, c)
    ca, cb, cc = da.CNDA(a, b), da.CNDA(b, a), da.CNDA()
    bench(p + "cmul", t_mul, ca, cb, cc)
    bench(p + "cexp", t_exp, ca, cb, cc)
    m = da.NDAList()
    n = da.NDAList()
    for _ in range(nvars):
        m.append(rand())
        n.append(rand())
    out = da.NDAList([da.NDA() for _ in range(nvars)])
    bench(p + "composition", t_composition, m, n, out)
    del a, b, c, ca, cb, cc, m, n, out
    da.clear()


def symbolic(nvars, order):
    da.init(order, nvars, 400, True)

    def rand(prefix):                        # every coefficient its own symbol
        v = da.SDA()
        for i, row in enumerate(da.exponents().tolist()):
            v.set_element(row, da.Expr(f"{prefix}{i}"))
        return v

    p = f"sym_n{nvars}o{order}/"
    a, b, c = rand("a"), rand("b"), da.SDA()
    bench(p + "mul", t_mul, a, b, c)
    bench(p + "exp", t_exp, a, b, c)
    del a, b, c
    da.clear()


IN_PLACE = ("iadd", "composition")       # composition writes into `out`


def gate(name, cpp, py):
    """Stage 2 gate (T2.5, second revision): one rule per case, chosen by its
    C++ time. Returns (rule text, passed). T5.6 sets ratio <= 1.05 for
    exp(SDA)."""
    if name == "sym_n3o3/exp":
        return "ratio <= 1.05", py / cpp <= 1.05
    if cpp < 1000:
        limit = 100 if name.split("/")[1] in IN_PLACE else 150
        return f"overhead <= {limit}", py - cpp <= limit
    return "ratio <= 1.10", py / cpp <= 1.10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=15)
    ap.add_argument("--cpp", default=Path(__file__).parents[2] / "build-bench/python/bench/bench_cpp",
                    help="bench_cpp binary (built with -DDA_BUILD_BENCH=ON)")
    args = ap.parse_args()
    if sys.platform == "linux":
        ctypes.CDLL(None).prctl(41, 1, 0, 0, 0)       # PR_SET_THP_DISABLE

    cpp_run = subprocess.Popen([str(args.cpp), "--driven"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    for r in range(args.rounds):
        cpp_run.stdin.write("\n")
        cpp_run.stdin.flush()
        if cpp_run.stdout.readline() != "done\n":
            sys.exit(f"{args.cpp} stopped in round {r + 1}")
        for nvars, order in SIZES:
            numeric(nvars, order)
        symbolic(3, 3)
        print(f"round {r + 1}/{args.rounds} done", file=sys.stderr)
    out, _ = cpp_run.communicate()
    if cpp_run.returncode:
        sys.exit(f"{args.cpp} exited with {cpp_run.returncode}")
    cpp_ns = json.loads(out)

    ok = True
    print(f"{'case':<20}{'C++ ns':>16}{'Python ns':>16}{'overhead ns':>14}{'ratio':>8}  gate")
    for name, py in best_ns.items():
        cpp = cpp_ns[name]
        g = gate(name, cpp, py)
        verdict = f"{g[0]}: {'pass' if g[1] else 'FAIL'}"
        ok = ok and g[1]
        print(f"{name:<20}{cpp:>16.1f}{py:>16.1f}{py - cpp:>14.1f}{py / cpp:>8.3f}  {verdict}")
    print("gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
