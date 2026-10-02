# Julia counterpart of bench_cpp (plan T3.4), measured the same way as python/bench/bench_ops.py.
#
# Each case first doubles its repetition count until two samples in a row take >= 10 ms (the
# warm-up; bench_ops.py needs one, but here a garbage collection can make one sample long).
# Then ROUNDS rounds each give every case 0.2 s of samples (at least one), and a case reports
# its minimum sample, in thread CPU time (garbage collections inside a sample count). The C++
# numbers come from a fresh `bench_cpp --driven` child (same CPU pinning, inherited), which
# runs one round before each Julia round. Prints the table and the gate of plan A.9, and exits
# 1 if the gate fails.
#
# With --scoped, the allocating cases also run inside `dascope`, one scope per batch of
# --batch N operations (default 16), as a user would wrap a step of a computation; their
# names start with `scoped_` and they get the strict rule (plan A.5b).
#
# Run: taskset -c 3 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl \
#          [--rounds N] [--cpp PATH] [--scoped] [--batch N]

using MiraDAC
using Printf
using Random

const SAMPLE_NS = 10_000_000
const PER_ROUND_NS = 200_000_000
# (nvars, order, pool size): plan A.9, a realistic pool (a slot at (6,10) is 64 KB).
const SIZES = [(3, 4, 10_000), (6, 6, 10_000), (6, 10, 2_000)]

const CLOCK_THREAD_CPUTIME_ID = Cint(3)

function cpu_ns()
    ts = Ref{NTuple{2,Clong}}()
    ccall(:clock_gettime, Cint, (Cint, Ptr{NTuple{2,Clong}}), CLOCK_THREAD_CPUTIME_ID, ts)
    return Int(ts[][1]) * 1_000_000_000 + Int(ts[][2])
end

# One sample: reps repetitions of the case; each returns the elapsed thread CPU time.
function t_add(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        c = a + b
    end
    return cpu_ns() - t0
end

function t_mul(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        c = a * b
    end
    return cpu_ns() - t0
end

# C++ `c += b`: in place.
function t_iadd(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        add!(c, c, b)
    end
    return cpu_ns() - t0
end

function t_mul_const(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        c = a * 2.0
    end
    return cpu_ns() - t0
end

function t_exp(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        c = exp(a)
    end
    return cpu_ns() - t0
end

function t_composition(m, n, out, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        out = compose(m, n)
    end
    return cpu_ns() - t0
end

function t_add!(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        add!(c, a, b)
    end
    return cpu_ns() - t0
end

function t_mul!(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        mul!(c, a, b)
    end
    return cpu_ns() - t0
end

function t_exp!(a, b, c, reps)
    t0 = cpu_ns()
    for _ in 1:reps
        exp!(c, a)
    end
    return cpu_ns() - t0
end

# The allocating cases inside dascope, one scope per batch of BATCH[] operations.
const BATCH = Ref(16)

function t_scoped(op::F, a, b, reps) where {F}
    t0 = cpu_ns()
    i = 0
    while i < reps
        n = min(BATCH[], reps - i)
        dascope() do
            for _ in 1:n
                op(a, b)
            end
        end
        i += n
    end
    return cpu_ns() - t0
end

const SCOPED = Ref(false)
const SCOPED_OPS = (("add", +), ("mul", *), ("mul_const", (a, b) -> a * 2.0), ("exp", (a, b) -> exp(a)))

const reps_of = Dict{String,Int}()
const best_ns = Dict{String,Float64}()
const names = String[]                 # in first-run order
const errors = Dict{String,String}()   # a case that threw, and why

function bench(name, fn, args...)
    haskey(errors, name) && return
    try
        reps = get(reps_of, name, 0)
        if reps == 0
            push!(names, name)
            # Compile first, and require two samples of >= 10 ms: one long garbage collection
            # must not set a small reps, whose samples then mostly skip the collections.
            fn(args..., 1)
            reps = 1
            while min(fn(args..., reps), fn(args..., reps)) < SAMPLE_NS
                reps *= 2
            end
            reps_of[name] = reps
            best_ns[name] = Inf
        end
        spent = 0
        while spent < PER_ROUND_NS
            t = fn(args..., reps)
            best_ns[name] = min(best_ns[name], t / reps)
            spent += t
        end
    catch e
        errors[name] = sprint(showerror, e)
    end
end

# One size: the unscoped cases, or (scoped = true) the scoped_* cases. A round runs every size
# unscoped first, then every size scoped, so the unscoped cases follow the same work as in a
# run without --scoped (interleaved, scoped cases made the next unscoped n6o6 cases ~40 ns slower).
function numeric(nvars, order, poolsize, scoped)
    init!(order, nvars, poolsize; table=true)
    rng = MersenneTwister(12345)       # same data every round
    rand_nda() = NDA(2 .* rand(rng, current_env().full_length) .- 1)
    p = "n$(nvars)o$(order)/"
    a, b, c = rand_nda(), rand_nda(), NDA(0.0)
    if !scoped                         # the cases and data in the order of a plain run
        bench(p * "add", t_add, a, b, c)
        bench(p * "mul", t_mul, a, b, c)
        bench(p * "iadd", t_iadd, a, b, c)
        bench(p * "mul_const", t_mul_const, a, b, c)
        bench(p * "add!", t_add!, a, b, c)
        bench(p * "mul!", t_mul!, a, b, c)
        bench(p * "exp!", t_exp!, a, b, c)
        bench(p * "exp", t_exp, a, b, c)
    else
        for (name, op) in SCOPED_OPS
            bench(p * "scoped_" * name, t_scoped, op, a, b)
        end
    end
    ca, cb, cc = CNDA(a, b), CNDA(b, a), CNDA(0.0)
    if !scoped
        bench(p * "cmul", t_mul, ca, cb, cc)
        bench(p * "cexp", t_exp, ca, cb, cc)
    else
        bench(p * "scoped_cmul", t_scoped, *, ca, cb)
        bench(p * "scoped_cexp", t_scoped, (a, b) -> exp(a), ca, cb)
    end
    m, n = NDAList(), NDAList()
    for _ in 1:nvars
        push!(m, rand_nda())
        push!(n, rand_nda())
    end
    scoped ? bench(p * "scoped_composition", t_scoped, compose, m, n) :
             bench(p * "composition", t_composition, m, n, nothing)
    clear!()
end

# The C++ case a Julia case is compared with.
cpp_case(name) = replace(name, "!" => "", "scoped_" => "")

const IN_PLACE = ("iadd", "add!", "mul!", "exp!")

# Plan A.9: in-place cases get one rule by C++ time; allocating cases pass if
# their overhead is <= 400 ns or their ratio is <= 1.5 (user decision 2026-10-01).
function gate(name, cpp, jl)
    # Plan A.5b: allocating operations inside dascope get the strict rule.
    if startswith(split(name, "/")[2], "scoped_")
        cpp >= 1000 && return "ratio <= 1.10", jl / cpp <= 1.10
        return "overhead <= 150", jl - cpp <= 150
    end
    if split(name, "/")[2] in IN_PLACE
        cpp >= 1000 && return "ratio <= 1.10", jl / cpp <= 1.10
        return "overhead <= 60", jl - cpp <= 60
    end
    return "overhead <= 400 or ratio <= 1.5", jl - cpp <= 400 || jl / cpp <= 1.5
end

function parse_json_numbers(s)
    Dict(m[1] => parse(Float64, m[2]) for m in eachmatch(r"\"([^\"]+)\":\s*([-+0-9.eE]+)", s))
end

function main(args)
    rounds, cpp = 15, normpath(joinpath(@__DIR__, "..", "..", "..", "build", "python", "bench", "bench_cpp"))
    i = 1
    while i <= length(args)
        args[i] == "--rounds" ? (rounds = parse(Int, args[i+1])) :
        args[i] == "--cpp" ? (cpp = args[i+1]) :
        args[i] == "--batch" ? (BATCH[] = parse(Int, args[i+1])) :
        args[i] == "--scoped" ? (SCOPED[] = true; i -= 1) : error("unknown argument $(args[i])")
        i += 2
    end
    Sys.islinux() && ccall(:prctl, Cint, (Cint, Culong, Culong, Culong, Culong), 41, 1, 0, 0, 0)

    child = open(`$cpp --driven`, "r+")
    for r in 1:rounds
        write(child, "\n")
        flush(child)
        readline(child) == "done" || error("$cpp stopped in round $r")
        for scoped in (false, true), (nvars, order, poolsize) in SIZES
            scoped && !SCOPED[] && continue
            numeric(nvars, order, poolsize, scoped)
        end
        println(stderr, "round $r/$rounds done")
    end
    close(child.in)
    out = read(child, String)
    wait(child)
    success(child) || error("$cpp exited with $(child.exitcode)")
    cpp_ns = parse_json_numbers(out)

    ok = Dict("A.9" => true, "A.5b" => true)   # unscoped cases, scoped_* cases
    @printf("%-20s%16s%16s%14s%8s  %s\n", "case", "C++ ns", "Julia ns", "overhead ns", "ratio", "gate")
    for name in names
        g = startswith(split(name, "/")[2], "scoped_") ? "A.5b" : "A.9"
        if haskey(errors, name)
            ok[g] = false
            @printf("%-20s%16.1f%16s%14s%8s  FAIL: %s\n", name, cpp_ns[cpp_case(name)], "error", "", "",
                    errors[name])
            continue
        end
        c, jl = cpp_ns[cpp_case(name)], best_ns[name]
        rule, pass = gate(name, c, jl)
        ok[g] &= pass
        @printf("%-20s%16.1f%16.1f%14.1f%8.3f  %s: %s\n", name, c, jl, jl - c, jl / c, rule,
                pass ? "pass" : "FAIL")
    end
    # Reported apart: the unscoped n6o6 cases have little A.9 margin and fail under load.
    SCOPED[] && println("gate A.5b (scoped_* cases): ", ok["A.5b"] ? "PASS" : "FAIL")
    println("gate A.9 (unscoped cases): ", ok["A.9"] ? "PASS" : "FAIL")
    pass = all(values(ok))
    println("gate: ", pass ? "PASS" : "FAIL")
    return pass ? 0 : 1
end

exit(main(ARGS))
