# Julia counterpart of bench_cpp (plan T3.4), measured the same way as python/bench/bench_ops.py.
#
# Each case first doubles its repetition count until two samples in a row take >= 10 ms (the
# warm-up; bench_ops.py needs one, but here a garbage collection can make one sample long).
# Then ROUNDS rounds each give every case 0.2 s of samples (at least one), and a case reports
# its minimum sample, in thread CPU time (garbage collections inside a sample count). The C++
# numbers come from a fresh `bench_cpp --driven` child (same CPU pinning, inherited), which
# runs one round before each Julia round. Prints the table and the gate of plan A.9: in-place
# cases are blocking, allocating ones informational until `dascope` lands. Exits 1 only if the
# blocking gate fails.
#
# Run: taskset -c 5 julia --project=julia/MiraDAC/bench julia/MiraDAC/bench/bench_ops.jl \
#          [--rounds N] [--cpp PATH]

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

function numeric(nvars, order, poolsize)
    init!(order, nvars, poolsize; table=true)
    rng = MersenneTwister(12345)       # same data every round
    rand_nda() = NDA(2 .* rand(rng, current_env().full_length) .- 1)
    p = "n$(nvars)o$(order)/"
    a, b, c = rand_nda(), rand_nda(), NDA(0.0)
    bench(p * "add", t_add, a, b, c)
    bench(p * "mul", t_mul, a, b, c)
    bench(p * "iadd", t_iadd, a, b, c)
    bench(p * "mul_const", t_mul_const, a, b, c)
    bench(p * "add!", t_add!, a, b, c)
    bench(p * "mul!", t_mul!, a, b, c)
    bench(p * "exp!", t_exp!, a, b, c)
    bench(p * "exp", t_exp, a, b, c)
    ca, cb, cc = CNDA(a, b), CNDA(b, a), CNDA(0.0)
    bench(p * "cmul", t_mul, ca, cb, cc)
    bench(p * "cexp", t_exp, ca, cb, cc)
    m, n = NDAList(), NDAList()
    for _ in 1:nvars
        push!(m, rand_nda())
        push!(n, rand_nda())
    end
    bench(p * "composition", t_composition, m, n, nothing)
    clear!()
end

# The C++ case a Julia case is compared with.
cpp_case(name) = replace(name, "!" => "")

const IN_PLACE = ("iadd", "add!", "mul!", "exp!")

is_in_place(name) = split(name, "/")[2] in IN_PLACE

# Plan A.9: in-place cases get one rule by C++ time; allocating cases pass if
# their overhead is <= 400 ns or their ratio is <= 1.5 (user decision 2026-10-01).
function gate(name, cpp, jl)
    if is_in_place(name)
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
        args[i] == "--cpp" ? (cpp = args[i+1]) : error("unknown argument $(args[i])")
        i += 2
    end
    Sys.islinux() && ccall(:prctl, Cint, (Cint, Culong, Culong, Culong, Culong), 41, 1, 0, 0, 0)

    child = open(`$cpp --driven`, "r+")
    for r in 1:rounds
        write(child, "\n")
        flush(child)
        readline(child) == "done" || error("$cpp stopped in round $r")
        for (nvars, order, poolsize) in SIZES
            numeric(nvars, order, poolsize)
        end
        println(stderr, "round $r/$rounds done")
    end
    close(child.in)
    out = read(child, String)
    wait(child)
    success(child) || error("$cpp exited with $(child.exitcode)")
    cpp_ns = parse_json_numbers(out)

    ok = Dict(true => true, false => true)   # keyed by is_in_place
    @printf("%-20s%16s%16s%14s%8s  %s\n", "case", "C++ ns", "Julia ns", "overhead ns", "ratio", "gate")
    for name in names
        if haskey(errors, name)
            ok[is_in_place(name)] = false
            @printf("%-20s%16.1f%16s%14s%8s  FAIL: %s\n", name, cpp_ns[cpp_case(name)], "error", "", "",
                    errors[name])
            continue
        end
        c, jl = cpp_ns[cpp_case(name)], best_ns[name]
        rule, pass = gate(name, c, jl)
        ok[is_in_place(name)] &= pass
        @printf("%-20s%16.1f%16.1f%14.1f%8.3f  %s: %s\n", name, c, jl, jl - c, jl / c, rule,
                pass ? "pass" : "FAIL")
    end
    println("gate (blocking): ", ok[true] ? "PASS" : "FAIL")
    println("allocating (informational): ", ok[false] ? "pass" : "fail",
            join([" $n" for n in names if !is_in_place(n) &&
                  (haskey(errors, n) || !gate(n, cpp_ns[cpp_case(n)], best_ns[n])[2])], ","))
    return ok[true] ? 0 : 1
end

exit(main(ARGS))
