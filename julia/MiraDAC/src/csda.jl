# csda.jl — complex symbolic DA vectors, CSDAList and cd_composition for SDA (plan T6.2).

"""
    CSDA(re::SDA, im::SDA)
    CSDA(re::SDA)
    CSDA()
    CSDA(z::Number)

A complex symbolic DA vector `re + im*i` (C++ `std::complex<SDA>`), with parts of one env;
`CSDA()` is zero and `z` gives a constant in the current env. Not a `Number`. Operators combine
it with `CSDA`, `SDA`, `SymExpr`, `Real` and `Complex` on either side (an `Integer` is an exact
`SymExpr`), and an `SDA` with a `Complex` gives a `CSDA`. Its math functions are those of `CNDA`:
`sqrt exp log asin acos atan asinh acosh atanh` and `^`; `abs` is the modulus
`sqrt(re^2 + im^2)` as an `SDA`. `evaluate` gives a `CNDA`, `promote_sda(::CNDA)` a `CSDA`.
"""
mutable struct CSDA
    ptr::Handle
    function CSDA(p::Handle)
        v = new(p)
        return finalizer(csda_finalizer, v)
    end
end

csda_finalizer(v::CSDA) = defer_free(CSDA_QUEUE, v)

Base.unsafe_convert(::Type{Handle}, v::CSDA) = v.ptr
op_env(v::CSDA) = mdac_csda_env(v)

# A new CSDA from the allocating C call f(out) in the env of x; d as in alloc_call.
@inline function new_csda(f, x, d=nothing)
    drain!(CSDA_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x, d)
    return CSDA(out[])
end

# An in-place C call f() into out; it takes temporary slots, so it needs the retry.
@inline function csda_into(f, out, d=nothing)
    drain!(CSDA_QUEUE, out)
    alloc_call(f, out, d)
    return out
end

CSDA(re::SDA, im::SDA) = new_csda(o -> mdac_csda_new(re, im, o), re)
CSDA(re::SDA) = new_csda(o -> mdac_csda_new(re, C_NULL, o), re)
CSDA(; env::Union{DAEnv,Nothing}=nothing) = CSDA(SDA(; env))

function CSDA(z::Number; env::Union{DAEnv,Nothing}=nothing)
    e = env_handle(env)
    return new_csda(o -> mdac_csda_new_z(e, real(z), imag(z), o), e)
end

"""`promote_sda(v::CNDA)`: `v` as a `CSDA` of its env (both parts promoted)."""
promote_sda(v::CNDA) = new_csda(o -> mdac_csda_promote(v, o), v)
promote_sda(e::DAEnv, v::CNDA) = CSDA(promote_sda(e, real(v)), promote_sda(e, imag(v)))

Base.copy(v::CSDA) = new_csda(o -> mdac_csda_copy(v, o), v)
Base.deepcopy_internal(v::CSDA, d::IdDict) = get!(() -> copy(v), d, v)::CSDA

env(v::CSDA) = DAEnv(mdac_csda_env(v))
import_vec(e::DAEnv, v::CSDA) = new_csda(o -> mdac_csda_import(e, v, o), e)

"""`real(v::CSDA)`: a copy of the real part."""
Base.real(v::CSDA) = new_sda(o -> mdac_csda_real(v, o), v)
"""`imag(v::CSDA)`: a copy of the imaginary part."""
Base.imag(v::CSDA) = new_sda(o -> mdac_csda_imag(v, o), v)

function to_string(v::CSDA)
    n = Ref{Csize_t}(0)
    check(mdac_csda_to_string(v, C_NULL, 0, n))
    buf = Vector{UInt8}(undef, n[])
    check(mdac_csda_to_string(v, buf, length(buf), n))
    return String(resize!(buf, length(buf) - 1))
end

function Base.show(io::IO, v::CSDA)
    e = env(v)
    if e.retired
        print(io, "CSDA(cleared env)")
    else
        print(io, "CSDA(order=", e.order, ", nvars=", e.nvars, ", nonzero=(",
              nterms(real(v)), ", ", nterms(imag(v)), "))")
    end
end

function Base.show(io::IO, m::MIME"text/plain", v::CSDA)
    env(v).retired ? show(io, v) : print(io, to_string(v))
end

# ---- Arithmetic -----------------------------------------------------------------------------

# Every operand pair of miradac.h, in the order of csda_op_shapes, as CNDA_OP_ARGS.
const CSDA_OP_ARGS = (
    ((:CSDA, :CSDA), (:a, :b)), ((:CSDA, :SDA), (:a, :b)), ((:SDA, :CSDA), (:a, :b)),
    ((:CSDA, :SymExpr), (:a, :b)), ((:SymExpr, :CSDA), (:a, :b)),
    ((:CSDA, :Real), (:a, :b)), ((:Real, :CSDA), (:a, :b)),
    ((:CSDA, :Complex), (:a, :(real(b)), :(imag(b)))), ((:Complex, :CSDA), (:(real(a)), :(imag(a)), :b)),
    ((:SDA, :Complex), (:a, :(real(b)), :(imag(b)))), ((:Complex, :SDA), (:(real(a)), :(imag(a)), :b)))

for (op, name, bang) in ((:+, :add, :add!), (:-, :sub, :sub!), (:*, :mul, :mul!), (:/, :div, :div!))
    target = bang === :mul! ? :(LinearAlgebra.mul!) : bang
    for ((f, _), ((A, B), cargs)) in zip(csda_op_shapes(name), CSDA_OP_ARGS)
        fi = Symbol(f, :_into)
        da = A in (:Real, :Complex, :SymExpr) ? :b : :a     # the DA operand gives the env
        @eval begin
            Base.$op(a::$A, b::$B) = new_csda(o -> $f($(cargs...), o), $da)
            $target(out::CSDA, a::$A, b::$B) = csda_into(() -> $fi(out, $(cargs...)), out)
        end
    end
    # An Integer is an exact SymExpr integer, as for SDA.
    @eval begin
        Base.$op(a::CSDA, x::Integer) = $op(a, SymExpr(x))
        Base.$op(x::Integer, a::CSDA) = $op(SymExpr(x), a)
        $target(out::CSDA, a::CSDA, x::Integer) = $target(out, a, SymExpr(x))
        $target(out::CSDA, x::Integer, a::CSDA) = $target(out, SymExpr(x), a)
    end
end

Base.:-(a::CSDA) = new_csda(o -> mdac_csda_neg(a, o), a)
Base.:^(a::CSDA, n::Integer) = new_csda(o -> mdac_csda_pow_i(a, n, o), a, a)
Base.:^(a::CSDA, x::Real) = new_csda(o -> mdac_csda_pow_d(a, x, o), a, a)

# ---- Math functions and evaluation ----------------------------------------------------------

for fn in CNDA_FUNCS
    f, fi, bang = Symbol(:mdac_csda_, fn), Symbol(:mdac_csda_, fn, :_into), Symbol(fn, :!)
    @eval begin
        Base.$fn(a::CSDA) = new_csda(o -> $f(a, o), a, a)
        $bang(out::CSDA, a::CSDA) = csda_into(() -> $fi(out, a), out, a)
    end
end

"""`abs(v::CSDA)`: the modulus `sqrt(real(v)^2 + imag(v)^2)`, an `SDA`."""
Base.abs(a::CSDA) = new_sda(o -> mdac_csda_abs(a, o), a, a)

"""
    evaluate(z::CSDA, d::AbstractDict)

The `CNDA` of `z` with the numbers `d[k]` for the symbols `k`, as `evaluate(::SDA, d)`.
"""
function evaluate(z::CSDA, d::AbstractDict{<:Any,<:Real})
    ks = SymExpr[SymExpr(k) for k in keys(d)]
    kp, vals = Handle[k.ptr for k in ks], Float64[v for v in values(d)]
    GC.@preserve ks new_cnda(o -> mdac_csda_evaluate(z, length(kp), kp, vals, o), z)
end

# ---- CSDAList and composition ---------------------------------------------------------------

"""
    CSDAList()
    CSDAList(v::AbstractVector{CSDA})

A list of `CSDA` held in C++; it behaves as `NDAList`.
"""
mutable struct CSDAList <: AbstractVector{CSDA}
    ptr::Handle
    function CSDAList(p::Handle)
        l = new(p)
        return finalizer(cslist_finalizer, l)
    end
end

cslist_finalizer(l::CSDAList) = defer_free(CSLIST_QUEUE, l)

Base.unsafe_convert(::Type{Handle}, l::CSDAList) = l.ptr
op_env(l::CSDAList) = mdac_csdalist_env(l)

function new_cslist(f, x)
    drain!(CSLIST_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x)
    return CSDAList(out[])
end

function CSDAList()
    out = Ref{Handle}(C_NULL)
    check(mdac_csdalist_new(out))
    return CSDAList(out[])
end

function CSDAList(v::AbstractVector{CSDA})
    isempty(v) && return CSDAList()
    ptrs = Handle[x.ptr for x in v]
    GC.@preserve v new_cslist(o -> mdac_csdalist_from(ptrs, length(ptrs), o), first(v))
end

aslist(l::CSDAList) = l
aslist(v::AbstractVector{CSDA}) = CSDAList(v)

Base.size(l::CSDAList) = (Int(mdac_csdalist_length(l)),)

function Base.getindex(l::CSDAList, i::Int)
    @boundscheck checkbounds(l, i)
    return new_csda(o -> mdac_csdalist_get(l, i - 1, o), l)
end

function Base.setindex!(l::CSDAList, v::CSDA, i::Int)
    @boundscheck checkbounds(l, i)
    check(mdac_csdalist_set(l, i - 1, v))
    return l
end

function Base.push!(l::CSDAList, v::CSDA)
    alloc_call(() -> mdac_csdalist_push(l, v), v)
    return l
end

function env(l::CSDAList)
    e = mdac_csdalist_env(l)
    e == C_NULL && throw(EnvError("an empty CSDAList has no env"))
    return DAEnv(e)
end

# C++ cd_composition for T = Expression, as a new CSDAList.
for (M, V, f) in ((:SDA, :CSDA, :mdac_sdalist_compose_c), (:CSDA, :CSDA, :mdac_csdalist_compose),
                  (:CSDA, :SDA, :mdac_csdalist_compose_s))
    @eval function compose(m::AbstractVector{$M}, v::AbstractVector{$V})
        isempty(m) && return CSDAList()
        ml, vl = aslist(m), aslist(v)
        return new_cslist(o -> $f(ml, vl, o), ml)
    end
end
