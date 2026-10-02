# cnda.jl — complex numeric DA vectors, CNDAList and cd_composition (plan T4.2).

"""
    CNDA(re::NDA, im::NDA)
    CNDA(re::NDA)
    CNDA(z::Number)

A complex numeric DA vector `re + im*i` (C++ `std::complex<NDA>`), with parts of one env; `z`
gives a constant in the current env. Not a `Number`. Operators combine it with `CNDA`, `NDA`,
`Real` and `Complex` on either side, and an `NDA` with a `Complex` gives a `CNDA`. Its math
functions are those C++ has for complex vectors: `sqrt exp log asin acos atan asinh acosh atanh`,
`^`, and `abs` (the larger of the parts' `norm`s). There is no `conj` (C++ has none).
"""
mutable struct CNDA
    ptr::Handle
    function CNDA(p::Handle)
        return adopt!(new(p))     # a finalizer, or the current dascope
    end
end

cnda_finalizer(v::CNDA) = defer_free(CNDA_QUEUE, v)

Base.unsafe_convert(::Type{Handle}, v::CNDA) = handle(v)
op_env(v::CNDA) = mdac_cnda_env(v)

# A new CNDA from the allocating C call f(out) in the env of x; d as in alloc_call.
@inline function new_cnda(f, x, d=nothing)
    drain!(CNDA_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x, d)
    return CNDA(out[])
end

# An in-place C call f() into out. The C++ complex operators and functions return new vectors,
# so these take temporary slots too and need the retry of alloc_call.
@inline function into_call(f, out, d=nothing)
    drain!(CNDA_QUEUE, out)
    alloc_call(f, out, d)
    return out
end

CNDA(re::NDA, im::NDA) = new_cnda(o -> mdac_cnda_new(re, im, o), re)
CNDA(re::NDA) = new_cnda(o -> mdac_cnda_new(re, C_NULL, o), re)

function CNDA(z::Number; env::Union{DAEnv,Nothing}=nothing)
    e = env_handle(env)
    return new_cnda(o -> mdac_cnda_new_z(e, real(z), imag(z), o), e)
end

Base.copy(v::CNDA) = new_cnda(o -> mdac_cnda_copy(v, o), v)
Base.deepcopy_internal(v::CNDA, d::IdDict) = get!(() -> copy(v), d, v)::CNDA

env(v::CNDA) = DAEnv(mdac_cnda_env(v))
import_vec(e::DAEnv, v::CNDA) = new_cnda(o -> mdac_cnda_import(e, v, o), e)

"""`real(v::CNDA)`: a copy of the real part."""
Base.real(v::CNDA) = new_nda(o -> mdac_cnda_real(v, o), v)
"""`imag(v::CNDA)`: a copy of the imaginary part."""
Base.imag(v::CNDA) = new_nda(o -> mdac_cnda_imag(v, o), v)

function to_string(v::CNDA)
    n = Ref{Csize_t}(0)
    check(mdac_cnda_to_string(v, C_NULL, 0, n))
    buf = Vector{UInt8}(undef, n[])
    check(mdac_cnda_to_string(v, buf, length(buf), n))
    return String(resize!(buf, length(buf) - 1))
end

function Base.show(io::IO, v::CNDA)
    getfield(v, :ptr) == C_NULL && return print(io, "CNDA(freed)")   # by its dascope
    e = env(v)
    if e.retired
        print(io, "CNDA(cleared env)")
    else
        print(io, "CNDA(order=", e.order, ", nvars=", e.nvars, ", nonzero=(",
              nterms(real(v)), ", ", nterms(imag(v)), "))")
    end
end

function Base.show(io::IO, m::MIME"text/plain", v::CNDA)
    getfield(v, :ptr) == C_NULL && return print(io, "CNDA(freed)")   # by its dascope
    env(v).retired ? show(io, v) : print(io, to_string(v))
end

# ---- Arithmetic -----------------------------------------------------------------------------

# Every operand pair of miradac.h, in the order of cnda_op_shapes: Julia argument types, and
# how each argument reaches C (a complex number as two doubles).
const CNDA_OP_ARGS = (
    ((:CNDA, :CNDA), (:a, :b)), ((:CNDA, :NDA), (:a, :b)), ((:NDA, :CNDA), (:a, :b)),
    ((:CNDA, :Real), (:a, :b)), ((:Real, :CNDA), (:a, :b)),
    ((:CNDA, :Complex), (:a, :(real(b)), :(imag(b)))), ((:Complex, :CNDA), (:(real(a)), :(imag(a)), :b)),
    ((:NDA, :Complex), (:a, :(real(b)), :(imag(b)))), ((:Complex, :NDA), (:(real(a)), :(imag(a)), :b)))

for (op, name, bang) in ((:+, :add, :add!), (:-, :sub, :sub!), (:*, :mul, :mul!), (:/, :div, :div!))
    target = bang === :mul! ? :(LinearAlgebra.mul!) : bang
    for ((f, _), ((A, B), cargs)) in zip(cnda_op_shapes(name), CNDA_OP_ARGS)
        fi = Symbol(f, :_into)
        da = A === :Real || A === :Complex ? :b : :a     # the DA operand gives the env
        @eval begin
            Base.$op(a::$A, b::$B) = new_cnda(o -> $f($(cargs...), o), $da)
            $target(out::CNDA, a::$A, b::$B) = into_call(() -> $fi(out, $(cargs...)), out)
        end
    end
end

Base.:-(a::CNDA) = new_cnda(o -> mdac_cnda_neg(a, o), a)
Base.:^(a::CNDA, n::Integer) = new_cnda(o -> mdac_cnda_pow_i(a, n, o), a, a)
Base.:^(a::CNDA, x::Real) = new_cnda(o -> mdac_cnda_pow_d(a, x, o), a, a)

# ---- Math functions -------------------------------------------------------------------------

for fn in CNDA_FUNCS
    f, fi, bang = Symbol(:mdac_cnda_, fn), Symbol(:mdac_cnda_, fn, :_into), Symbol(fn, :!)
    @eval begin
        Base.$fn(a::CNDA) = new_cnda(o -> $f(a, o), a, a)
        $bang(out::CNDA, a::CNDA) = into_call(() -> $fi(out, a), out, a)
    end
end

"""`abs(v::CNDA)`: the larger of `norm(real(v))` and `norm(imag(v))`, as in C++."""
function Base.abs(a::CNDA)
    r = Ref{Cdouble}()
    check(mdac_cnda_abs(a, r))
    return r[]
end

# ---- CNDAList and composition ---------------------------------------------------------------

"""
    CNDAList()
    CNDAList(v::AbstractVector{CNDA})

A list of `CNDA` held in C++; it behaves as `NDAList`.
"""
mutable struct CNDAList <: AbstractVector{CNDA}
    ptr::Handle
    function CNDAList(p::Handle)
        return adopt!(new(p))     # a finalizer, or the current dascope
    end
end

clist_finalizer(l::CNDAList) = defer_free(CLIST_QUEUE, l)

Base.unsafe_convert(::Type{Handle}, l::CNDAList) = handle(l)
op_env(l::CNDAList) = mdac_cndalist_env(l)

function new_clist(f, x)
    drain!(CLIST_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x)
    return CNDAList(out[])
end

function CNDAList()
    out = Ref{Handle}(C_NULL)
    check(mdac_cndalist_new(out))
    return CNDAList(out[])
end

function CNDAList(v::AbstractVector{CNDA})
    isempty(v) && return CNDAList()
    ptrs = Handle[handle(x) for x in v]
    GC.@preserve v new_clist(o -> mdac_cndalist_from(ptrs, length(ptrs), o), first(v))
end

aslist(l::CNDAList) = l
aslist(v::AbstractVector{CNDA}) = CNDAList(v)

Base.size(l::CNDAList) = (Int(mdac_cndalist_length(l)),)

function Base.getindex(l::CNDAList, i::Int)
    @boundscheck checkbounds(l, i)
    return new_cnda(o -> mdac_cndalist_get(l, i - 1, o), l)
end

function Base.setindex!(l::CNDAList, v::CNDA, i::Int)
    @boundscheck checkbounds(l, i)
    check(mdac_cndalist_set(l, i - 1, v))
    return l
end

function Base.push!(l::CNDAList, v::CNDA)
    alloc_call(() -> mdac_cndalist_push(l, v), v)
    return l
end

function env(l::CNDAList)
    e = mdac_cndalist_env(l)
    e == C_NULL && throw(EnvError("an empty CNDAList has no env"))
    return DAEnv(e)
end

# C++ cd_composition, as a new CNDAList.
for (M, V, f) in ((:NDA, :CNDA, :mdac_ndalist_compose_c), (:CNDA, :CNDA, :mdac_cndalist_compose),
                  (:CNDA, :NDA, :mdac_cndalist_compose_n))
    @eval function compose(m::AbstractVector{$M}, v::AbstractVector{$V})
        isempty(m) && return CNDAList()
        ml, vl = aslist(m), aslist(v)
        return new_clist(o -> $f(ml, vl, o), ml)
    end
end
