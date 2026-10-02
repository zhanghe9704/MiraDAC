# symbolic.jl — symbolic scalars (SymExpr), symbolic DA vectors (SDA), SDAList and the SDA
# algorithms (plan T5.3, T5.4). Needs a library built with symbolic support (HAS_SYMBOLIC);
# otherwise every call throws.

# ---- SymExpr --------------------------------------------------------------------------------

"""
    SymExpr(x::Real)
    SymExpr(s::AbstractString)

A symbolic scalar (C++ `SymEngine::Expression`), the coefficient type of `SDA`. An `Integer`
gives an exact integer, another `Real` a floating-point number, a string the parsed expression
(`SymExpr("a + 2*b")`; `ArgumentError` for bad input). Not a `Number`; `+ - * / ^` combine it
with `SymExpr` and `Real`. See `dasymbols`, `Float64`, `MiraDAC.subs`, `MiraDAC.expand`,
`diff`, `MiraDAC.free_symbols`, `simplify`.
"""
mutable struct SymExpr
    ptr::Handle
    function SymExpr(p::Handle)
        x = new(p)
        return finalizer(expr_finalizer, x)
    end
end

expr_finalizer(x::SymExpr) = defer_free(EXPR_QUEUE, x)

Base.unsafe_convert(::Type{Handle}, x::SymExpr) = x.ptr

function new_expr(f)
    drain!(EXPR_QUEUE, C_NULL)
    out = Ref{Handle}(C_NULL)
    check(f(out))
    return SymExpr(out[])
end

SymExpr(x::SymExpr) = x
SymExpr(x::Integer) = typemin(Int64) <= x <= typemax(Int64) ?
    new_expr(o -> mdac_expr_new_i(Int64(x), o)) : SymExpr(string(x))
SymExpr(x::Real) = new_expr(o -> mdac_expr_new_d(Float64(x), o))
SymExpr(s::AbstractString) = new_expr(o -> mdac_expr_parse(s, o))

"""
    dasymbols(names)

The symbols named in `names`, separated by spaces or commas, as a tuple of `SymExpr`:
`a, b = dasymbols("a b")`.
"""
dasymbols(names::AbstractString) =
    Tuple(new_expr(o -> mdac_expr_symbol(n, o)) for n in split(names, (' ', ',', '\t', '\n'); keepempty=false))

Base.copy(x::SymExpr) = new_expr(o -> mdac_expr_copy(x, o))
Base.deepcopy_internal(x::SymExpr, d::IdDict) = get!(() -> copy(x), d, x)::SymExpr

function to_string(x::SymExpr)
    n = Ref{Csize_t}(0)
    check(mdac_expr_to_string(x, C_NULL, 0, n))
    buf = Vector{UInt8}(undef, n[])
    check(mdac_expr_to_string(x, buf, length(buf), n))
    return String(resize!(buf, length(buf) - 1))
end

Base.print(io::IO, x::SymExpr) = print(io, to_string(x))
Base.show(io::IO, x::SymExpr) = print(io, "SymExpr(", repr(to_string(x)), ")")

"""`Float64(x::SymExpr)`: the numeric value; `ArgumentError` if free symbols remain."""
function Base.Float64(x::SymExpr)
    r = Ref{Cdouble}()
    check(mdac_expr_to_double(x, r))
    return r[]
end

function Base.:(==)(a::SymExpr, b::SymExpr)
    r = Ref{Cint}()
    check(mdac_expr_eq(a, b, r))
    return r[] != 0
end
Base.:(==)(a::SymExpr, b::Real) = a == SymExpr(b)
Base.:(==)(a::Real, b::SymExpr) = SymExpr(a) == b

function Base.hash(x::SymExpr, h::UInt)
    r = Ref{UInt64}()
    check(mdac_expr_hash(x, r))
    return hash(r[], hash(SymExpr, h))
end

"""`iszero(x::SymExpr)`: true if `x` is zero (the test `SDA` uses: simplify and expand)."""
function Base.iszero(x::SymExpr)
    r = Ref{Cint}()
    check(mdac_expr_is_zero(x, r))
    return r[] != 0
end

for (op, name) in ((:+, :add), (:-, :sub), (:*, :mul), (:/, :div), (:^, :pow))
    f, fd, df = Symbol(:mdac_expr_, name), Symbol(:mdac_expr_, name, :_d), Symbol(:mdac_expr_d, name)
    @eval begin
        Base.$op(a::SymExpr, b::SymExpr) = new_expr(o -> $f(a, b, o))
        Base.$op(a::SymExpr, x::Integer) = $op(a, SymExpr(x))
        Base.$op(x::Integer, a::SymExpr) = $op(SymExpr(x), a)
        Base.$op(a::SymExpr, x::Real) = new_expr(o -> $fd(a, x, o))
        Base.$op(x::Real, a::SymExpr) = new_expr(o -> $df(x, a, o))
    end
end

Base.:-(a::SymExpr) = new_expr(o -> mdac_expr_neg(a, o))
Base.:+(a::SymExpr) = a

# Keys and values as SymExpr handles, for the C substitution functions.
function subs_args(d::AbstractDict)
    ks, vs = SymExpr[SymExpr(k) for k in keys(d)], SymExpr[SymExpr(v) for v in values(d)]
    return ks, vs, Handle[k.ptr for k in ks], Handle[v.ptr for v in vs]
end

"""
    MiraDAC.subs(x, d::AbstractDict)

`x` (a `SymExpr` or an `SDA`, in every coefficient) with `d`'s values substituted for its keys
(symbols); keys and values are `SymExpr`s or numbers. Not exported (SymEngine.jl exports
`subs`).
"""
function subs(x::SymExpr, d::AbstractDict)
    ks, vs, kp, vp = subs_args(d)
    GC.@preserve ks vs new_expr(o -> mdac_expr_subs(x, length(kp), kp, vp, o))
end

"""`MiraDAC.expand(x)`: `x` (a `SymExpr` or every coefficient of an `SDA`) expanded."""
expand(x::SymExpr) = new_expr(o -> mdac_expr_expand(x, o))

"""`simplify(x)`: `x` (a `SymExpr` or every coefficient of an `SDA`) simplified."""
simplify(x::SymExpr) = new_expr(o -> mdac_expr_simplify(x, o))

"""`diff(x::SymExpr, s::SymExpr)`: the derivative with respect to the symbol `s`."""
Base.diff(x::SymExpr, s::SymExpr) = new_expr(o -> mdac_expr_diff(x, s, o))

"""`MiraDAC.free_symbols(x::SymExpr)`: the set of symbols in `x`."""
function free_symbols(x::SymExpr)
    n = Ref{Csize_t}(0)
    check(mdac_expr_free_symbols(x, C_NULL, 0, n))
    buf = fill(C_NULL, n[])
    check(mdac_expr_free_symbols(x, buf, length(buf), n))
    return Set{SymExpr}(SymExpr(p) for p in buf)
end

# ---- SDA ------------------------------------------------------------------------------------

"""
    SDA()
    SDA(x)

A symbolic DA vector of the current env: a truncated power series whose coefficients are
`SymExpr`s. `SDA()` is zero; `SDA(x)` the constant `x` (a `SymExpr`, an `Integer` as an exact
integer, or another `Real`). Not a `Number`. Operators combine it with `SDA`, `NDA`,
`SymExpr` and `Real` on either side (an `NDA` with a `SymExpr` gives an `SDA` too), with
in-place forms (`add!`, `exp!`, ...). Its math functions: `sqrt exp log sin cos tan asin acos
atan sinh cosh tanh erf` and `^` (C++ has no SDA `asinh acosh atanh abs`).
"""
mutable struct SDA
    ptr::Handle
    function SDA(p::Handle)
        v = new(p)
        return adopt!(v)     # a finalizer, or the current dascope
    end
end

sda_finalizer(v::SDA) = defer_free(SDA_QUEUE, v)

Base.unsafe_convert(::Type{Handle}, v::SDA) = handle(v)
op_env(v::SDA) = mdac_sda_env(v)

# A new SDA from the allocating C call f(out) in the env of x; d as in alloc_call.
@inline function new_sda(f, x, d=nothing)
    drain!(SDA_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x, d)
    return SDA(out[])
end

# An in-place C call f() into out; it takes temporary slots, so it needs the retry.
@inline function sda_into(f, out, d=nothing)
    drain!(SDA_QUEUE, out)
    alloc_call(f, out, d)
    return out
end

function SDA(; env::Union{DAEnv,Nothing}=nothing)
    e = env_handle(env)
    return new_sda(o -> mdac_sda_new(e, C_NULL, o), e)
end

function SDA(x::SymExpr; env::Union{DAEnv,Nothing}=nothing)
    e = env_handle(env)
    return new_sda(o -> mdac_sda_new(e, x, o), e)
end

SDA(x::Integer; env::Union{DAEnv,Nothing}=nothing) = SDA(SymExpr(x); env)

function SDA(x::Real; env::Union{DAEnv,Nothing}=nothing)
    e = env_handle(env)
    return new_sda(o -> mdac_sda_new_d(e, x, o), e)
end

"""
    sdavar(i; env=current_env())

The symbolic base vector of variable `i` (1-based) of `env`: `promote_sda(davar(i))`.
"""
function sdavar(i::Integer; env::Union{DAEnv,Nothing}=nothing)
    e = env === nothing ? current_env() : env
    n = e.nvars
    1 <= i <= n || throw(BoundsError(1:n, i))
    h = op_env(e)
    return new_sda(o -> mdac_sda_var(h, i - 1, o), h)
end

"""
    promote_sda(v::NDA)

`v` as an `SDA` of its env; integer-valued coefficients become exact integers.
"""
promote_sda(v::NDA) = new_sda(o -> mdac_sda_promote(v, o), v)

"""
    promote_sda(env, v)

`promote_sda(v)` (`NDA` to `SDA`, `CNDA` to `CSDA`) with the result in `env`; the current env
does not change. `ArgumentError` if the two envs have different variables or orders.
"""
promote_sda(e::DAEnv, v::NDA) = new_sda(o -> mdac_nda_promote_to(e, v, o), e)

Base.copy(v::SDA) = new_sda(o -> mdac_sda_copy(v, o), v)
Base.deepcopy_internal(v::SDA, d::IdDict) = get!(() -> copy(v), d, v)::SDA

env(v::SDA) = DAEnv(mdac_sda_env(v))
import_vec(e::DAEnv, v::SDA) = new_sda(o -> mdac_sda_import(e, v, o), e)

"""The constant term of `v`, a `SymExpr`."""
con(v::SDA) = new_expr(o -> mdac_sda_con(v, o))

"""The stored coefficients of `v` in monomial order, as `SymExpr`s."""
function coeffs(v::SDA)
    n = Ref{Csize_t}(0)
    check(mdac_sda_coeffs(v, C_NULL, 0, n))
    drain!(EXPR_QUEUE, C_NULL)
    buf = fill(C_NULL, n[])
    check(mdac_sda_coeffs(v, buf, length(buf), n))
    return SymExpr[SymExpr(p) for p in buf]
end

function coeff(v::SDA, exps)
    e = collect(Cint, exps)
    return new_expr(o -> mdac_sda_coeff(v, e, length(e), o))
end

function nterms(v::SDA)
    r = Ref{Csize_t}()
    check(mdac_sda_nterms(v, r))
    return Int(r[])
end

"""`iszero(v::SDA)`: true if every coefficient is zero (`iszero(::SymExpr)`)."""
function Base.iszero(v::SDA)
    r = Ref{Cint}()
    check(mdac_sda_iszero(v, r))
    return r[] != 0
end

function to_string(v::SDA)
    n = Ref{Csize_t}(0)
    check(mdac_sda_to_string(v, C_NULL, 0, n))
    buf = Vector{UInt8}(undef, n[])
    check(mdac_sda_to_string(v, buf, length(buf), n))
    return String(resize!(buf, length(buf) - 1))
end

function Base.show(io::IO, v::SDA)
    getfield(v, :ptr) == C_NULL && return print(io, "SDA(freed)")   # by its dascope
    e = env(v)
    if e.retired
        print(io, "SDA(cleared env)")
    else
        print(io, "SDA(order=", e.order, ", nvars=", e.nvars, ", nonzero=", nterms(v), ")")
    end
end

function Base.show(io::IO, m::MIME"text/plain", v::SDA)
    getfield(v, :ptr) == C_NULL && return print(io, "SDA(freed)")   # by its dascope
    env(v).retired ? show(io, v) : print(io, to_string(v))
end

# ---- SDA arithmetic -------------------------------------------------------------------------

# Every operand pair of miradac.h, in the order of sda_op_shapes.
const SDA_OP_ARGS = ((:SDA, :SDA), (:SDA, :NDA), (:NDA, :SDA), (:SDA, :SymExpr), (:SymExpr, :SDA),
                     (:SDA, :Real), (:Real, :SDA), (:NDA, :SymExpr), (:SymExpr, :NDA))

for (op, name, bang) in ((:+, :add, :add!), (:-, :sub, :sub!), (:*, :mul, :mul!), (:/, :div, :div!))
    target = bang === :mul! ? :(LinearAlgebra.mul!) : bang
    for ((f, _), (A, B)) in zip(sda_op_shapes(name), SDA_OP_ARGS)
        fi = Symbol(f, :_into)
        da = A === :Real || A === :SymExpr ? :b : :a     # the DA operand gives the env
        @eval begin
            Base.$op(a::$A, b::$B) = new_sda(o -> $f(a, b, o), $da)
            $target(out::SDA, a::$A, b::$B) = sda_into(() -> $fi(out, a, b), out)
        end
    end
    # An Integer is an exact SymExpr integer, as in SDA(x).
    @eval begin
        Base.$op(a::SDA, x::Integer) = $op(a, SymExpr(x))
        Base.$op(x::Integer, a::SDA) = $op(SymExpr(x), a)
        $target(out::SDA, a::SDA, x::Integer) = $target(out, a, SymExpr(x))
        $target(out::SDA, x::Integer, a::SDA) = $target(out, SymExpr(x), a)
    end
end

Base.:-(a::SDA) = new_sda(o -> mdac_sda_neg(a, o), a)
Base.:^(a::SDA, n::Integer) = new_sda(o -> mdac_sda_pow_i(a, n, o), a, a)
Base.:^(a::SDA, x::Real) = new_sda(o -> mdac_sda_pow_d(a, x, o), a, a)

for fn in SDA_FUNCS
    f, fi, bang = Symbol(:mdac_sda_, fn), Symbol(:mdac_sda_, fn, :_into), Symbol(fn, :!)
    target = fn === :erf ? fn : :(Base.$fn)
    @eval begin
        $target(a::SDA) = new_sda(o -> $f(a, o), a, a)
        $bang(out::SDA, a::SDA) = sda_into(() -> $fi(out, a), out, a)
    end
end

# ---- SDA per-coefficient operations and evaluation -------------------------------------------

simplify(v::SDA) = new_sda(o -> mdac_sda_simplify(v, o), v)
expand(v::SDA) = new_sda(o -> mdac_sda_expand(v, o), v)

function subs(v::SDA, d::AbstractDict)
    ks, vs, kp, vp = subs_args(d)
    GC.@preserve ks vs new_sda(o -> mdac_sda_subs(v, length(kp), kp, vp, o), v)
end

"""
    evaluate(s::SDA, d::AbstractDict)

The `NDA` of `s` with the numbers `d[k]` for the symbols `k` (`SymExpr`s); every symbol of `s`
needs a value (`ArgumentError` otherwise).
"""
function evaluate(s::SDA, d::AbstractDict{<:Any,<:Real})
    ks = SymExpr[SymExpr(k) for k in keys(d)]
    kp, vals = Handle[k.ptr for k in ks], Float64[v for v in values(d)]
    GC.@preserve ks new_nda(o -> mdac_sda_evaluate(s, length(kp), kp, vals, o), s)
end

# ---- SDAList and algorithms -----------------------------------------------------------------

"""
    SDAList()
    SDAList(v::AbstractVector{SDA})

A list of `SDA` held in C++; it behaves as `NDAList`.
"""
mutable struct SDAList <: AbstractVector{SDA}
    ptr::Handle
    function SDAList(p::Handle)
        l = new(p)
        return adopt!(l)     # a finalizer, or the current dascope
    end
end

slist_finalizer(l::SDAList) = defer_free(SLIST_QUEUE, l)

Base.unsafe_convert(::Type{Handle}, l::SDAList) = handle(l)
op_env(l::SDAList) = mdac_sdalist_env(l)

function new_slist(f, x)
    drain!(SLIST_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x)
    return SDAList(out[])
end

function SDAList()
    out = Ref{Handle}(C_NULL)
    check(mdac_sdalist_new(out))
    return SDAList(out[])
end

function SDAList(v::AbstractVector{SDA})
    isempty(v) && return SDAList()
    ptrs = Handle[handle(x) for x in v]
    GC.@preserve v new_slist(o -> mdac_sdalist_from(ptrs, length(ptrs), o), first(v))
end

aslist(l::SDAList) = l
aslist(v::AbstractVector{SDA}) = SDAList(v)

Base.size(l::SDAList) = (Int(mdac_sdalist_length(l)),)

function Base.getindex(l::SDAList, i::Int)
    @boundscheck checkbounds(l, i)
    return new_sda(o -> mdac_sdalist_get(l, i - 1, o), l)
end

function Base.setindex!(l::SDAList, v::SDA, i::Int)
    @boundscheck checkbounds(l, i)
    check(mdac_sdalist_set(l, i - 1, v))
    return l
end

function Base.push!(l::SDAList, v::SDA)
    alloc_call(() -> mdac_sdalist_push(l, v), v)
    return l
end

function env(l::SDAList)
    e = mdac_sdalist_env(l)
    e == C_NULL && throw(EnvError("an empty SDAList has no env"))
    return DAEnv(e)
end

der(v::SDA, i::Integer) = new_sda(o -> mdac_sda_der(v, cvar(v, i), o), v)
integ(v::SDA, i::Integer) = new_sda(o -> mdac_sda_integ(v, cvar(v, i), o), v)

function substitute(v::SDA, i::Integer, x::Real)
    c = cvar(v, i)
    return new_sda(o -> mdac_sda_substitute_d(v, c, x, o), v)
end

function substitute(v::SDA, i::Integer, x::SDA)
    c = cvar(v, i)
    return new_sda(o -> mdac_sda_substitute(v, c, x, o), v)
end

function substitute(v::SDA, ids::AbstractVector{<:Integer}, xs::AbstractVector{SDA})
    c, l = Cuint[cvar(v, i) for i in ids], aslist(xs)
    return new_sda(o -> mdac_sda_substitute_multi(v, c, length(c), l, o), v)
end

function substitute(m::AbstractVector{SDA}, ids::AbstractVector{<:Integer}, xs::AbstractVector{SDA})
    isempty(m) && return SDAList()
    ml, l = aslist(m), aslist(xs)
    c = Cuint[cvar(ml, i) for i in ids]
    return new_slist(o -> mdac_sdalist_substitute(ml, c, length(c), l, o), ml)
end

function compose(m::AbstractVector{SDA}, v::AbstractVector{SDA})
    isempty(m) && return SDAList()
    ml, vl = aslist(m), aslist(v)
    return new_slist(o -> mdac_sdalist_compose(ml, vl, o), ml)
end

# At a point: the values, as SymExprs.
function compose(m::AbstractVector{SDA}, pt::AbstractVector{<:Real})
    ml, p = aslist(m), convert(Vector{Float64}, pt)
    isempty(ml) && return SymExpr[]
    drain!(EXPR_QUEUE, C_NULL)
    drain!(SDA_QUEUE, ml)
    buf = fill(C_NULL, length(ml))
    # The kernel takes nvars + length(m) temporary slots, so it needs the retry.
    alloc_call(() -> mdac_sdalist_compose_d(ml, p, length(p), buf), ml)
    return SymExpr[SymExpr(q) for q in buf]
end
