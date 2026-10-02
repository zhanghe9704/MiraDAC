# nda.jl — numeric DA vectors (plan T2.2, T2.3).

"""
    NDA(x::Real)
    NDA(coeffs::AbstractVector{<:Real})

A numeric DA vector of the current env: the constant `x`, or the coefficients `coeffs` in
monomial order (`coeffs(v)` gives them back). Not a `Number`; the operators and math functions
are defined explicitly, with in-place forms (`add!`, `exp!`, ...) for hot loops.
"""
mutable struct NDA
    ptr::Handle
    function NDA(p::Handle)
        return adopt!(new(p))     # a finalizer, or the current dascope
    end
end

nda_finalizer(v::NDA) = defer_free(NDA_QUEUE, v)

Base.unsafe_convert(::Type{Handle}, v::NDA) = handle(v)
op_env(v::NDA) = mdac_nda_env(v)
op_env(e::DAEnv) = getfield(e, :ptr)

# A new NDA from the allocating C call f(out) in the env of x; d as in alloc_call.
@inline function new_nda(f, x, d=nothing)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x, d)
    return NDA(out[])
end

function NDA(x::Real)
    e = current_env_handle()
    return new_nda(o -> mdac_nda_new(e, x, o), e)
end

function NDA(c::AbstractVector{<:Real})
    e = current_env_handle()
    c64 = convert(Vector{Float64}, c)
    return new_nda(o -> mdac_nda_from_coeffs(e, c64, length(c64), o), e)
end

"""
    davar(i)

The `i`-th variable (1-based) of the current env: the vector with coefficient 1 for `x_i`.
"""
function davar(i::Integer)
    e = current_env()
    n = e.nvars
    1 <= i <= n || throw(BoundsError(1:n, i))
    h = op_env(e)
    return new_nda(o -> mdac_nda_var(h, i - 1, o), h)
end

Base.copy(v::NDA) = new_nda(o -> mdac_nda_copy(v, o), v)
Base.deepcopy_internal(v::NDA, d::IdDict) = get!(() -> copy(v), d, v)::NDA

"""The env of `v`."""
env(v::NDA) = DAEnv(mdac_nda_env(v))

# ---- Inspection -----------------------------------------------------------------------------

"""The constant term of `v`."""
function con(v::NDA)
    r = Ref{Cdouble}()
    check(mdac_nda_con(v, r))
    return r[]
end

"""The stored coefficients of `v` in monomial order (up to the last non-zero one)."""
function coeffs(v::NDA)
    n = Ref{Csize_t}(0)
    check(mdac_nda_coeffs(v, C_NULL, 0, n))
    buf = Vector{Float64}(undef, n[])
    check(mdac_nda_coeffs(v, buf, length(buf), n))
    return buf
end

"""
    coeff(v, exps)

The coefficient of the monomial with exponents `exps` (one per variable; missing trailing
ones are 0).
"""
function coeff(v::NDA, exps)
    e = collect(Cint, exps)
    r = Ref{Cdouble}()
    check(mdac_nda_coeff(v, e, length(e), r))
    return r[]
end

"""The number of non-zero coefficients of `v`."""
function nterms(v::NDA)
    r = Ref{Csize_t}()
    check(mdac_nda_nterms(v, r))
    return Int(r[])
end

"""The largest coefficient magnitude of `v`."""
function LinearAlgebra.norm(v::NDA)
    r = Ref{Cdouble}()
    check(mdac_nda_norm(v, r))
    return r[]
end

"""
    iszero(v::NDA; eps=get_eps())

True if every coefficient of `v` is below `eps` in magnitude.
"""
function Base.iszero(v::NDA; eps::Real=get_eps())
    r = Ref{Cint}()
    check(mdac_nda_iszero(v, eps, r))
    return r[] != 0
end

# The C++ operator<< text.
function to_string(v::NDA)
    n = Ref{Csize_t}(0)
    check(mdac_nda_to_string(v, C_NULL, 0, n))
    buf = Vector{UInt8}(undef, n[])
    check(mdac_nda_to_string(v, buf, length(buf), n))
    return String(resize!(buf, length(buf) - 1))
end

function Base.show(io::IO, v::NDA)
    e = env(v)
    if e.retired
        print(io, "NDA(cleared env)")
    else
        print(io, "NDA(order=", e.order, ", nvars=", e.nvars, ", nonzero=", nterms(v), ")")
    end
end

function Base.show(io::IO, m::MIME"text/plain", v::NDA)
    env(v).retired ? show(io, v) : print(io, to_string(v))
end

# ---- Arithmetic -----------------------------------------------------------------------------

for (op, name, bang) in ((:+, :add, :add!), (:-, :sub, :sub!), (:*, :mul, :mul!), (:/, :div, :div!))
    f, fd, df = Symbol(:mdac_nda_, name), Symbol(:mdac_nda_, name, :_d), Symbol(:mdac_nda_d, name)
    fi, fdi, dfi = Symbol(f, :_into), Symbol(fd, :_into), Symbol(df, :_into)
    target = bang === :mul! ? :(LinearAlgebra.mul!) : bang
    @eval begin
        Base.$op(a::NDA, b::NDA) = new_nda(o -> $f(a, b, o), a)
        Base.$op(a::NDA, x::Real) = new_nda(o -> $fd(a, x, o), a)
        Base.$op(x::Real, a::NDA) = new_nda(o -> $df(x, a, o), a)
        $target(out::NDA, a::NDA, b::NDA) = (check($fi(out, a, b)); out)
        $target(out::NDA, a::NDA, x::Real) = (check($fdi(out, a, x)); out)
        $target(out::NDA, x::Real, a::NDA) = (check($dfi(out, x, a)); out)
    end
end

@doc "`add!(out, a, b)`: `out = a + b` in `out`'s own storage (no allocation); `a` or `b` may be a `Real`." add!
@doc "`sub!(out, a, b)`: `out = a - b` in place; see `add!`." sub!
@doc "`div!(out, a, b)`: `out = a / b` in place; see `add!`." div!

Base.:-(a::NDA) = new_nda(o -> mdac_nda_neg(a, o), a)
Base.:^(a::NDA, n::Integer) = new_nda(o -> mdac_nda_pow_i(a, n, o), a, a)
Base.:^(a::NDA, x::Real) = new_nda(o -> mdac_nda_pow_d(a, x, o), a, a)

# ---- Math functions -------------------------------------------------------------------------

"""
    erf(v::NDA)

The error function of a DA vector (MiraDAC's own; no SpecialFunctions dependency).
"""
function erf end

for fn in NDA_FUNCS
    f, fi, bang = Symbol(:mdac_nda_, fn), Symbol(:mdac_nda_, fn, :_into), Symbol(fn, :!)
    target = fn === :erf ? fn : :(Base.$fn)
    @eval begin
        $target(a::NDA) = new_nda(o -> $f(a, o), a, a)
        @doc $("`$bang(out, a)`: `out = $fn(a)` in `out`'s own storage (no allocation).")
        $bang(out::NDA, a::NDA) = (check($fi(out, a), a); out)
    end
end

"""`abs(v::NDA)`: the largest coefficient magnitude, as in C++ (a `Float64`)."""
function Base.abs(a::NDA)
    r = Ref{Cdouble}()
    check(mdac_nda_abs(a, r))
    return r[]
end
