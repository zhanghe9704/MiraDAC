# algorithms.jl — NDAList and the NDA algorithms (plan T3.2).

"""
    NDAList()
    NDAList(v::AbstractVector{NDA})

A list of `NDA` held in C++, for the map-level functions (`compose`, `inv_map`,
`evaluate_map`, ...). Built from a vector, it holds a copy of each element; `l[i]` returns a
copy and `l[i] = v` copies `v` into the list. All elements belong to one env. Functions taking
a map accept a plain `Vector{NDA}` too, copied into a temporary `NDAList`.
"""
mutable struct NDAList <: AbstractVector{NDA}
    ptr::Handle
    function NDAList(p::Handle)
        return adopt!(new(p))     # a finalizer, or the current dascope
    end
end

list_finalizer(l::NDAList) = defer_free(LIST_QUEUE, l)

Base.unsafe_convert(::Type{Handle}, l::NDAList) = handle(l)
op_env(l::NDAList) = mdac_ndalist_env(l)

# A new NDAList from the allocating C call f(out) in the env of x.
function new_list(f, x)
    drain!(LIST_QUEUE, x)
    out = Ref{Handle}(C_NULL)
    alloc_call(() -> f(out), x)
    return NDAList(out[])
end

function NDAList()
    out = Ref{Handle}(C_NULL)
    check(mdac_ndalist_new(out))
    return NDAList(out[])
end

function NDAList(v::AbstractVector{NDA})
    isempty(v) && return NDAList()
    ptrs = Handle[handle(x) for x in v]
    GC.@preserve v new_list(o -> mdac_ndalist_from(ptrs, length(ptrs), o), first(v))
end

aslist(l::NDAList) = l
aslist(v::AbstractVector{NDA}) = NDAList(v)

Base.size(l::NDAList) = (Int(mdac_ndalist_length(l)),)

function Base.getindex(l::NDAList, i::Int)
    @boundscheck checkbounds(l, i)
    return new_nda(o -> mdac_ndalist_get(l, i - 1, o), l)
end

function Base.setindex!(l::NDAList, v::NDA, i::Int)
    @boundscheck checkbounds(l, i)
    check(mdac_ndalist_set(l, i - 1, v))
    return l
end

function Base.push!(l::NDAList, v::NDA)
    alloc_call(() -> mdac_ndalist_push(l, v), v)
    return l
end

"""The env of the list's vectors (`EnvError` for an empty list)."""
function env(l::NDAList)
    e = mdac_ndalist_env(l)
    e == C_NULL && throw(EnvError("an empty NDAList has no env"))
    return DAEnv(e)
end

# 1-based variable i of x's env as a 0-based C index.
function cvar(x, i::Integer)
    n = env(x).nvars
    1 <= i <= n || throw(BoundsError(1:n, i))
    return Cuint(i - 1)
end

"""
    der(v, i)

The derivative of `v` with respect to variable `i` (1-based).
"""
der(v::NDA, i::Integer) = new_nda(o -> mdac_nda_der(v, cvar(v, i), o), v)

"""
    integ(v, i)

The integral of `v` with respect to variable `i` (1-based).
"""
integ(v::NDA, i::Integer) = new_nda(o -> mdac_nda_integ(v, cvar(v, i), o), v)

"""
    substitute(v, i, x)
    substitute(v, ids, xs)

`v` with `x` (a number or an `NDA`) substituted for variable `i` (1-based), or with `xs[j]`
substituted for variable `ids[j]` for each `j`. `v` may be a map (an `NDAList` or a
`Vector{NDA}`); the result is then an `NDAList`.
"""
function substitute(v::NDA, i::Integer, x::Real)
    c = cvar(v, i)
    return new_nda(o -> mdac_nda_substitute_d(v, c, x, o), v)
end

function substitute(v::NDA, i::Integer, x::NDA)
    c = cvar(v, i)
    return new_nda(o -> mdac_nda_substitute(v, c, x, o), v)
end

function substitute(v::NDA, ids::AbstractVector{<:Integer}, xs::AbstractVector{NDA})
    c, l = Cuint[cvar(v, i) for i in ids], aslist(xs)
    return new_nda(o -> mdac_nda_substitute_multi(v, c, length(c), l, o), v)
end

function substitute(m::AbstractVector{NDA}, ids::AbstractVector{<:Integer}, xs::AbstractVector{NDA})
    isempty(m) && return NDAList()
    ml, l = aslist(m), aslist(xs)
    c = Cuint[cvar(ml, i) for i in ids]
    return new_list(o -> mdac_ndalist_substitute(ml, c, length(c), l, o), ml)
end

"""
    compose(m, v)

The map `m` (an `NDAList` or `Vector{NDA}`) composed with `v`: `m[i](v[1], …, v[nvars])`.
With DA vectors `v`, an `NDAList`; with a point `v` (real or complex numbers), the values as
a `Vector{Float64}` or `Vector{ComplexF64}`. When `m` or `v` holds `CNDA`s (C++
`cd_composition`; `m` and `v` may be `CNDAList`s), a `CNDAList`.
"""
function compose(m::AbstractVector{NDA}, v::AbstractVector{NDA})
    isempty(m) && return NDAList()
    ml, vl = aslist(m), aslist(v)
    return new_list(o -> mdac_ndalist_compose(ml, vl, o), ml)
end

function compose(m::AbstractVector{NDA}, pt::AbstractVector{<:Real})
    ml, p = aslist(m), convert(Vector{Float64}, pt)
    out = Vector{Float64}(undef, length(ml))
    isempty(ml) || check(mdac_ndalist_compose_d(ml, p, length(p), out))
    return out
end

function compose(m::AbstractVector{NDA}, pt::AbstractVector{<:Complex})
    ml, p = aslist(m), convert(Vector{ComplexF64}, pt)
    out = Vector{ComplexF64}(undef, length(ml))
    isempty(ml) || check(mdac_ndalist_compose_z(ml, p, length(p), out))
    return out
end

"""
    inv_map(m, dim=length(m))

The inverse of the map `m[1:dim]` (zero constant parts), an `NDAList` of `dim` vectors.
`ArgumentError` if the linear part of the map is singular.
"""
function inv_map(m::AbstractVector{NDA}, dim::Integer=length(m))
    isempty(m) && throw(ArgumentError("inv_map: empty map"))
    ml = aslist(m)
    return new_list(o -> mdac_ndalist_inv_map(ml, dim, o), ml)
end

"""
    evaluate_map(m, pts)

The map `m` at many points in one C call. `pts` is an `nvars × N` matrix, one point per
column; the result is the `length(m) × N` matrix whose column `p` is `compose(m, pts[:, p])`.
"""
function evaluate_map(m::AbstractVector{NDA}, pts::AbstractMatrix{<:Real})
    ml, n = aslist(m), size(pts, 2)
    out = Matrix{Float64}(undef, length(ml), n)
    isempty(ml) && return out
    nv = env(ml).nvars
    size(pts, 1) == nv || throw(DimensionMismatch("pts must have nvars = $nv rows"))
    p = convert(Matrix{Float64}, pts)
    check(mdac_ndalist_evaluate_map(ml, p, n, out))
    return out
end

"""
    (v::NDA)(pt)

The value of `v` at the point `pt` (`nvars` numbers).
"""
function (v::NDA)(pt::AbstractVector{<:Real})
    p = convert(Vector{Float64}, pt)
    r = Ref{Cdouble}()
    check(mdac_nda_eval(v, p, length(p), r))
    return r[]
end

"""
    exponents(env=current_env())

The exponents of every monomial, an `nvars × full_length` `Matrix{Int32}`: column `i` belongs
to coefficient `i` of `coeffs`.
"""
function exponents(e::DAEnv=current_env())
    nv, len = e.nvars, e.full_length
    out = Matrix{Int32}(undef, nv, len)
    n = Ref{Csize_t}(0)
    check(mdac_exponents(e, out, length(out), n))
    return out
end
