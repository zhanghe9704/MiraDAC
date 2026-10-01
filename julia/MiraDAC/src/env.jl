# env.jl — DA environments, the default env and the truncation order (plan T2.4).

"""
    DAEnv

A DA environment: variables, order and a pool of vectors. Envs belong to the library; a
`DAEnv` is a handle. Properties: `order` (the current truncation order), `max_order` (the
order at creation), `nvars`, `full_length`, `poolsize`, `count` (vectors in use), `remain`,
`retired`.
"""
struct DAEnv
    ptr::Handle
end

Base.unsafe_convert(::Type{Handle}, e::DAEnv) = getfield(e, :ptr)

function Base.getproperty(e::DAEnv, s::Symbol)
    s === :order && return env_query(mdac_env_order, e)
    s === :max_order && return env_query(mdac_env_max_order, e)
    s === :nvars && return env_query(mdac_env_nvars, e)
    s === :full_length && return env_query(mdac_env_full_length, e)
    s === :poolsize && return env_query(mdac_env_poolsize, e)
    s === :count && return env_query(mdac_env_count, e)
    s === :remain && return env_query(mdac_env_remain, e)
    s === :retired && return mdac_env_retired(e) != 0
    return getfield(e, s)
end

Base.propertynames(::DAEnv) =
    (:order, :max_order, :nvars, :full_length, :poolsize, :count, :remain, :retired)

function env_query(f, e::DAEnv)
    r = Ref{Cuint}()
    check(f(e, r))
    return Int(r[])
end

function Base.show(io::IO, e::DAEnv)
    if e.retired
        print(io, "DAEnv(retired)")
    else
        print(io, "DAEnv(order=", e.max_order, ", nvars=", e.nvars, ", poolsize=", e.poolsize, ")")
    end
end

"""
    init!(order, nvars, poolsize; table=false)

Create the default env (C++ `da_init`) and make it current. Vectors have `nvars` variables,
are truncated at `order`, and at most `poolsize` of them exist at a time.
"""
function init!(order::Integer, nvars::Integer, poolsize::Integer; table::Bool=false)
    check(mdac_init(order, nvars, poolsize, table))
    return default_env()
end

"""
    clear!()

Release the default env (C++ `da_clear`); its vectors then throw `EnvError` when used.
"""
clear!() = check(mdac_clear())

function current_env_handle()
    r = Ref{Handle}(C_NULL)
    check(mdac_env_current(r))
    return r[]
end

"""The env that new vectors belong to (`EnvError` if there is none)."""
current_env() = DAEnv(current_env_handle())

"""The env of the last `init!` (`EnvError` after `clear!`)."""
function default_env()
    r = Ref{Handle}(C_NULL)
    check(mdac_env_default(r))
    return DAEnv(r[])
end

"""
    with_order(f, n)

Run `f()` with the current env truncating at order `n` (at most its `max_order`); the order
before the call is restored afterwards, also when `f` throws. Nests correctly.
"""
function with_order(f, n::Integer)
    e = current_env()
    saved = e.order
    change_order(e, n)
    try
        return f()
    finally
        saved == e.max_order ? check(mdac_env_restore_order(e)) : change_order(e, saved)
    end
end

function change_order(e::DAEnv, n::Integer)
    n >= 0 || throw(ArgumentError("order must be non-negative, got $n"))
    ok = Ref{Cint}(0)
    check(mdac_env_change_order(e, n, ok))
    ok[] == 0 && throw(ArgumentError("order $n exceeds the env's order $(e.max_order)"))
    return nothing
end

"""Coefficients below `get_eps()` in magnitude are dropped (shared by all envs)."""
get_eps() = mdac_get_eps()

"""Set the threshold of `get_eps()`; `x` must be positive."""
set_eps!(x::Real) = check(mdac_set_eps(x))
