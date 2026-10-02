# scope.jl — explicit scopes that free temporaries promptly (DEVELOPMENT_PLAN_JULIA.md A.5b).
#
# An object created while a scope is active on the current task gets no finalizer: it is
# registered in the innermost scope, and freed through the C API when the scope ends unless it
# is kept (reachable from the block's value, or passed to keep!). A kept object moves to the
# enclosing scope, or at the outermost scope gets its deferred-free finalizer. Freeing at once,
# in creation order, lets the next result reuse the warm slot just freed (the C++ pool hands
# out the most recently freed slot first).

"""
    FreedObjectError

A MiraDAC object was used after the `dascope` that created it ended and freed it. Return the
object from the `dascope` block, or call `keep!(x)` inside the block, to use it afterwards.
"""
struct FreedObjectError <: MiraDACError
    msg::String
end

@noinline throw_freed(x) = throw(FreedObjectError(
    "this $(nameof(typeof(x))) was freed when the dascope that created it ended; return it " *
    "from the dascope block or call keep!(x) inside the block to use it afterwards"))

# The checked handle of every object type: one pointer comparison.
@inline handle(x) = (p = getfield(x, :ptr); p == C_NULL ? throw_freed(x) : p)

const DAObject = Union{NDA,CNDA,NDAList,CNDAList}

mutable struct Scope
    const ndas::Vector{NDA}
    const cndas::Vector{CNDA}
    const lists::Vector{NDAList}
    const clists::Vector{CNDAList}
    kept::Union{Nothing,Base.IdSet{Any}}     # created on first use (Base.IdSet: not exported before 1.11)
end

Scope() = Scope(NDA[], CNDA[], NDAList[], CNDAList[], nothing)

store(s::Scope, ::NDA) = s.ndas
store(s::Scope, ::CNDA) = s.cndas
store(s::Scope, ::NDAList) = s.lists
store(s::Scope, ::CNDAList) = s.clists

# Per task: scopes[1:depth] are active; deeper ones are kept for reuse (no allocation per scope).
mutable struct ScopeStack
    const scopes::Vector{Scope}
    depth::Int
end

const SCOPE_KEY = :MiraDAC_scopes

# Active scopes over all tasks: while 0, creating an object skips the task-local lookup.
const ACTIVE_SCOPES = Threads.Atomic{Int}(0)

@inline current_scope() = ACTIVE_SCOPES[] == 0 ? nothing : task_scope()

function task_scope()
    st = get(task_local_storage(), SCOPE_KEY, nothing)::Union{Nothing,ScopeStack}
    return st === nothing || st.depth == 0 ? nothing : @inbounds st.scopes[st.depth]
end

attach_finalizer(v::NDA) = finalizer(nda_finalizer, v)
attach_finalizer(v::CNDA) = finalizer(cnda_finalizer, v)
attach_finalizer(l::NDAList) = finalizer(list_finalizer, l)
attach_finalizer(l::CNDAList) = finalizer(clist_finalizer, l)

free_now(v::NDA) = mdac_nda_free(getfield(v, :ptr))
free_now(v::CNDA) = mdac_cnda_free(getfield(v, :ptr))
free_now(l::NDAList) = mdac_ndalist_free(getfield(l, :ptr))
free_now(l::CNDAList) = mdac_cndalist_free(getfield(l, :ptr))

# Called by every inner constructor: register in the current scope, or attach the finalizer.
@inline function adopt!(v)
    s = current_scope()
    s === nothing ? attach_finalizer(v) : push!(store(s, v), v)
    return v
end

"""
    dascope(f)
    dascope() do ... end

Run `f()` and return its value. MiraDAC objects (`NDA`, `CNDA`, `NDAList`, `CNDAList`)
created meanwhile on this task are freed as soon as `f` returns or throws, instead of after
a garbage collection, except those reachable from the returned value (directly or inside
`Tuple`, `NamedTuple`, `AbstractArray` and `AbstractDict` values; a lazy wrapper such as a
`view` or `v'` keeps the array it wraps) and those passed to `keep!`. Kept objects belong to the enclosing `dascope`, or, at the outermost one, are left to
the garbage collector as usual. Using a freed object throws `FreedObjectError`.

Use it around a step of a computation that creates many temporaries: they cost less (no
finalizer), and each new result reuses the cache-warm pool slot of a freed one, so the pool
never fills up with garbage. Scopes nest; `with_env` and `with_order` may be used inside.
In-place operations (`add!`, `exp!`, ...) create no objects and need no scope.

```julia
y = dascope() do
    t = exp(x) * sin(x)       # temporaries: freed at the end of the block
    t + 1.0                   # returned: kept
end
```
"""
function dascope(f)
    st = get!(ScopeStack, task_local_storage(), SCOPE_KEY)::ScopeStack
    d = st.depth + 1
    d > length(st.scopes) && push!(st.scopes, Scope())
    s = @inbounds st.scopes[d]
    st.depth = d
    Threads.atomic_add!(ACTIVE_SCOPES, 1)
    try
        r = f()
        keep!(s, r)
        return r
    finally
        st.depth = d - 1
        Threads.atomic_sub!(ACTIVE_SCOPES, 1)
        close!(s, d > 1 ? @inbounds(st.scopes[d-1]) : nothing)
    end
end

ScopeStack() = ScopeStack(Scope[], 0)

"""
    keep!(x)

Keep `x` (a MiraDAC object, or a `Tuple`, `NamedTuple`, `AbstractArray` or `AbstractDict`
holding some) alive when the innermost `dascope` ends: it then moves to the enclosing scope,
or is left to the garbage collector. Returns `x`; outside a `dascope` it does nothing.
"""
function keep!(x)
    s = current_scope()
    s === nothing || keep!(s, x)
    return x
end

kept_set(s::Scope) = (k = s.kept; k === nothing ? (s.kept = Base.IdSet{Any}()) : k)

keep!(s::Scope, x) = nothing
keep!(s::Scope, x::DAObject) = (push!(kept_set(s), x); nothing)
# Lists are arrays too; this is more specific than both.
keep!(s::Scope, x::Union{NDAList,CNDAList}) = (push!(kept_set(s), x); nothing)
keep!(s::Scope, x::Union{Tuple,NamedTuple}) = foreach(y -> keep!(s, y), x)

function keep!(s::Scope, x::Union{AbstractArray,AbstractDict})
    x isa AbstractArray && isbitstype(eltype(x)) && return nothing
    k = kept_set(s)
    x in k && return nothing           # visited container (cycles)
    push!(k, x)
    if x isa AbstractDict
        foreach(y -> keep!(s, y), values(x))
    elseif !(x isa Array) && (p = parent(x)) !== x
        # A lazy wrapper (view, adjoint, reshape, ...): keep the array it wraps. Indexing the
        # wrapper could create objects (a view of an NDAList copies) or fail (adjoint of NDA).
        keep!(s, p)
    else
        for i in eachindex(x)          # skip #undef slots of a partly filled array
            isassigned(x, i) && keep!(s, @inbounds x[i])
        end
    end
end

# End of scope s: free what is not kept; move the rest to `parent` or give it a finalizer.
function close!(s::Scope, parent)
    k = s.kept
    s.kept = nothing
    close!(s.ndas, k, parent)
    close!(s.cndas, k, parent)
    close!(s.lists, k, parent)
    close!(s.clists, k, parent)
end

function close!(vs::Vector, k, parent)
    for v in vs
        if k !== nothing && v in k
            parent === nothing ? attach_finalizer(v) : push!(store(parent, v), v)
        else
            free_now(v)
            setfield!(v, :ptr, C_NULL)
        end
    end
    empty!(vs)
    return nothing
end
