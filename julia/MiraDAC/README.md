# MiraDAC.jl

Julia binding of the MiraDAC differential algebra library, through its C API (`libmiradac_c`).
Development setup: build with `-DDA_BUILD_CAPI=ON`, then `julia julia/dev_setup.jl`. Design and
API: `DEVELOPMENT_PLAN_JULIA.md` (Part A).

```julia
using MiraDAC
init!(10, 6, 2000)             # order, nvars, pool size
x = davar(1) + 0.5
y = exp(x) * sin(x)
```

## Memory

Every DA vector takes one slot of its env's fixed-size pool. Julia frees an object only after a
garbage collection, so by default a vector is returned to the pool late: when the pool is full,
MiraDAC runs a collection and retries. This works, but each new object pays for a finalizer,
and a new result never reuses the cache-warm slot of a temporary that just died. Two tools
avoid that cost:

- **In-place operations** (`add!(out, a, b)`, `mul!`, `exp!(out, a)`, ...) write into an
  existing vector and create nothing. Use them in the innermost loops.
- **Scopes** (below) free the temporaries of a block as soon as it ends.

## Scopes

```julia
y = dascope() do
    t = exp(x) * sin(x)        # temporaries: freed when the block ends
    t + 1.0                    # the block's value: kept
end
```

`dascope(f)` runs `f()` and returns its value. Every MiraDAC object (`NDA`, `CNDA`, `NDAList`,
`CNDAList`) created inside the block on the current task is freed when the block ends, normally
or by an exception, except:

- objects reachable from the returned value (the value itself, or inside a `Tuple`,
  `NamedTuple`, `AbstractArray` or `AbstractDict`, recursively; a `view`, `v'` or `reshape`
  keeps the array it wraps), and
- objects passed to `keep!(x)` inside the block (use it for objects stored elsewhere, for
  example in a field of your own struct or a global).

A kept object belongs to the enclosing `dascope` if there is one; otherwise it is left to the
garbage collector as usual. Using a freed object throws `FreedObjectError`; nothing crashes.

When to use one: around a step of a computation that creates many temporaries — one
integration step, one function evaluation, a batch of a loop — especially when vectors are
large (a slot is 64 KB at 6 variables, order 10) or the pool is small. Inside a scope, objects
cost less (no finalizer), the next result reuses the warm slot of a freed one, and the pool never
fills with garbage, so no garbage collection is needed. Keep a scope's lifetime short: everything
created in it stays allocated until it ends, so a scope around a loop of a million operations
needs a pool of a million slots — put the scope inside the loop instead, and carry results across
iterations in place:

```julia
x = NDA(0.0)
for k in 1:100_000
    dascope() do
        add!(x, x * 0.5, davar(1))     # temporaries freed at once; result written into x
    end
end
```

Outside any scope nothing changes. Scopes nest and are per task (a task spawned inside a scope
is not in it).
