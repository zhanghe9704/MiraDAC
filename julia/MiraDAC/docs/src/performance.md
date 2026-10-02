# Performance

A call into the library costs a few tens of nanoseconds; most of the cost of DA arithmetic in
Julia comes from creating and freeing vectors, and the in-place API avoids both.

## The in-place API

`add!(out, a, b)`, `sub!`, `mul!` (a method of `LinearAlgebra.mul!`) and `div!` write `a ⊕ b`
into `out`, and `exp!(out, a)`, `sin!(out, a)`, … (one for every math function) write `f(a)`
into `out`. `out` keeps its pool slot, so an `NDA` loop allocates nothing:

```julia
init!(4, 3, 10_000)
x = davar(1); y = NDA(0.0); t = NDA(0.0)
for k in 1:10_000
    mul!(t, x, 0.5)                 # t = 0.5x
    exp!(y, t)                      # y = exp(t)
    add!(x, x, y)                   # x += y
end
```

`a` or `b` may also be a number, and `out` may be one of the operands. The CNDA, SDA and CSDA
forms compute into temporary slots and then copy, so they may also need free pool slots. Measured
against C++ (`bench/REPORT.md`), in-place operations below 1 µs add at most a few tens of
nanoseconds.

## Pool and garbage collector

Every new DA vector takes a slot in its env's fixed-size pool. Julia frees memory by garbage
collection, not by reference counting, so a vector that is no longer used keeps its slot until
the next collection runs its finalizer. The finalizer does not call the library; it puts the
vector on a free queue, which the next operation in that env drains.

When an operation finds the pool full, it drains the queues, runs `GC.gc(false)`, drains again
and retries; if the pool is still full it runs a full `GC.gc(true)` and retries once more; then
it throws `PoolExhaustedError`. `PoolExhaustedError` therefore means that more vectors are
*reachable* than the pool holds.

Consequences:

- **Size the pool generously.** Each retry costs a young collection (~100 µs). With a pool of
  400, a loop of allocating operations pays about 290 ns per operation for these collections;
  with 10 000 slots, about 10 ns. A slot costs `full_length * 8` bytes (`full_length` is the
  number of monomials), for example 64 KB at 6 variables and order 10.
- **Use the in-place forms in hot loops.** They take no new slot, so they never trigger a
  collection.
- `MiraDAC.drain!()` frees every queued vector now.
- An allocating operation pays for a finalizer registration and, because its result is freed
  only after a collection, a cold pool slot: about 200–400 ns over C++ for small vectors.
