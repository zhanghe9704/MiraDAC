# CNDA

A `CNDA` is a complex numeric DA vector `re + im*i` (C++ `std::complex<NDA>`).

```julia
init!(4, 3, 10_000)
x = davar(1) + 2davar(2)
z = CNDA(cos(x), sin(x))            # from two NDA parts of one env
w = CNDA(1.0 + 2.0im)               # a constant
real(z), imag(z)                    # copies of the parts
```

Operators combine `CNDA` with `CNDA`, `NDA`, `Real` and `Complex` on either side, and an `NDA`
with a `Complex` gives a `CNDA`. The functions are those C++ has for complex vectors: `sqrt exp
log asin acos atan asinh acosh atanh` and `^`; `abs(z)` is the larger of `norm(real(z))` and
`norm(imag(z))`. There is no `conj` (C++ has none). The in-place forms (`add!`, `mul!`, `exp!`,
…) take `CNDA` outputs.

A `CNDAList` holds a complex map. `compose(m, v)` with `CNDA`s in the map `m` or in the
arguments `v` is C++'s `cd_composition` and returns a `CNDAList`; `compose(m, pt)` with an
`NDA` map and a complex point returns a `Vector{ComplexF64}`.

```julia
m = NDAList([cos(x), sin(x)])
compose(m, [0.1 + 0.2im, 0.3im, 1.0])          # Vector{ComplexF64}
compose(m, CNDAList([z, z * z, CNDA(x)]))      # CNDAList
```
