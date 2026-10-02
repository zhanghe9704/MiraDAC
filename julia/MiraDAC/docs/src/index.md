# MiraDAC.jl

Julia binding of MiraDAC, a C++ library for differential algebra (DA):
truncated multivariate Taylor polynomials with numeric (`NDA`), complex (`CNDA`), symbolic
(`SDA`, coefficients are SymEngine expressions) and complex symbolic (`CSDA`) coefficients.

The package calls the C API library `libmiradac_c` with `ccall`. DA vectors live in fixed-size
C++ pools that belong to a *DA environment* (`DAEnv`): the number of variables, the maximum
order and the pool size are fixed when the env is made.

```julia
using MiraDAC

init!(4, 3, 10_000)                 # order 4, 3 variables, 10 000 vectors
x = 1.0 + davar(1) + 2davar(2)
y = exp(x)
coeff(y, [1, 0, 0])                 # ∂y/∂x₁ at 0 = e
```

- [Getting started](@ref): building the C library and pointing the package at it.
- [NDA](@ref), [CNDA](@ref), [Symbolic](@ref): the vector types.
- [Multiple envs](@ref): several environments in one session.
- [Performance](@ref): the in-place API and how the pool and the garbage collector interact.
- [API reference](@ref): every exported name.

Indices are 1-based everywhere: `davar(1)` is the first variable and coefficient 1 of `coeffs(v)`
is the constant term. DA types are not subtypes of `Number`; their operators are defined
explicitly.
