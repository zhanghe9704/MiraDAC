# Getting started

## Building the C API library

The package needs `libmiradac_c`, built from a MiraDAC checkout with the CMake option
`DA_BUILD_CAPI=ON`. Symbolic support (`SymExpr`, `SDA`, `CSDA`) needs the pinned SymEngine:

```bash
scripts/setup_symengine.sh                      # once: builds the pinned SymEngine
eval "$(scripts/setup_symengine.sh --print-env)"
cmake -S . -B build -G Ninja -DDA_BUILD_CAPI=ON -DSymEngine_DIR="$SymEngine_DIR"
cmake --build build                             # -> build/capi/libmiradac_c.so
```

A numeric-only library (`NDA`, `CNDA`, envs and the map functions; no SymEngine needed) is built
with `-DWITH_SYMBOLIC=OFF` instead; then `MiraDAC.HAS_SYMBOLIC` is `false` and the symbolic
functions throw.

## Pointing the package at the library

```bash
julia julia/dev_setup.jl
```

instantiates `julia/MiraDAC` and stores the absolute path of `build/capi/libmiradac_c.so` as the
`libmiradac` preference of the package (in `LocalPreferences.toml`). Then

```bash
julia --project=julia/MiraDAC -e 'using MiraDAC; println(MiraDAC.c_version())'
julia --project=julia/MiraDAC -e 'using Pkg; Pkg.test()'
```

`MiraDAC.set_library!(path)` changes the path of a working setup; it takes effect after a
restart. When the package loads it checks the ABI version of the library and fails with the
path and both versions if they differ.

## A first session

```julia
using MiraDAC

init!(4, 3, 10_000)            # the default env: order 4, 3 variables, 10 000 slots
x = 1.0 + davar(1) + 2davar(2) # NDA
y = exp(x)
y                              # the full coefficient table (C++ operator<<)
con(y), coeff(y, [1, 1, 0])    # constant term, coefficient of x₁x₂
clear!()                       # retire the default env
```

`examples/` holds Julia ports of the C++ examples; run one with
`julia --project=julia/MiraDAC julia/MiraDAC/examples/examples.jl`.

## Threads

One env may be used by one Julia thread at a time (the C++ pools have no locks); different
threads may use different envs. Symbolic objects are for one thread at a time, since SymEngine's
reference counts are not thread safe. Neither is enforced.
