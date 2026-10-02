# MiraDAC_jll recipe

`build_tarballs.jl` builds `libmiradac_c` (the C API) of the tagged MiraDAC commit with the pinned
SymEngine (`cmake/symengine_pin.txt`) linked statically and hidden, and GMP from `GMP_jll`
(DEVELOPMENT_PLAN_JULIA.md, T9.1). The tarball holds `lib/libmiradac_c.so`, `include/miradac.h`
and the license.

## Local build (x86_64-linux-gnu)

```sh
S=/tmp/miradac-jll                     # any scratch directory
mkdir -p $S/bbenv $S/depot
export JULIA_DEPOT_PATH=$S/depot:$HOME/.julia   # shards, artifacts and the JLL go to $S/depot
export BINARYBUILDER_RUNNER=userns     # needs unprivileged user namespaces; "docker" also
                                       # works but runs `sudo chown` after each step
julia --project=$S/bbenv -e 'using Pkg; Pkg.add("BinaryBuilder")'
cd $S && julia --project=$S/bbenv /path/to/MiraDAC/julia/binarybuilder/build_tarballs.jl \
    --verbose --deploy=local
```

Tarballs land in `$S/products/`, the JLL package in `$S/depot/dev/MiraDAC_jll` (the first depot's
`dev/`; `JULIA_PKG_DEVDIR` is not used).

## Test MiraDAC.jl against the JLL (no `libmiradac` preference)

`Pkg.test` copies the package's `LocalPreferences.toml` (written by `dev_setup.jl`), whose
preference overrides the JLL, so move it away first. From the MiraDAC checkout:

```sh
mv julia/MiraDAC/LocalPreferences.toml $S/                 # if present
julia --project=$S/testenv -e 'using Pkg;
    Pkg.develop(path=joinpath(DEPOT_PATH[1], "dev", "MiraDAC_jll"));
    Pkg.develop(path="julia/MiraDAC"); Pkg.test("MiraDAC")'
```

The test output names the library in use (`@info "MiraDAC: C API library"`): here the file in
`$S/depot/artifacts/<tree hash>/lib/`.

## Back to the development library

Until `MiraDAC_jll` is registered (T9.2), every environment with MiraDAC.jl needs the local JLL,
which `dev_setup.jl` adds when `MIRADAC_JLL` is set:

```sh
MIRADAC_JLL=$S/depot/dev/MiraDAC_jll julia julia/dev_setup.jl
julia --project=julia/MiraDAC -e 'using Pkg; Pkg.test()'
```

## Actions runtime tests

Upload the `products/MiraDAC.v*.tar.gz` files as artifacts of an Actions run
**in this repository**, for example with `actions/upload-artifact@v4` using
`path: products/MiraDAC.v*.tar.gz`. Dispatch **JLL runtime tests** with that
run's numeric ID (the number at the end of its Actions URL) as `run_id`.
The matrix tests available platforms on Julia 1.10 and 1, logging `SKIP` for
missing tarballs. The workflow creates a temporary local JLL from each tarball;
it does not require registration or publish a JLL. Keep
`julia/MiraDAC/LocalPreferences.toml` absent so the preference cannot override it.
