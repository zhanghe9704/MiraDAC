# MiraDAC_jll recipe

`build_tarballs.jl` builds `libmiradac_c` (the C API) with the pinned SymEngine
(`cmake/symengine_pin.txt`) linked statically and hidden, and GMP from `GMP_jll`
(DEVELOPMENT_PLAN_JULIA.md, T9.1). Targets: `aarch64-linux-gnu`, `x86_64-linux-gnu`,
`x86_64-linux-musl`, `aarch64-linux-musl`, `armv7l-linux-gnueabihf`,
`powerpc64le-linux-gnu`, `riscv64-linux-gnu`, `x86_64-unknown-freebsd`,
`aarch64`/`x86_64` Apple and `x86_64`/`i686` Windows. The tarballs hold
`lib/libmiradac_c.{so,dylib}` / `bin/libmiradac_c.dll`, `include/miradac.h` and the license.

## Local build (all platforms)

```sh
S=/tmp/miradac-jll                     # any scratch directory
mkdir -p $S/bbenv $S/depot
export JULIA_DEPOT_PATH=$S/depot:$HOME/.julia   # shards, artifacts and the JLL go to $S/depot
export BINARYBUILDER_RUNNER=userns     # needs unprivileged user namespaces; "docker" also
                                       # works but runs `sudo chown` after each step
export BINARYBUILDER_AUTOMATIC_APPLE=true  # accept the Apple SDK terms; needed for the
                                           # macOS legs (first darwin build downloads the SDK)
julia --project=$S/bbenv -e 'using Pkg; Pkg.add("BinaryBuilder")'
cd $S && julia --project=$S/bbenv /path/to/MiraDAC/julia/binarybuilder/build_tarballs.jl \
    --verbose --deploy=local           # add a platform triplet to build one leg only, e.g.
                                       # aarch64-apple-darwin
```

If `/tmp` is tight, also `export TMPDIR=<dir-on-a-big-disk>`: the first darwin build stages
gigabytes of toolchain and SDK downloads through it.

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

Create a published prerelease with the built tarballs attached (the
`MiraDAC-logs...` tarballs are not included by this pattern):

```sh
gh release create <tag> --prerelease --title "MiraDAC JLL test binaries" \
  products/MiraDAC.v*.tar.gz
```

Draft releases are not supported. The workflow downloads matching tarballs from
that release and tests available platforms on Julia 1.10 and 1, logging `SKIP`
for missing platforms. It creates a temporary local JLL from each tarball, so no
JLL registration is needed. Keep `julia/MiraDAC/LocalPreferences.toml` absent so
the preference cannot override the JLL.

`workflow_dispatch` is available only after this workflow file exists on the
default branch (`main`). Then dispatch the test branch with:

```sh
gh workflow run jll-tests.yml --ref jll -f tag=<tag>
```

Replace `<tag>` with the tag used when creating the prerelease above.
