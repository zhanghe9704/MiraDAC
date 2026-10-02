"""
    MiraDAC

Julia binding of the MiraDAC differential algebra library: numeric (`NDA`), complex (`CNDA`),
symbolic (`SDA`) and complex symbolic (`CSDA`) DA vectors, their maps and DA environments
(`DAEnv`), through the C API library `libmiradac_c` (see `set_library!`).
"""
module MiraDAC

using Libdl
using LinearAlgebra
using Preferences

export MiraDACError, EnvError, PoolExhaustedError
export DAEnv, init!, clear!, current_env, default_env, with_env, with_order, get_eps, set_eps!
export import_vec
export NDA, davar, con, coeffs, coeff, nterms, erf, add!, sub!, mul!, div!
export NDAList, der, integ, substitute, compose, inv_map, evaluate_map, exponents
export CNDA, CNDAList
export SymExpr, dasymbols, simplify, SDA, sdavar, promote_sda, SDAList, evaluate
export CSDA, CSDAList
export dascope, keep!, FreedObjectError

# Bumped together with MDAC_ABI_VERSION in capi/include/miradac.h.
const ABI_VERSION = 1

# A constant, so `ccall` resolves it once; changing the preference recompiles the package.
const libmiradac = @load_preference("libmiradac", "")

"""
    set_library!(path)

Use the C API library at `path` (`libmiradac_c.so`). Takes effect after restarting Julia.
"""
function set_library!(path::AbstractString)
    @set_preferences!("libmiradac" => abspath(path))
    @info "MiraDAC: library set to $(abspath(path)); restart Julia to use it."
end

"""
    c_version()

Version string of the loaded C API library (the MiraDAC project version).
"""
c_version() = mdac_version()

include("capi.jl")
include("env.jl")
include("nda.jl")
include("algorithms.jl")
include("cnda.jl")
include("symbolic.jl")
include("csda.jl")
include("scope.jl")

# A literal exponent (`x^-1`) would otherwise call `inv`, which these types do not have.
Base.literal_pow(::typeof(^), x::Union{NDA,CNDA,SymExpr,SDA,CSDA}, ::Val{p}) where {p} = x^p

"""True if the C API library was built with symbolic support (`SymExpr`, `SDA`)."""
HAS_SYMBOLIC::Bool = false

for fn in NDA_FUNCS
    @eval export $(Symbol(fn, :!))
end

function __init__()
    isempty(libmiradac) && error("MiraDAC: no C API library configured. Run " *
        "`julia julia/dev_setup.jl` in a MiraDAC checkout after building with " *
        "-DDA_BUILD_CAPI=ON, or set the \"libmiradac\" preference of MiraDAC to the path " *
        "of libmiradac_c.")
    Libdl.dlopen(libmiradac; throw_error=false) === nothing &&
        error("MiraDAC: cannot load the C API library $libmiradac")
    abi = mdac_abi_version()
    abi == ABI_VERSION || error("MiraDAC: the C API library $libmiradac has ABI version " *
        "$abi, but this package needs ABI version $ABI_VERSION")
    global HAS_SYMBOLIC = mdac_has_symbolic() != 0
    return nothing
end

end # module
