# Points the MiraDAC.jl package at the C API library of this checkout's build
# (cmake -DDA_BUILD_CAPI=ON, see DEVELOPMENT_PLAN_JULIA.md), in the package's own
# environment and in the benchmark environment julia/MiraDAC/bench.
# Run: julia julia/dev_setup.jl  (any active project)
using Pkg
Pkg.activate(joinpath(@__DIR__, "MiraDAC"))
Pkg.instantiate()

using Preferences, UUIDs

lib = normpath(joinpath(@__DIR__, "..", "build", "capi", "libmiradac_c.so"))
isfile(lib) || error("$lib not found; build with -DDA_BUILD_CAPI=ON first")
uuid = UUID(Pkg.TOML.parsefile(joinpath(@__DIR__, "MiraDAC", "Project.toml"))["uuid"])
set_preferences!(uuid, "libmiradac" => lib; force=true)
# The bench environment reads its own LocalPreferences.toml.
Pkg.activate(joinpath(@__DIR__, "MiraDAC", "bench"))
Pkg.instantiate()
set_preferences!(uuid, "libmiradac" => lib; force=true)
println("MiraDAC: libmiradac = $lib")
