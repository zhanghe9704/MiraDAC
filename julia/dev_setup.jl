# Points the MiraDAC.jl package at the C API library of this checkout's build
# (cmake -DDA_BUILD_CAPI=ON, see DEVELOPMENT_PLAN_JULIA.md), in the package's own
# environment and in its benchmark and docs environments julia/MiraDAC/{bench,docs}.
# Run: julia julia/dev_setup.jl  (any active project)
using Pkg
pkg = joinpath(@__DIR__, "MiraDAC")
Pkg.activate(pkg)
Pkg.instantiate()

using Preferences, UUIDs

lib = normpath(joinpath(@__DIR__, "..", "build", "capi", "libmiradac_c.so"))
isfile(lib) || error("$lib not found; build with -DDA_BUILD_CAPI=ON first")
uuid = UUID(Pkg.TOML.parsefile(joinpath(pkg, "Project.toml"))["uuid"])
set_preferences!(uuid, "libmiradac" => lib; force=true)
# The bench and docs environments read their own LocalPreferences.toml. Their [sources] entry
# needs Julia >= 1.11; `develop` gives Julia 1.10 the same path dependency.
for sub in ("bench", "docs")
    Pkg.activate(joinpath(pkg, sub))
    VERSION < v"1.11" && Pkg.develop(path=pkg)
    Pkg.instantiate()
    set_preferences!(uuid, "libmiradac" => lib; force=true)
end
println("MiraDAC: libmiradac = $lib")
