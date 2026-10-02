# Points the MiraDAC.jl package at the C API library of this checkout's build
# (cmake -DDA_BUILD_CAPI=ON, see DEVELOPMENT_PLAN_JULIA.md), in the package's own
# environment and in its benchmark and docs environments julia/MiraDAC/{bench,docs}.
# Run: julia julia/dev_setup.jl  (any active project)
using Pkg
pkg = joinpath(@__DIR__, "MiraDAC")
# Until MiraDAC_jll is registered (plan T9.2), its local build (julia/binarybuilder/README.md)
# must be put in each environment: MIRADAC_JLL=/path/to/MiraDAC_jll julia julia/dev_setup.jl
# Only the (untracked) Manifest.toml keeps it: on Julia >= 1.11 `develop` also writes a
# machine-local [sources] path into the tracked Project.toml, which is restored here.
function jll()
    haskey(ENV, "MIRADAC_JLL") || return
    project = Base.active_project()
    saved = read(project, String)
    Pkg.develop(path=ENV["MIRADAC_JLL"])
    write(project, saved)
end
Pkg.activate(pkg)
jll()
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
    jll()
    VERSION < v"1.11" && Pkg.develop(path=pkg)
    Pkg.instantiate()
    set_preferences!(uuid, "libmiradac" => lib; force=true)
end
println("MiraDAC: libmiradac = $lib")
