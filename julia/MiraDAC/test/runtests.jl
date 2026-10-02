using MiraDAC
using LinearAlgebra
using Test

@testset "MiraDAC" begin
    @testset "library" begin
        version = match(r"^version = \"(.*)\"$"m,
                        read(joinpath(@__DIR__, "..", "Project.toml"), String))[1]
        @test MiraDAC.c_version() == version
        # T8.2: a docstring on every exported name (Docs.undocumented_names is Julia >= 1.11).
        isdefined(Docs, :undocumented_names) && @test isempty(Docs.undocumented_names(MiraDAC))
    end
    include("test_capi.jl")
    include("test_nda.jl")
    include("test_algorithms.jl")
    include("test_cnda.jl")
    if MiraDAC.HAS_SYMBOLIC
        include("test_expr.jl")
        include("test_sda.jl")
        include("test_symbolic.jl")
        include("test_csda.jl")
        include("test_symengine_ext.jl")
        include("test_coexistence.jl")
    else
        @info "MiraDAC built without symbolic support: symbolic tests skipped"
    end
    include("test_env.jl")
    include("test_multienv.jl")
    include("test_memory.jl")
    include("test_examples.jl")
end
