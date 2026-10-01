using MiraDAC
using LinearAlgebra
using Test

@testset "MiraDAC" begin
    @testset "library" begin
        version = match(r"^version = \"(.*)\"$"m,
                        read(joinpath(@__DIR__, "..", "Project.toml"), String))[1]
        @test MiraDAC.c_version() == version
    end
    include("test_capi.jl")
    include("test_nda.jl")
    include("test_algorithms.jl")
    include("test_cnda.jl")
    include("test_env.jl")
    include("test_memory.jl")
end
