# T8.1: every script in examples/ runs to completion, each in its own Julia process (each calls
# init! and clear!) and in a temporary directory (examples.jl writes da_output.txt).
const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const SYMBOLIC_EXAMPLES = ("example_1_symbolic.jl", "example_complex_symbolic.jl",
                           "example_interop.jl")

@testset "examples" begin
    names = sort(filter(endswith(".jl"), readdir(EXAMPLES_DIR)))
    @test names == ["example_1_symbolic.jl", "example_complex_da.jl",
                    "example_complex_symbolic.jl", "example_interop.jl", "examples.jl"]
    for name in names
        name in SYMBOLIC_EXAMPLES && !MiraDAC.HAS_SYMBOLIC && continue
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $(joinpath(EXAMPLES_DIR, name))`
        out = mktempdir() do dir
            p = run(pipeline(ignorestatus(Cmd(cmd; dir)); stdout=(io = IOBuffer()), stderr=io))
            (p.exitcode, String(take!(io)))
        end
        @test (out[1] == 0 && !occursin("FAIL", out[2])) || (println(name, ":\n", out[2]); false)
    end
end
