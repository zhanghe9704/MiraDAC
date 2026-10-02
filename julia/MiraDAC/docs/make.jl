# Builds the MiraDAC.jl documentation into docs/build (no deploy).
# Run: julia --project=julia/MiraDAC/docs julia/MiraDAC/docs/make.jl (after julia/dev_setup.jl)
using Documenter, MiraDAC

makedocs(;
    sitename = "MiraDAC.jl",
    modules = [MiraDAC],
    # No public repository: no repository or edit links.
    format = Documenter.HTML(; prettyurls = false, repolink = nothing, edit_link = nothing),
    remotes = nothing,
    warnonly = false,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "NDA" => "nda.md",
        "CNDA" => "cnda.md",
        "Symbolic" => "symbolic.md",
        "Multiple envs" => "envs.md",
        "Performance" => "performance.md",
        "API reference" => "api.md",
    ],
)
