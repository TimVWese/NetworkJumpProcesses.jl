using Documenter
using NetworkJumpProcesses

makedocs(
    sitename = "NetworkJumpProcesses",
    format = Documenter.HTML(),
    modules = [NetworkJumpProcesses],
    pages = [
        "Home" => "index.md",
        "Examples" => Any[
            "examples/sis.md",
            "examples/reaction_diffusion.md",
        ],
        "Guide" => "guide.md",
        "Reference" => "reference.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/TimVWese/NetworkJumpProcesses.jl",
)
