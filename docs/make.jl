using Documenter
using QuasiStableColors

makedocs(
    sitename="QuasiStableColors.jl",
    format=Documenter.HTML(),
    modules=[QuasiStableColors],
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "API Reference" => "api.md",
        "Internals" => "internals.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="github.com/mkyl/QuasiStableColors.jl.git"
)
