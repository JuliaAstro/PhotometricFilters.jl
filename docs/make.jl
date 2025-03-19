using PhotometricFilters
using Documenter
using Documenter.Remotes: GitHub

makedocs(;
    modules=[PhotometricFilters],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo=GitHub("JuliaAstro/PhotometricFilters.jl"),
    sitename="PhotometricFilters.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/PhotometricFilters.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/PhotometricFilters.jl",
)
