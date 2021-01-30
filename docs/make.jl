using PhotometricFilters
using Documenter

makedocs(;
    modules=[PhotometricFilters],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/juliaastro/PhotometricFilters.jl/blob/{commit}{path}#L{line}",
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
    repo="github.com/juliaastro/PhotometricFilters.jl",
)
