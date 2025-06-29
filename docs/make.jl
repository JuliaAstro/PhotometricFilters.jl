using PhotometricFilters
using Documenter: makedocs, HTML, deploydocs
using Documenter.Remotes: GitHub

makedocs(;
         sitename = "PhotometricFilters.jl",
         modules = [PhotometricFilters],
         authors = "Miles Lucas <mdlucas@hawaii.edu>, Chris Garling, Lucas Valenzuela, and contributors",
         repo = GitHub("JuliaAstro/PhotometricFilters.jl"),
         format = HTML(;
                       prettyurls = get(ENV, "CI", nothing) == "true",
                       canonical = "https://juliaastro.org/PhotometricFilters/stable/",
                       assets = String[],
                       ),
         pages = [
                  "Home" => "index.md",
                 ],
         doctest = false,
         linkcheck = true,
         warnonly = [:missing_docs, :linkcheck]
         )

deploydocs(;
           repo = "github.com/JuliaAstro/PhotometricFilters.jl",
           versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
           push_preview = true
           )
