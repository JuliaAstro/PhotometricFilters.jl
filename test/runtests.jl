using PhotometricFilters
using Test, SafeTestsets

# Run doctests
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(PhotometricFilters, :DocTestSetup, :(using PhotometricFilters); recursive=true)
doctest(PhotometricFilters)

@testset "PhotometricFilters.jl" begin
    # Write your tests here.
end
