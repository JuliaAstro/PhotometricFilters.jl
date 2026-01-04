using PhotometricFilters
using Test, SafeTestsets

# Run doctests
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(PhotometricFilters, :DocTestSetup, :(using PhotometricFilters; using Logging; disable_logging(Logging.Info)); recursive=true)
doctest(PhotometricFilters)

@testset verbose=true "PhotometricFilters.jl" begin
    @safetestset "Library" include("library_tests.jl")
    @safetestset "Iteration" include("iteration_tests.jl")
    @safetestset "SVO tests" include("svo_tests.jl")
end
