import PhotometricFilters
using Test

@testset "Library" begin
    for filter_name in PhotometricFilters.FILTER_NAMES
        F = Symbol(filter_name)
        filt = @eval PhotometricFilters.$F()
        @test filt isa PhotometricFilters.PhotometricFilter
    end
end
