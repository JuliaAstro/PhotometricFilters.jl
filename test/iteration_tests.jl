using PhotometricFilters
using Test

for f in (PhotometricFilters.SDSS_u(), get_filter("2MASS/2MASS.J", :Vega))
    @test f[firstindex(f)] isa Tuple{<:Number, <:Number}
    @test f[lastindex(f)] isa Tuple{<:Number, <:Number}
    @test eltype(f) == typeof(f[1])
    @test length(f) == length(wavelength(f)) == length(throughput(f))
    @test size(f) == size(wavelength(f)) == size(throughput(f))
    for (i, val) in enumerate(f)
        @test val == (wavelength(f)[i], throughput(f)[i])
    end
end