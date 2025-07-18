using PhotometricFilters
using Test

for f in (PhotometricFilters.SDSS_u(), get_filter("2MASS/2MASS.J", :Vega))
    @test f[firstindex(f)] isa Tuple{<:Number, <:Number}
    @test f[lastindex(f)] isa Tuple{<:Number, <:Number}
    @test eltype(f) == typeof(f[1])
    @test length(f) == length(wave(f)) == length(throughput(f))
    @test size(f) == size(wave(f)) == size(throughput(f))
    for (i, val) in enumerate(f)
        @test val == (wave(f)[i], throughput(f)[i])
    end
end