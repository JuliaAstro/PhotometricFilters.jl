using RecipesBase

@recipe function f(filt::PhotometricFilter)
    label --> filtername(filt)
    yguide --> "throughput"
    xguide --> "wavelength"

    wave(filt), throughput(filt)
end

@recipe function f(filts::Vector{<:PhotometricFilter})
    yguide --> "throughput"
    xguide --> "wavelength"

    for filt in filts
        @series begin
            label --> filtername(filt)
            wave(filt), throughput(filt)
        end
    end
end
