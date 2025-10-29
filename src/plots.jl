using RecipesBase: @recipe, @series

@recipe function f(filt::AbstractFilter)
    label --> filtername(filt)
    yguide --> "throughput"
    xguide --> "wavelength"

    wavelength(filt), throughput(filt)
end

@recipe function f(filts::Vector{<:AbstractFilter})
    yguide --> "throughput"
    xguide --> "wavelength"

    for filt in filts
        @series begin
            label --> filtername(filt)
            wavelength(filt), throughput(filt)
        end
    end
end
