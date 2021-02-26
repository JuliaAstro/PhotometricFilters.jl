using Interpolations

abstract type AbstractFilter end

abstract type DetectorType end
struct Photon <: DetectorType end
struct Energy <: DetectorType end

function Base.parse(::Type{DetectorType}, s::AbstractString)
    tok = lowercase(s)
    if tok === "photon"
        return Photon()
    elseif tok === "energy"
        return Energy()
    else
        throw(ArgumentError("failed to parse DetectorType from $s"))
    end
end

struct PhotometricFilter{WLT,TT,DT<:DetectorType,ST<:Union{String,Nothing}} <: AbstractFilter
    wave::WLT
    throughput::TT
    detector::DT
    name::ST
end

function PhotometricFilter(wave::AbstractVector, throughput::AbstractVector; detector=Photon(), name=nothing)
    length(wave) == length(throughput) || error("wavelength and throughput arrays must match")
    return PhotometricFilter(wave, throughput, detector, name)
end

Base.show(io::IO, f::PhotometricFilter) = print(io, f.name)

function Base.show(io::IO, ::MIME"text/plain", f::PhotometricFilter)
    print(io, "PhotometricFilter: $(f.name)\n wave: ", f.wave, "\n throughput: ", f.throughput)
end

wave(f::PhotometricFilter) = f.wave
throughput(f::PhotometricFilter) = f.throughput

Base.size(f::PhotometricFilter) = size(throughput(f))
Base.length(f::PhotometricFilter) = prod(size(f))

Base.iterate(f::PhotometricFilter, args...) = iterate(wave(f), args...), iterate(throughput(f), args...)


function central_wavelength()

end

# """
#     effective_wavelength(::PhotometricFilter)

# Return the effective wavelength
# """
# function effective_wavelength()
# end

function pivot_wavelength()
end

function norm()

end

"""
    apply(::PhotometricFilter, wave, flux)

Use linear interpolation to map the wavelenghts of the photometric filter to the given wavelengths and apply the filter throughput to the `flux`. The wavelengths of the filter and `wave` need to be compatible. This means if one has units, the other one needs units, too.

# See also
[`interpolator`](@ref)
"""
apply(filt::PhotometricFilter, wave, flux) = apply!(filt, wave, flux, similar(flux))

"""
    apply!(::PhotometricFilter, wave, flux, out)

In-place version of [`apply`](@ref) which modifies `out`. It should have a compatible element type with `flux`.
"""
function apply!(filt::PhotometricFilter, wave, flux, out)
    etp = interpolator(filt)
    @. out = flux * etp(wave)
    return out
end

function interpolator(filt::PhotometricFilter)
    x = wave(filt) |> ustrip
    y = throughput(filt)
    bc = zero(eltype(y))
    return LinearInterpolation(x, y; extrapolation_bc=bc)
end

function fwhm(filt::PhotometricFilter)
    y = throughput(filt)
    Δ = diff(sign.(y ./ maximum(y) .- 0.5))
    i1 = findfirst(!iszero, Δ)
    i2 = findnext(!iszero, Δ, i1 + 1)
    return wave(filt)[i2] - wave(filt)[i1]
end

function Base.:+(f1::PhotometricFilter, f2::PhotometricFilter)
    # find total extent and log-spacing
    w1 = wave(f1)
    w2 = wave(f2)
    e1 = extrema(w1)
    e2 = extrema(w2)
    min_wl = min(e1[1], e2[1])
    max_wl = max(e1[2], e2[2])
    δwl1 = (e1[2] - e1[1]) / length(w1)
    δwl2 = (e2[2] - e2[1]) / length(w2)
    δ = min(δwl1, δwl2)

    wl = min_wl:δ:max_wl

    etp1 = interpolator(f1)
    etp2 = interpolator(f2)
    through = @. etp1(ustrip(wl)) + etp2(ustrip(wl))
    name = "$(f1.name) + $(f2.name)"
    return PhotometricFilter(wl, through, f1.detector, name)
end
