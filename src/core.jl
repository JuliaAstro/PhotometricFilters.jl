using Interpolations
using Trapz
using Unitful

abstract type AbstractFilter{T} <: AbstractVector{T} end

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

struct PhotometricFilter{T,WT,DT<:DetectorType,TT<:AbstractVector{T},ST<:Union{String,Nothing},ET} <: AbstractFilter{T}
    wave::WT
    throughput::TT
    detector::DT
    name::ST
    etp::ET
end

function PhotometricFilter(wave::AbstractVector, throughput::AbstractVector{T}; detector=Photon(), name=nothing) where T
    length(wave) == length(throughput) || error("wavelength and throughput arrays must match")
    bc = zero(T)
    etp = LinearInterpolation(ustrip(wave), throughput; extrapolation_bc=bc)
    return PhotometricFilter(wave, throughput, detector, name, etp)
end

Base.getindex(f::PhotometricFilter, idx...) = throughput(f)[idx...]

Base.show(io::IO, f::PhotometricFilter) = print(io, f.name)

function Base.show(io::IO, ::MIME"text/plain", f::PhotometricFilter{T}) where T
    N = length(f)
    min_wl = min_wave(f)
    max_wl = max_wave(f)
    piv_wl = pivot_wavelength(f)
    cen_wl = central_wavelength(f)
    eff_wl = effective_wavelength(f)
    eff_width = width(f)
    Γ = fwhm(f)
    println(io, "$N-element PhotometricFilter{$T}: ", f.name)
    println(io, " min. wave.: ", min_wl)
    println(io, " max. wave.: ", max_wl)
    println(io, " effective wave.: ", eff_wl)
    println(io, " central wave.: ", cen_wl)
    println(io, " pivot wave.: ", piv_wl)
    println(io, " eff. width: ", eff_width)
    print(io,   " fwhm: ", Γ)
end

wave(f::PhotometricFilter) = f.wave
throughput(f::PhotometricFilter) = f.throughput

(f::PhotometricFilter)(wave) = f.etp(wave)

function (f::PhotometricFilter)(wave::Q) where Q <: Unitful.Length
    # strip units because Interpolations.jl doesn't work nicely with them
    wl = ustrip(unit(eltype(f.wave)), wave)
    return f.etp(wl)
end

Base.size(f::PhotometricFilter) = size(throughput(f))

"""
    effective_wavelength(::PhotometricFilter)

Return the effective wavelength using the Vega spectrum as a standard
"""
function effective_wavelength(f::PhotometricFilter)
    wvega, fvega = Vega()
    u = unit(eltype(wave(f)))
    wl = ustrip.(u, wvega)
    filt = f.(wvega) .* fvega
    norm = trapz(wl, filt)
    leff = trapz(wl, wl .* filt)
    return leff / norm * u
end

"""
    pivot_wavelength(::PhotometricFilter)

Returns the pivot wavelength of the filter, described by the equation below. Internally integration is carried out using trapezoidal integration. It can be convenient to think of this as the "center of mass" of the filter.
"""
function pivot_wavelength(f::PhotometricFilter{T,S,<:Photon}) where {T,S}
    wl = ustrip(wave(f))
    y = throughput(f) ./ wl
    norm = trapz(wl, wl .* throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2) * unit(eltype(wave(f)))
end

function pivot_wavelength(f::PhotometricFilter{T,S,<:Energy}) where {T,S}
    wl = ustrip(wave(f))
    y = throughput(f) ./ wl.^2
    norm = trapz(wl, throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2) * unit(eltype(wave(f)))
end


function central_wavelength(f::PhotometricFilter)
    wl = ustrip(wave(f))
    norm = trapz(wl, throughput(f))
    lt = trapz(wl, wl .* throughput(f))
    return  lt / norm * unit(eltype(wave(f)))
end

function min_wave(f::PhotometricFilter; level=0.01)
    y = throughput(f)
    idx = findfirst(q -> q / maximum(y) > level, y)
    return wave(f)[idx]
end

function max_wave(f::PhotometricFilter; level=0.01)
    y = throughput(f)
    idx = findlast(q -> q / maximum(y) > level, y)
    return wave(f)[idx]
end

"""
    apply(::PhotometricFilter, wave, flux)

Use linear interpolation to map the wavelenghts of the photometric filter to the given wavelengths and apply the filter throughput to the `flux`. The wavelengths of the filter and `wave` need to be compatible. This means if one has units, the other one needs units, too.
"""
apply(filt::PhotometricFilter, wave, flux) = apply!(filt, wave, flux, similar(flux))

"""
    apply!(::PhotometricFilter, wave, flux, out)

In-place version of [`apply`](@ref) which modifies `out`. It should have a compatible element type with `flux`.
"""
function apply!(filt::PhotometricFilter, wave, flux, out)
    @. out = flux * filt(wave)
    return out
end

function fwhm(filt::PhotometricFilter)
    Δ = diff(sign.(filt ./ maximum(filt) .- 0.5))
    nonzeros = findall(!iszero, Δ)
    i1, i2 = first(nonzeros), last(nonzeros)
    return wave(filt)[i2] - wave(filt)[i1]
end

function width(f::PhotometricFilter)
    norm = trapz(ustrip(wave(f)), throughput(f))
    return norm / maximum(throughput(f)) * unit(eltype(wave(f)))
end

function Base.:*(f1::PhotometricFilter, f2::PhotometricFilter)
    # find total extent and log-spacing
    w1 = wave(f1)
    w2 = wave(f2)
    e1 = extrema(w1)
    e2 = extrema(w2)
    min_wl = max(e1[1], e2[1])
    max_wl = min(e1[2], e2[2])
    max_wl ≤ min_wl && error("no overlap between the two filters")
    δwl1 = (e1[2] - e1[1]) / length(w1)
    δwl2 = (e2[2] - e2[1]) / length(w2)
    δ = min(δwl1, δwl2)

    wl = min_wl:δ:max_wl

    through = @. f1(wl) * f2(wl)
    name = "$(f1.name) * $(f2.name)"
    return PhotometricFilter(wl, through; detector=f1.detector, name=name)
end
