using Interpolations: linear_interpolation, deduplicate_knots!
using Trapz: trapz
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

struct PhotometricFilter{T, WT, DT <: DetectorType, TT <: AbstractVector{T},
                         ST <: Union{String, Nothing}, ET} <: AbstractFilter{T}
    wave::WT
    throughput::TT
    detector::DT
    name::ST
    etp::ET
end

"""Default unit for wavelengths when not otherwise specified."""
const wave_unit = Unitful.angstrom

# Dispatch method for PhotometricFilter constructor
_convert_wave(w) = w
_convert_wave(w::Unitful.Length) = ustrip(wave_unit, w)
_convert_wave(w::Quantity) = throw(ArgumentError("Provided wavelengths have incompatible dimensions -- expected length (𝐋), received $(dimension(w))."))

"""
    PhotometricFilter(wave::AbstractVector, throughput::AbstractVector{T};
                      detector::DetectorType=Photon(), name::Union{String, Nothing}=nothing)
Struct representing a photometric filter, defined by vectors of wavelengths (`wave`) and filter throughputs (`throughput`).
`wave` can have `Unitful` units attached, otherwise they are assumed to be $wave_unit.
Optional keyword arguments define the detector type for which the filter is valid and a name to identify the filter.
```jldoctest
julia> using PhotometricFilters: PhotometricFilter, Photon, wave, throughput

julia> using Unitful

julia> f = PhotometricFilter(1000:2000, vcat(fill(0.25, 250), fill(0.5, 500), fill(0.25, 251))) # Specify only wavelength and throughput
1001-element PhotometricFilter{Float64}: nothing
 min. wave.: 1000 Å
 max. wave.: 2000 Å
 effective wave.: 1603.6927025575474 Å
 mean wave.: 1499.8333333333333 Å
 central wave.: 1499.5 Å
 pivot wave.: 1478.1028279485677 Å
 eff. width: 750.0 Å
 fwhm: 501.0 Å

julia> f == PhotometricFilter(uconvert.(Unitful.nm, wave(f)), throughput(f)) # Can also specify wavelength argument with Unitful units
true

julia> f[10] # Indexing into the filter as `f[i]` returns `throughput(f)[i]`
0.25

julia> f(1001.1) # Calling `f` like a function interpolates the throughput
0.25

julia> f(100.11 * Unitful.nm) # Can also specify wavelength with units
0.25
```
"""
function PhotometricFilter(wave::AbstractVector, throughput::AbstractVector{T};
                           detector::DetectorType=Photon(), name::Union{String, Nothing}=nothing) where T
    if length(wave) != length(throughput)
        throw(ArgumentError("Wavelength and throughput arrays must have equal length"))
    end
    bc = zero(T)
    wv = _convert_wave.(wave)
    deduplicate_knots!(wv; move_knots=true) # Ensure wave entries are all unique
    etp = linear_interpolation(wv, throughput; extrapolation_bc=bc)
    return PhotometricFilter(wv, throughput, detector, name, etp)
end

Base.getindex(f::PhotometricFilter, idx...) = throughput(f)[idx...]

Base.show(io::IO, f::PhotometricFilter) = print(io, f.name)

function Base.show(io::IO, ::MIME"text/plain", f::PhotometricFilter{T}) where T
    N = length(f)
    println(io, "$N-element PhotometricFilter{$T}: ", f.name)
    println(io, " min. wave.: ", min_wave(f))
    println(io, " max. wave.: ", max_wave(f))
    println(io, " effective wave.: ", effective_wavelength(f))
    println(io, " mean wave.: ", mean_wavelength(f))
    println(io, " central wave.: ", central_wavelength(f))
    println(io, " pivot wave.: ", pivot_wavelength(f))
    println(io, " eff. width: ", width(f))
    print(io,   " fwhm: ", fwhm(f))
end

wave(f::PhotometricFilter) = f.wave .* wave_unit
throughput(f::PhotometricFilter) = f.throughput

(f::PhotometricFilter)(wave) = f.etp(wave)

function (f::PhotometricFilter)(wave::Q) where Q <: Unitful.Length
    wl = ustrip.(wave_unit, wave)
    return f(wl)
end

Base.size(f::PhotometricFilter) = size(throughput(f))

"""
    effective_wavelength(f::PhotometricFilter)

Returns the effective wavelength of the filter `f` using the Vega spectrum as a standard. Defined as

```math
\\frac{\\int \\lambda \\, T(\\lambda) \\text{Vg}(\\lambda) \\, d\\lambda}{\\int T(\\lambda) \\text{Vg}(\\lambda) \\, d\\lambda}
```
where ``T(\\lambda)`` is the filter transmission at wavelength ``\\lambda`` and ``\\text{Vg}(\\lambda)`` is the spectrum of Vega.
"""
function effective_wavelength(f::PhotometricFilter)
    wvega, fvega = Vega()
    wl = uconvert.(wave_unit, wvega)
    filt = f.(wl) .* fvega
    norm = trapz(wl, filt)
    leff = trapz(wl, wl .* filt)
    return leff / norm
end

"""
    pivot_wavelength(f::PhotometricFilter)

Returns the pivot wavelength of the filter `f`, defined for filters with `Energy` detector types as

```math
\\sqrt{ \\frac{\\int T(\\lambda) \\, d\\lambda}{\\int T(\\lambda) / \\lambda^2 \\, d\\lambda} }
```

For filters with `Photon` detector types, ``\\lambda \\, T(\\lambda)`` is substituted for ``T(\\lambda)`` in the above expression.

Internally integration is carried out using trapezoidal integration. It can be convenient to think of this as the "center of mass" of the filter.
"""
function pivot_wavelength(f::PhotometricFilter{T, S, <:Energy}) where {T, S}
    wl = wave(f)
    y = throughput(f) ./ wl.^2
    norm = trapz(wl, throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end
function pivot_wavelength(f::PhotometricFilter{T, S, <:Photon}) where {T, S}
    wl = wave(f)
    y = throughput(f) ./ wl
    norm = trapz(wl, wl .* throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end


"""
    mean_wavelength(f::PhotometricFilter)

Returns the mean wavelength of the filter `f`, defined as

```math
\\frac{\\int \\lambda \\, T(\\lambda) \\, d\\lambda}{\\int T(\\lambda) \\, d\\lambda}
```
"""
function mean_wavelength(f::PhotometricFilter)
    wl = wave(f)
    norm = trapz(wl, throughput(f))
    lt = trapz(wl, wl .* throughput(f))
    return  lt / norm
end

"""
    central_wavelength(f::PhotometricFilter)

Returns the central wavelength of the filter `f`, defined as the central wavelength between the two wavelengths used for the FWHM ([`fwhm`](@ref)).
"""
function central_wavelength(f::PhotometricFilter)
    y = throughput(f)
    wl = wave(f)
    thresh = maximum(y) / 2

    # lower FWHM wavelength
    idx_low = findfirst(>(thresh), y)
    if idx_low == firstindex(y)
        wl_low = wl[idx_low]
    else
        frac_low = (thresh - y[idx_low-1]) / (y[idx_low] - y[idx_low-1])
        wl_low = frac_low * wl[idx_low] + (1 - frac_low) * wl[idx_low-1]
    end

    # upper FWHM wavelength
    idx_high = findlast(>(thresh), y)
    if idx_high == lastindex(y)
        wl_high = wl[idx_high]
    else
        frac_high = (thresh - y[idx_high+1]) / (y[idx_high] - y[idx_high+1])
        wl_high = frac_high * wl[idx_high] + (1 - frac_high) * wl[idx_high+1]
    end

    return (wl_low + wl_high) / 2
end

"""
    min_wave(f::PhotometricFilter; level=0.01)

Returns the shortest wavelength at which the filter transmission is equal to `level * maximum(transmission)`.
"""
function min_wave(f::PhotometricFilter; level=0.01)
    y = throughput(f)
    wl = wave(f)
    thresh = level * maximum(y)
    idx = findfirst(>(thresh), y)
    if idx == firstindex(y)
        return wl[idx]
    end
    frac = (thresh - y[idx-1]) / (y[idx] - y[idx-1])
    return frac * wl[idx] + (1 - frac) * wl[idx-1]
end

"""
    max_wave(f::PhotometricFilter; level=0.01)

Returns the longest wavelength at which the filter transmission is equal to `level * maximum(transmission)`.
"""
function max_wave(f::PhotometricFilter; level=0.01)
    y = throughput(f)
    wl = wave(f)
    thresh = level * maximum(y)
    idx = findlast(>(thresh), y)
    if idx == lastindex(y)
        return wl[idx]
    end
    frac = (thresh - y[idx+1]) / (y[idx] - y[idx+1])
    return frac * wl[idx] + (1 - frac) * wl[idx+1]
end

"""
    fwhm(f::PhotometricFilter)

Returns the difference between the furthest two wavelengths for which the filter transmission is equal to half its maximum value.
"""
function fwhm(f::PhotometricFilter)
    y = throughput(f)
    wl = wave(f)
    thresh = maximum(y) / 2

    # lower FWHM wavelength
    idx_low = findfirst(>(thresh), y)
    if idx_low == firstindex(y)
        wl_low = wl[idx_low]
    else
        frac_low = (thresh - y[idx_low-1]) / (y[idx_low] - y[idx_low-1])
        wl_low = frac_low * wl[idx_low] + (1 - frac_low) * wl[idx_low-1]
    end

    # upper FWHM wavelength
    idx_high = findlast(>(thresh), y)
    if idx_high == lastindex(y)
        wl_high = wl[idx_high]
    else
        frac_high = (thresh - y[idx_high+1]) / (y[idx_high] - y[idx_high+1])
        wl_high = frac_high * wl[idx_high] + (1 - frac_high) * wl[idx_high+1]
    end

    return wl_high - wl_low
end

"""
    width(f::PhotometricFilter)

Returns the effective width of the filter, defined as the horizontal size of a rectangle with height equal to the maximum transmission of the filter such that the area of the rectangle is equal to the area under the filter transmission curve. This is calculated as

```math
\\frac{\\int T(\\lambda) \\, d\\lambda}{\\text{max}(T(\\lambda))}
```
"""
function width(f::PhotometricFilter)
    norm = trapz(wave(f), throughput(f))
    return norm / maximum(throughput(f))
end

"""
    apply(f::PhotometricFilter, wave, flux)

Use linear interpolation to map the wavelengths of the photometric filter `f` to the given wavelengths `wave` and apply the filter throughput to the `flux`. The wavelengths of the filter and `wave` need to be compatible. This means if one has units, the other one needs units, too.

```jldoctest
julia> using PhotometricFilters: SDSS_u, wave_unit, apply

julia> f = SDSS_u();

julia> λ = 3000:4000
3000:4000

julia> flux = fill(1.0, length(λ)); # If `flux` is all `1`, `apply` reduces to `f` interpolated at `λ`

julia> apply(f, λ, flux) == f(λ)
true

julia> λ_u = λ .* wave_unit # Can also put units on λ
(3000:4000) Å

julia> apply(f, λ_u, flux) == f.(λ_u)
true
```
"""
apply(filt::PhotometricFilter, wave, flux) = apply!(filt, wave, flux, similar(flux))

"""
    apply!(f::PhotometricFilter, wave, flux, out)

In-place version of [`apply`](@ref) which modifies `out`. It should have a compatible element type with `flux`.
"""
function apply!(filt::PhotometricFilter, wave, flux, out)
    @. out = flux * filt(wave)
    return out
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
