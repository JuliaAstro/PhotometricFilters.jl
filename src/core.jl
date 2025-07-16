using Interpolations: linear_interpolation, deduplicate_knots!
using Trapz: trapz
using Unitful
using UnitfulAstro

"""Default unit for wavelengths."""
const wave_unit = Unitful.angstrom

# Dispatch methods for AbstractFilter constructors
_convert_wave(w) = w * wave_unit
_convert_wave(w::Unitful.Length) = uconvert(wave_unit, w)
_convert_wave(w::Quantity) = throw(ArgumentError("Provided wavelengths have incompatible dimensions -- expected length (𝐋), received $(dimension(w))."))

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

"""
    AbstractFilter{T}
Abstract supertype for representing photometric filters. Most functions provided by this package (e.g., [`effective_wavelength`](@ref PhotometricFilters.effective_wavelength) and similar methods) are designed to work with any subtype of `AbstractFilter` so long as a minimal API is defined for new subtypes. The methods that should be implemented for new types to conform to this API are summarized below:

 - [`name(f::NewType)`](@ref name) should return a string indicating a human-readable name for the filter (e.g., "SDSS_u").
 - [`wave(f::NewType)`](@ref wave) should return the wavelength vector of the filter transmission curve with proper `Unitful.jl` units.
 - [`throughput(f::NewType)`](@ref throughput) should return the throughput vector of the filter transmission curve (no units).
 - [`detector_type(f::NewType)`](@ref detector_type) should return an instance of `PhotometricFilters.Energy` if the filter is defined for energy-counting detectors or `PhotometricFilters.Photon` for photon-counting detectors.

Additionally, all subtypes should support filter interpolation at user-defined wavelengths with a call signature `(f::NewType)(wavelengths)`. To support this, new types should implement a method like `(f::PhotometricFilter)(wave::Q) where Q <: Unitful.Length`. A generic fallback for inputs without units is already defined.
"""
abstract type AbstractFilter{T} end

# Generic methods for all AbstractFilter
# Methods to implement AbstractVector interface
Base.getindex(f::AbstractFilter, i::Int) = (wave(f)[i], throughput(f)[i])
Base.length(f::AbstractFilter) = length(throughput(f))
Base.size(f::AbstractFilter) = size(throughput(f))
Base.firstindex(f::AbstractFilter) = firstindex(throughput(f))
Base.lastindex(f::AbstractFilter) = lastindex(throughput(f))
Base.eltype(f::AbstractFilter) = Tuple{eltype(wave(f)), eltype(throughput(f))}
function Base.iterate(f::AbstractFilter, state=0)
    state == length(f) && return nothing
    return f[begin + state], state + 1
end
# Interpolation should be a generic feature of all AbstractFilter
# Concrete subtypes should implement (f::NewType)(wave::Q) where Q <: Unitful.Length
(f::AbstractFilter)(wave) = @. f(wave * wave_unit)
# Concrete subtypes should implement name(::NewType)::String, detector_type(::NewType)
Base.show(io::IO, f::AbstractFilter) = print(io, name(f))
function Base.show(io::IO, ::MIME"text/plain", f::T) where T <: AbstractFilter
    N = length(f)
    # println(io, "$N-element $T: ", name(f))
    s1 = split(string(T), ",")[1]
    println(io, "$N-element $(s1 * repeat("}", count("{", s1))): ", name(f))
    println(io, " min. wave.: ", min_wave(f))
    println(io, " max. wave.: ", max_wave(f))
    println(io, " effective wave.: ", effective_wavelength(f))
    println(io, " mean wave.: ", mean_wavelength(f))
    println(io, " central wave.: ", central_wavelength(f))
    println(io, " pivot wave.: ", pivot_wavelength(f))
    println(io, " eff. width: ", width(f))
    print(io,   " fwhm: ", fwhm(f))
end
"""
    wave(f::AbstractFilter)
Returns the wavelength vector of the filter transmission curve with proper `Unitful.jl` units.
```jldoctest
julia> using PhotometricFilters: SDSS_u, wave

julia> using Unitful: Quantity

julia> wave(SDSS_u()) isa Vector{<:Quantity}
true
```
"""
function wave(::AbstractFilter) end

"""
    throughput(f::AbstractFilter)
Returns the throughput vector of the filter transmission curve (no units).
```jldoctest
julia> using PhotometricFilters: SDSS_u, throughput

julia> throughput(SDSS_u()) isa Vector{<:Number}
true
```
"""
function throughput(::AbstractFilter) end

"""
    detector_type(f::AbstractFilter)
Return an instance of `PhotometricFilters.Energy` if the filter is defined for energy-counting detectors or `PhotometricFilters.Photon` for photon-counting detectors.
```jldoctest
julia> using PhotometricFilters: SDSS_u, detector_type, Photon

julia> detector_type(SDSS_u()) === Photon()
true
```
"""
function detector_type(::AbstractFilter) end

"""
    name(f::AbstractFilter)
Returns a string indicating a human-readable name for the filter (e.g., "SDSS_u").
```jldoctest
julia> using PhotometricFilters: SDSS_u, name

julia> name(SDSS_u())
"SDSS_u"
```
"""
function name(::AbstractFilter) end

# Statistics
"""
    effective_wavelength(f::AbstractFilter)

Returns the effective wavelength of the filter `f` using the Vega spectrum as a standard. Defined as

```math
\\frac{\\int \\lambda \\, T(\\lambda) \\text{Vg}(\\lambda) \\, d\\lambda}{\\int T(\\lambda) \\text{Vg}(\\lambda) \\, d\\lambda}
```
where ``T(\\lambda)`` is the filter transmission at wavelength ``\\lambda`` and ``\\text{Vg}(\\lambda)`` is the spectrum of Vega.
"""
function effective_wavelength(f::AbstractFilter)
    wvega, fvega = Vega()
    wl = uconvert.(wave_unit, wvega)
    filt = f.(wl) .* fvega
    norm = trapz(wl, filt)
    leff = trapz(wl, wl .* filt)
    return leff / norm
end

"""
    pivot_wavelength(f::AbstractFilter)

Returns the pivot wavelength of the filter `f`, defined for filters with `Energy` detector types as

```math
\\sqrt{ \\frac{\\int T(\\lambda) \\, d\\lambda}{\\int T(\\lambda) / \\lambda^2 \\, d\\lambda} }
```

For filters with `Photon` detector types, ``\\lambda \\, T(\\lambda)`` is substituted for ``T(\\lambda)`` in the above expression.

Internally integration is carried out using trapezoidal integration. It can be convenient to think of this as the "center of mass" of the filter.
"""
pivot_wavelength(f::AbstractFilter) = pivot_wavelength(f, detector_type(f))
function pivot_wavelength(f::AbstractFilter, ::Energy)
    wl = wave(f)
    y = throughput(f) ./ wl.^2
    norm = trapz(wl, throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end
function pivot_wavelength(f::AbstractFilter, ::Photon)
    wl = wave(f)
    y = throughput(f) ./ wl
    norm = trapz(wl, wl .* throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end

"""
    mean_wavelength(f::AbstractFilter)

Returns the mean wavelength of the filter `f`, defined as

```math
\\frac{\\int \\lambda \\, T(\\lambda) \\, d\\lambda}{\\int T(\\lambda) \\, d\\lambda}
```
"""
function mean_wavelength(f::AbstractFilter)
    wl = wave(f)
    norm = trapz(wl, throughput(f))
    lt = trapz(wl, wl .* throughput(f))
    return  lt / norm
end

"""
    central_wavelength(f::AbstractFilter)

Returns the central wavelength of the filter `f`, defined as the central wavelength between the two wavelengths used for the FWHM ([`fwhm`](@ref)).
"""
function central_wavelength(f::AbstractFilter)
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
    min_wave(f::AbstractFilter; level=0.01)

Returns the shortest wavelength at which the filter transmission is equal to `level * maximum(transmission)`.
"""
function min_wave(f::AbstractFilter; level=0.01)
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
    max_wave(f::AbstractFilter; level=0.01)

Returns the longest wavelength at which the filter transmission is equal to `level * maximum(transmission)`.
"""
function max_wave(f::AbstractFilter; level=0.01)
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
    fwhm(f::AbstractFilter)

Returns the difference between the furthest two wavelengths for which the filter transmission is equal to half its maximum value.
"""
function fwhm(f::AbstractFilter)
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
    width(f::AbstractFilter)

Returns the effective width of the filter, defined as the horizontal size of a rectangle with height equal to the maximum transmission of the filter such that the area of the rectangle is equal to the area under the filter transmission curve. This is calculated as

```math
\\frac{\\int T(\\lambda) \\, d\\lambda}{\\text{max}(T(\\lambda))}
```
"""
function width(f::AbstractFilter)
    norm = trapz(wave(f), throughput(f))
    return norm / maximum(throughput(f))
end

"""
    apply(f::AbstractFilter, wave, flux)

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
apply(filt::AbstractFilter, wave, flux) = apply!(filt, wave, flux, similar(flux))

"""
    apply!(f::AbstractFilter, wave, flux, out)

In-place version of [`apply`](@ref) which modifies `out`. It should have a compatible element type with `flux`.
"""
function apply!(filt::AbstractFilter, wave, flux, out)
    @. out = flux * filt(wave)
    return out
end

"""
    integrate(filt::AbstractFilter, wave, flux)
Returns the mean flux density of a spectrum (defined by wavelengths `wave` and fluxes `flux`) when integrated over the provided filter `filt`.

For photon counting detectors, this is

```math
\\overline{f_\\lambda} = \\frac{\\int_\\lambda \\lambda \\, f_\\lambda \\, T(\\lambda) \\, d\\lambda}{\\int_\\lambda \\lambda \\, T(\\lambda) \\, d\\lambda}
```

which can also be interpreted as the mean photon rate density, while for energy counting detectors, this is

```math
\\overline{f_\\lambda} = \\frac{\\int_\\lambda f_\\lambda \\, T(\\lambda) \\, d\\lambda}{\\int_\\lambda T(\\lambda) \\, d\\lambda}
```

which is essentially just the mean flux weighted by the filter throughput.
"""
function integrate(wave, flux, throughput, ::Energy)
    return trapz(wave, flux .* throughput) / trapz(wave, throughput)
end
function integrate(wave, flux, throughput, ::Photon)
    t1 = wave .* throughput
    return trapz(wave, flux .* t1) / trapz(wave, t1)
end
function integrate(filt::AbstractFilter, wave, flux)
    return integrate(wave, flux, filt.(wave), detector_type(filt))
end
@derived_dimension SpectralFluxDensity Unitful.𝐌 / Unitful.𝐋 / Unitful.𝐓^3
@derived_dimension SpectralEnergyDensity Unitful.𝐌 / Unitful.𝐓^2
"""
    F_nu(F_lambda::SpectralFluxDensity, λpivot)
    F_nu(F_lambda::SpectralFluxDensity, f::AbstractFilter)
Convert a spectral flux density `F_lambda` into a spectral energy density. Assuming `F_lambda` in *erg / s / cm^2 / Angstrom*, and `F_nu` in *Jy*, this conversion is

```math
F_\\nu = \\frac{10^5}{10^{-8} \\, c} \\, \\lambda^2_p \\, F_\\lambda
```

where ``c`` is the speed of light in *m/s* and ``\\lambda_p`` is the pivot wavelength (`λpivot`) in *Angstroms*. If providing an [`AbstractFilter`](@ref PhotometricFilters.AbstractFilter) as the second argument, the pivot wavelength will be automatically computed with [`pivot_wavelength`](@ref PhotometricFilters.pivot_wavelength).
"""
F_nu(Fλ::SpectralFluxDensity, λpivot) = 25370985150//760603 * _ustrip(u"angstrom", λpivot)^2 * ustrip(u"erg/s/cm^2/angstrom", Fλ) * u"Jy" # Prefactor is 1e5 / (c / 10^8), c = 2.99792458e8 m s^-1, see PR #22
F_nu(Fλ::SpectralFluxDensity, f::AbstractFilter) = F_nu(Fλ, pivot_wavelength(f))

"""
    F_lambda(F_nu::SpectralEnergyDensity, λpivot)
    F_lambda(F_nu::SpectralEnergyDensity, f::AbstractFilter)
Convert a spectral energy density `F_nu` into a spectral flux density. Assuming `F_nu` in *Jy* and `F_lambda` in *erg / s / cm^2 / Angstrom*, this conversion is

```math
F_\\lambda = \\frac{10^{-8} \\, c}{10^5} \\, \\lambda^{-2}_p \\, F_\\nu
```

where ``c`` is the speed of light in *m/s* and ``\\lambda_p`` is the pivot wavelength (`λpivot`) in *Angstroms*. If providing an [`AbstractFilter`](@ref PhotometricFilters.AbstractFilter) as the second argument, the pivot wavelength will be automatically computed with [`pivot_wavelength`](@ref PhotometricFilters.pivot_wavelength).
"""
F_lambda(F_nu::SpectralEnergyDensity, λpivot) = 760603//25370985150 / _ustrip(u"angstrom", λpivot)^2 * ustrip(u"Jy", F_nu) * u"erg/s/cm^2/angstrom" # Prefactor is (c / 10^8) / 1e5, c = 2.99792458e8 m s^-1, see PR #22
F_lambda(F_nu::SpectralEnergyDensity, f::AbstractFilter) = F_lambda(F_nu, pivot_wavelength(f))

############################################################
# Definition and methods for PhotometricFilter concrete type

struct PhotometricFilter{T, WT, DT <: DetectorType, TT <: AbstractVector{T},
                         ST <: Union{String, Nothing}, ET} <: AbstractFilter{T}
    wave::WT
    throughput::TT
    detector::DT
    name::ST
    etp::ET
end

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

julia> f[10] # Indexing into the filter as `f[i]` returns `(wave(f)[i], throughput(f)[i])`
(1009 Å, 0.25)

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

name(f::PhotometricFilter) = f.name
detector_type(f::PhotometricFilter) = f.detector
wave(f::PhotometricFilter) = f.wave
throughput(f::PhotometricFilter) = f.throughput
function (f::PhotometricFilter)(wave::Q) where Q <: Unitful.Length
    wl = uconvert.(wave_unit, wave)
    return f.etp(wl)
end

function Base.:*(f1::PhotometricFilter, f2::PhotometricFilter)
    if detector_type(f1) !== detector_type(f2)
        throw(ArgumentError("When multiplying two `PhotometricFilter` instances, they must have identical detector types -- filter $f1 has `detector_type` $(detector_type(f1)), while filter $f2 has `detector_type` $(detector_type(f2))."))
    end

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
    return PhotometricFilter(wl, through; detector=detector_type(f1), name=name)
end
