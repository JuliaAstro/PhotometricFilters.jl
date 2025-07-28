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

 - [`filtername(f::NewType)`](@ref filtername) should return a string indicating a human-readable name for the filter (e.g., "SDSS_u").
 - [`wavelength(f::NewType)`](@ref wavelength) should return the wavelength vector of the filter transmission curve with proper `Unitful.jl` units.
 - [`throughput(f::NewType)`](@ref throughput) should return the throughput vector of the filter transmission curve (no units).
 - [`detector_type(f::NewType)`](@ref detector_type) should return an instance of `PhotometricFilters.Energy` if the filter is defined for energy-counting detectors or `PhotometricFilters.Photon` for photon-counting detectors.

Additionally, all subtypes should support filter interpolation at user-defined wavelengths with a call signature `(f::NewType)(wavelengths)`. To support this, new types should implement a method like `(f::PhotometricFilter)(wavelength::Q) where Q <: Unitful.Length`. A generic fallback for inputs without units is already defined.
"""
abstract type AbstractFilter{T} end

# Generic methods for all AbstractFilter
# Methods to implement AbstractVector interface
Base.getindex(f::AbstractFilter, i::Int) = (wavelength(f)[i], throughput(f)[i])
Base.length(f::AbstractFilter) = length(throughput(f))
Base.size(f::AbstractFilter) = size(throughput(f))
Base.firstindex(f::AbstractFilter) = firstindex(throughput(f))
Base.lastindex(f::AbstractFilter) = lastindex(throughput(f))
Base.eltype(f::AbstractFilter) = Tuple{eltype(wavelength(f)), eltype(throughput(f))}
function Base.iterate(f::AbstractFilter, state=0)
    state == length(f) && return nothing
    return f[begin + state], state + 1
end
Base.:(==)(f1::AbstractFilter, f2::AbstractFilter) = wavelength(f1) == wavelength(f2) && throughput(f1) == throughput(f2)
# Interpolation should be a generic feature of all AbstractFilter
# Concrete subtypes should implement (f::NewType)(wavelength::Q) where Q <: Unitful.Length
(f::AbstractFilter)(wavelength) = @. f(wavelength * wave_unit)
# Concrete subtypes should implement filtername(::NewType)::String, detector_type(::NewType)
Base.show(io::IO, f::AbstractFilter) = print(io, filtername(f))
function Base.show(io::IO, ::MIME"text/plain", f::T) where T <: AbstractFilter
    N = length(f)
    # println(io, "$N-element $T: ", filtername(f))
    s1 = split(string(T), ",")[1]
    println(io, "$N-element $(s1 * repeat("}", count("{", s1))): ", filtername(f))
    println(io, " reference wave.: ", reference_wavelength(f))
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
    wavelength(f::AbstractFilter)
Returns the wavelength vector of the filter transmission curve with proper `Unitful.jl` units.
```jldoctest
julia> using PhotometricFilters: SDSS_u, wavelength

julia> using Unitful: Quantity

julia> wavelength(SDSS_u()) isa Vector{<:Quantity}
true
```
"""
function wavelength(::AbstractFilter) end

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
    filtername(f::AbstractFilter)
Returns a string indicating a human-readable name for the filter (e.g., "SDSS_u").
```jldoctest
julia> using PhotometricFilters: SDSS_u, filtername

julia> filtername(SDSS_u())
"SDSS_u"
```
"""
function filtername(::AbstractFilter) end

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
    wl = wavelength(f)
    y = throughput(f) ./ wl.^2
    norm = trapz(wl, throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end
function pivot_wavelength(f::AbstractFilter, ::Photon)
    wl = wavelength(f)
    y = throughput(f) ./ wl
    norm = trapz(wl, wl .* throughput(f))
    lp2 = norm / trapz(wl, y)
    return sqrt(lp2)
end

"""
    reference_wavelength(f::AbstractFilter)

Returns the reference wavelength of the filter `f`, used for conversions of the flux and for determination of magnitudes.

By default the pivot wavelength is returned ([`pivot_wavelength`](@ref)), but filter providers sometimes provide their own specified values.
"""
reference_wavelength(f::AbstractFilter) = pivot_wavelength(f)

"""
    mean_wavelength(f::AbstractFilter)

Returns the mean wavelength of the filter `f`, defined as

```math
\\frac{\\int \\lambda \\, T(\\lambda) \\, d\\lambda}{\\int T(\\lambda) \\, d\\lambda}
```
"""
function mean_wavelength(f::AbstractFilter)
    wl = wavelength(f)
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
    wl = wavelength(f)
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
    wl = wavelength(f)
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
    wl = wavelength(f)
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
    wl = wavelength(f)
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
    norm = trapz(wavelength(f), throughput(f))
    return norm / maximum(throughput(f))
end

"""
    apply_throughput(f::AbstractFilter, wavelengths, flux)

Use linear interpolation to map the wavelengths of the photometric filter `f` to the given `wavelengths` and apply the filter throughput to the `flux`. The provided `wavelengths` and those of the filter must be compatible. This means if one has units, the other one needs units, too.

```jldoctest
julia> using PhotometricFilters: SDSS_u, wave_unit, apply_throughput

julia> f = SDSS_u();

julia> λ = 3000:4000
3000:4000

julia> flux = fill(1.0, length(λ)); # If `flux` is all `1`, `apply_throughput` reduces to `f` interpolated at `λ`

julia> apply_throughput(f, λ, flux) == f(λ)
true

julia> λ_u = λ .* wave_unit # Can also put units on λ
(3000:4000) Å

julia> apply_throughput(f, λ_u, flux) == f.(λ_u)
true
```
"""
apply_throughput(filt::AbstractFilter, wavelengths, flux) = apply_throughput!(filt, wavelengths, flux, similar(flux))

"""
    apply_throughput!(f::AbstractFilter, wavelengths, flux, out)

In-place version of [`apply_throughput`](@ref) which modifies `out`. It should have a compatible element type with `flux`.
"""
function apply_throughput!(filt::AbstractFilter, wavelengths, flux, out)
    @. out = flux * filt(wavelengths)
    return out
end

"""
    mean_flux_density(filt::AbstractFilter, wavelengths, flux)
Returns the mean flux density of a spectrum (defined by wavelengths `wavelengths` and fluxes `flux`) when integrated over the provided filter `filt`.

For photon counting detectors, this is

```math
\\overline{f_\\lambda} = \\frac{\\int_\\lambda \\lambda \\, f_\\lambda \\, T(\\lambda) \\, d\\lambda}{\\int_\\lambda \\lambda \\, T(\\lambda) \\, d\\lambda}
```

which can also be interpreted as the mean photon rate density, while for energy counting detectors, this is

```math
\\overline{f_\\lambda} = \\frac{\\int_\\lambda f_\\lambda \\, T(\\lambda) \\, d\\lambda}{\\int_\\lambda T(\\lambda) \\, d\\lambda}
```

which is essentially just the mean flux weighted by the filter throughput.

Below we show example usage that can be compared against [this example](https://github.com/mfouesneau/pyphot/blob/master/examples/Sun_Vega.ipynb) from pyphot.

```jldoctest
julia> using PhotometricFilters: mean_flux_density, HST_WFC3_F110W, Vega

julia> using Unitful, UnitfulAstro

julia> mfd = mean_flux_density(HST_WFC3_F110W(), Vega()...);

julia> isapprox(mfd, 4.082289e-10 * u"erg/s/cm^2/angstrom"; rtol=1e-3)
true
```
"""
function mean_flux_density(wavelengths, flux, throughput, ::Energy)
    return trapz(wavelengths, flux .* throughput) / trapz(wavelengths, throughput)
end
function mean_flux_density(wavelengths, flux, throughput, ::Photon)
    t1 = wavelengths .* throughput
    return trapz(wavelengths, flux .* t1) / trapz(wavelengths, t1)
end
function mean_flux_density(filt::AbstractFilter, wavelengths, flux)
    return mean_flux_density(wavelengths, flux, filt.(wavelengths), detector_type(filt))
end
@derived_dimension SpectralFluxDensity Unitful.𝐌 / Unitful.𝐋 / Unitful.𝐓^3
@derived_dimension SpectralEnergyDensity Unitful.𝐌 / Unitful.𝐓^2
"""
    F_nu(F_lambda::SpectralFluxDensity, λref)
    F_nu(F_lambda::SpectralFluxDensity, f::AbstractFilter)
Convert a spectral flux density `F_lambda` into a spectral energy density. Assuming `F_lambda` in *erg / s / cm^2 / Angstrom*, and `F_nu` in *Jy*, this conversion is

```math
F_\\nu = \\frac{10^5}{10^{-8} \\, c} \\, \\lambda^2_r \\, F_\\lambda
```

where ``c`` is the speed of light in *m/s* and ``\\lambda_r`` is the reference wavelength (`λref`) in *Angstroms*. If providing an [`AbstractFilter`](@ref PhotometricFilters.AbstractFilter) as the second argument, the reference wavelength will be automatically computed with [`reference_wavelength`](@ref PhotometricFilters.reference_wavelength).
"""
F_nu(Fλ::SpectralFluxDensity, λref) = 25370985150//760603 * _ustrip(u"angstrom", λref)^2 * ustrip(u"erg/s/cm^2/angstrom", Fλ) * u"Jy" # Prefactor is 1e5 / (c / 10^8), c = 2.99792458e8 m s^-1, see PR #22
F_nu(Fλ::SpectralFluxDensity, f::AbstractFilter) = F_nu(Fλ, reference_wavelength(f))
F_nu(Fν::SpectralEnergyDensity, ::AbstractFilter) = uconvert(u"Jy", Fν)

"""
    F_lambda(F_nu::SpectralEnergyDensity, λref)
    F_lambda(F_nu::SpectralEnergyDensity, f::AbstractFilter)
Convert a spectral energy density `F_nu` into a spectral flux density. Assuming `F_nu` in *Jy* and `F_lambda` in *erg / s / cm^2 / Angstrom*, this conversion is

```math
F_\\lambda = \\frac{10^{-8} \\, c}{10^5} \\, \\lambda^{-2}_r \\, F_\\nu
```

where ``c`` is the speed of light in *m/s* and ``\\lambda_r`` is the reference wavelength (`λref`) in *Angstroms*. If providing an [`AbstractFilter`](@ref PhotometricFilters.AbstractFilter) as the second argument, the reference wavelength will be automatically computed with [`reference_wavelength`](@ref PhotometricFilters.reference_wavelength).
"""
F_lambda(F_nu::SpectralEnergyDensity, λref) = 760603//25370985150 / _ustrip(u"angstrom", λref)^2 * ustrip(u"Jy", F_nu) * u"erg/s/cm^2/angstrom" # Prefactor is (c / 10^8) / 1e5, c = 2.99792458e8 m s^-1, see PR #22
F_lambda(F_nu::SpectralEnergyDensity, f::AbstractFilter) = F_lambda(F_nu, reference_wavelength(f))
F_lambda(Fλ::SpectralFluxDensity, ::AbstractFilter) = uconvert(u"erg/s/cm^2/angstrom", Fλ)

"""
    Vega_zeropoint_flux(f::AbstractFilter)
Returns the Vega flux zero point of the filter `f` in units of erg / s / cm^2 / Angstrom.

```jldoctest
julia> using PhotometricFilters: Vega_zeropoint_flux, HST_WFC3_F110W

julia> using Unitful

julia> isapprox(Vega_zeropoint_flux(HST_WFC3_F110W()), 4.082289e-10 * u"erg/s/cm^2/angstrom"; rtol=1e-3)
true
```
"""
Vega_zeropoint_flux(f::AbstractFilter) = mean_flux_density(f, Vega()...)
"""
    Vega_zeropoint_mag(f::AbstractFilter)
Returns the Vega magnitude zero point of the filter `f`. Can be used to calculate Vega magnitudes from spectra in units of [`F_lambda`](@ref) as
```math
m_{\\text{Vega}} = -2.5 * \\text{log} \\left( f_\\lambda \\right) - \\text{Zpt}
```

```jldoctest
julia> using PhotometricFilters: Vega_zeropoint_mag, HST_WFC3_F110W

julia> isapprox(Vega_zeropoint_mag(HST_WFC3_F110W()), 23.4727487; rtol=1e-3)
true
```
"""
function Vega_zeropoint_mag(f::AbstractFilter)
    flux = Vega_zeropoint_flux(f)
    return -25//10 * log10(ustrip(flux))
end

"""
    Vega_zeropoint_Jy(f::AbstractFilter)
Returns the Vega flux zeropoint of the filter `f` in Jansky.

```jldoctest
julia> using PhotometricFilters: Vega_zeropoint_Jy, HST_WFC3_F110W

julia> using Unitful, UnitfulAstro

julia> isapprox(Vega_zeropoint_Jy(HST_WFC3_F110W()), 1816.43597 * u"Jy"; rtol=1e-3)
true
```
"""
Vega_zeropoint_Jy(f::AbstractFilter) = F_nu(Vega_zeropoint_flux(f), f)

"""
    Vega_mag(f::AbstractFilter, wavelengths, flux)
Calculates the Vega magnitude in the given filter `f` from a spectrum defined by arrays `wavelengths` and `flux`, both of which must have valid Unitful units.

By definition, this method will return 0 for Vega.

```jldoctest
julia> using PhotometricFilters: Vega_mag, Vega, HST_WFC3_F110W

julia> isapprox(Vega_mag(HST_WFC3_F110W(), Vega()...), 0; rtol=1e-3)
true
```
"""
function Vega_mag(f::AbstractFilter, wavelengths, flux)
    fbar = ustrip(u"erg/s/cm^2/angstrom", mean_flux_density(f, wavelengths, F_lambda.(flux, Ref(f))))
    return -25//10 * log10(fbar) - Vega_zeropoint_mag(f)
end

"""
    ST_zeropoint_mag(::AbstractFilter)
Returns the ST magnitude zero point, which is always equal to 21.1.
```math
m_{\\text{ST}} = -2.5 * \\text{log} \\left( f_\\lambda \\right) - 21.1
```

```jldoctest
julia> using PhotometricFilters: ST_zeropoint_mag, SDSS_u

julia> float(ST_zeropoint_mag(SDSS_u()))
21.1
```
"""
ST_zeropoint_mag(::AbstractFilter) = 211//10

"""
    ST_zeropoint_flux(f::AbstractFilter)
Returns the ST flux zero point of the filter `f` in units of erg / s / cm^2 / Angstrom.

```jldoctest
julia> using PhotometricFilters: ST_zeropoint_flux, HST_WFC3_F110W

julia> using Unitful

julia> isapprox(ST_zeropoint_flux(HST_WFC3_F110W()), 3.6307805e-9 * u"erg/s/cm^2/angstrom"; rtol=1e-3)
true
```
"""
ST_zeropoint_flux(f::AbstractFilter) = exp10(-4//10 * ST_zeropoint_mag(f)) * u"erg/s/cm^2/angstrom"

"""
    ST_zeropoint_Jy(f::AbstractFilter)
Returns the ST flux zeropoint of the filter `f` in Jansky.

```jldoctest
julia> using PhotometricFilters: ST_zeropoint_Jy, HST_WFC3_F110W

julia> using Unitful, UnitfulAstro

julia> isapprox(ST_zeropoint_Jy(HST_WFC3_F110W()), 16155.46954* u"Jy"; rtol=1e-3)
true
```
"""
ST_zeropoint_Jy(f::AbstractFilter) = F_nu(ST_zeropoint_flux(f), f)

"""
    ST_mag(f::AbstractFilter, wavelengths, flux)
Calculates the ST magnitude in the given filter `f` from a spectrum defined by arrays `wavelengths` and `flux`, both of which must have valid Unitful units.

```jldoctest
julia> using PhotometricFilters: ST_mag, Vega, HST_WFC3_F110W

julia> isapprox(ST_mag(HST_WFC3_F110W(), Vega()...), 2.372748728; rtol=1e-3)
true
```
"""
function ST_mag(f::AbstractFilter, wavelengths, flux)
    fbar = ustrip(u"erg/s/cm^2/angstrom", mean_flux_density(f, wavelengths, F_lambda.(flux, Ref(f))))
    return -25//10 * log10(fbar) - ST_zeropoint_mag(f)
end

"""
    AB_zeropoint_Jy(::AbstractFilter{T})
Returns the AB flux zeropoint in Jansky. It is often approximated that this is 3631 Jy, following from the definition ``m_\\text{AB} = -2.5 \\text{log} f_\\nu - 48.6`` where ``f_\\nu`` is in units of erg / s / cm^2 / Hz. This can be solved for ``m_\\text{AB} = 0`` to give ``f_{\\nu, 0} = 10^{\\frac{48.6}{-2.5}}`` which is approximately ``3.631 \\times 10^{-20}`` erg / s / cm^2 / Hz, or ≈ 3631 Jy. This function returns the exact value.

```jldoctest
julia> using PhotometricFilters: AB_zeropoint_Jy, HST_WFC3_F110W

julia> using Unitful, UnitfulAstro

julia> isapprox(AB_zeropoint_Jy(HST_WFC3_F110W()), 3630.78054 * u"Jy"; rtol=1e-3)
true
```
"""
AB_zeropoint_Jy(::AbstractFilter) = exp10(48.6 / -2.5 + 23) * u"Jy"
"""
    AB_zeropoint_flux(f::AbstractFilter)
Returns the AB flux zero point of the filter `f` in units of erg / s / cm^2 / Angstrom.

```jldoctest
julia> using PhotometricFilters: AB_zeropoint_flux, HST_WFC3_F110W

julia> using Unitful

julia> isapprox(AB_zeropoint_flux(HST_WFC3_F110W()), 8.159816925e-10 * u"erg/s/cm^2/angstrom"; rtol=1e-3)
true
```
"""
AB_zeropoint_flux(f::AbstractFilter) = F_lambda(AB_zeropoint_Jy(f), f)

"""
    AB_zeropoint_mag(f::AbstractFilter)
Returns the AB magnitude zero point of the filter `f`. Can be used to calculate AB magnitudes from spectra in units of [`F_lambda`](@ref) as
```math
m_{\\text{AB}} = -2.5 * \\text{log} \\left( f_\\lambda \\right) - \\text{Zpt}
```

```jldoctest
julia> using PhotometricFilters: AB_zeropoint_mag, HST_WFC3_F110W

julia> isapprox(AB_zeropoint_mag(HST_WFC3_F110W()), 22.7207989; rtol=1e-3)
true
```
"""
AB_zeropoint_mag(f::AbstractFilter) = -25//10 * log10(ustrip(AB_zeropoint_flux(f)))

"""
    AB_mag(f::AbstractFilter, wavelengths, flux)
Calculates the AB magnitude in the given filter `f` from a spectrum defined by arrays `wavelengths` and `flux`, both of which must have valid Unitful units.

```jldoctest
julia> using PhotometricFilters: AB_mag, Vega, HST_WFC3_F110W

julia> isapprox(AB_mag(HST_WFC3_F110W(), Vega()...), 0.7519497; rtol=1e-3)
true
```
"""
function AB_mag(f::AbstractFilter, wavelengths, flux)
    fbar = ustrip(u"erg/s/cm^2/angstrom", mean_flux_density(f, wavelengths, F_lambda.(flux, Ref(f))))
    return -25//10 * log10(fbar) - AB_zeropoint_mag(f)
end

############################################################
# Definition and methods for PhotometricFilter concrete type

struct PhotometricFilter{T, WT, DT <: DetectorType, TT <: AbstractVector{T},
                         ST <: Union{String, Nothing}, ET} <: AbstractFilter{T}
    wavelength::WT
    throughput::TT
    detector::DT
    filtername::ST
    etp::ET
end

"""
    PhotometricFilter(wavelength::AbstractVector, throughput::AbstractVector{T};
                      detector::DetectorType=Photon(), filtername::Union{String, Nothing}=nothing)
Struct representing a photometric filter, defined by vectors of wavelengths (`wavelength`) and filter throughputs (`throughput`).
`wavelength` can have `Unitful` units attached, otherwise they are assumed to be $wave_unit.
Optional keyword arguments define the detector type for which the filter is valid and a name to identify the filter.
```jldoctest
julia> using PhotometricFilters: PhotometricFilter, Photon, wavelength, throughput

julia> using Unitful

julia> f = PhotometricFilter(1000:2000, vcat(fill(0.25, 250), fill(0.5, 500), fill(0.25, 251))) # Specify only wavelength and throughput
1001-element PhotometricFilter{Float64}: nothing
 reference wave.: 1478.1028279485677 Å
 min. wave.: 1000 Å
 max. wave.: 2000 Å
 effective wave.: 1603.6927025575474 Å
 mean wave.: 1499.8333333333333 Å
 central wave.: 1499.5 Å
 pivot wave.: 1478.1028279485677 Å
 eff. width: 750.0 Å
 fwhm: 501.0 Å

julia> f == PhotometricFilter(uconvert.(Unitful.nm, wavelength(f)), throughput(f)) # Can also specify wavelength argument with Unitful units
true

julia> f[10] # Indexing into the filter as `f[i]` returns `(wavelength(f)[i], throughput(f)[i])`
(1009 Å, 0.25)

julia> f(1001.1) # Calling `f` like a function interpolates the throughput
0.25

julia> f(100.11 * Unitful.nm) # Can also specify wavelength with units
0.25
```
"""
function PhotometricFilter(wavelength::AbstractVector, throughput::AbstractVector{T};
                           detector::DetectorType=Photon(), filtername::Union{String, Nothing}=nothing) where T
    if length(wavelength) != length(throughput)
        throw(ArgumentError("Wavelength and throughput arrays must have equal length"))
    end
    bc = zero(T)
    wv = _convert_wave.(wavelength)
    deduplicate_knots!(wv; move_knots=true) # Ensure wavelength entries are all unique
    etp = linear_interpolation(wv, throughput; extrapolation_bc=bc)
    return PhotometricFilter(wv, throughput, detector, filtername, etp)
end

filtername(f::PhotometricFilter) = f.filtername
detector_type(f::PhotometricFilter) = f.detector
wavelength(f::PhotometricFilter) = f.wavelength
throughput(f::PhotometricFilter) = f.throughput
function (f::PhotometricFilter)(wavelength::Q) where Q <: Unitful.Length
    wl = uconvert.(wave_unit, wavelength)
    return f.etp(wl)
end

function Base.:*(f1::PhotometricFilter, f2::PhotometricFilter)
    if detector_type(f1) !== detector_type(f2)
        throw(ArgumentError("When multiplying two `PhotometricFilter` instances, they must have identical detector types -- filter $f1 has `detector_type` $(detector_type(f1)), while filter $f2 has `detector_type` $(detector_type(f2))."))
    end

    # find total extent and log-spacing
    w1 = wavelength(f1)
    w2 = wavelength(f2)
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
    filtername = "$(f1.filtername) * $(f2.filtername)"
    return PhotometricFilter(wl, through; detector=detector_type(f1), filtername=filtername)
end
