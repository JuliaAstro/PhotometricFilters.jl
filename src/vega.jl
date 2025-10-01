import HDF5
import FITSIO

"""
    Vega(name::String = "alpha_lyr_stis_011")
Loads the reference spectrum with filename `name` and returns an appropriate instance of `Vega` that can be used to compute zeropoints and magnitudes in the Vega magnitude system.

If the provided `name` is the full path to the existing file on disk, the spectrum is loaded from that file. Otherwise, it is downloaded from the [CALSPEC database](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec) of standard stars maintained by STScI. Files downloaded this way are cached for future use. Specifically, we draw from the extended catalog [here](https://ssb.stsci.edu/cdbs/calspec/). The CALSPEC catalog with the most recent reference spectrum for each star is located [here](https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/). Standard Vega spectra start with `"alpha_lyr"`. `name`s like `"alpha_lyr_stis_XXX"` are based on composites of calibrated stellar models and HST/STIS data, while `"alpha_lyr_mod_XXX"` are based on stellar models only.

Sometimes Vega is not used as the standard star for photometric systems even when the system follows the Vega magnitude convention. For example, in the near-IR Sirius is often used as the standard reference spectrum rather than Vega. This is the case for the definition of the [JWST zeropoints](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#gsc.tab=0), which presently use the `"sirius_stis_005.fits"` CALSPEC spectrum as their standard. To load a different standard, simply provide the corresponding `name` to `Vega`. For example, to load the Sirius spectrum, use `Vega("sirius_stis_005")` (the `".fits"` extension is optional).
"""
function Vega(name::String = "alpha_lyr_stis_011")
    name = splitext(name)[1] * ".fits" # ensure .fits extension
    if isfile(name)
        return FITSIO.FITS(name, "r") do f
            wavelength = read(f[2], "WAVELENGTH") .* u"angstrom"
            flux = read(f[2], "FLUX") .* u"erg/s/cm^2/angstrom"
            Vega(wavelength, flux, basename(name))
        end
    end

    # Load from CALSPEC database or local cache
    fname = joinpath(vega_cache, name)
    if !isfile(fname)
        HTTP.download("https://ssb.stsci.edu/cdbs/calspec/" * name, fname)
    end
    return Vega(fname)
end
Vega(name::AbstractString) = Vega(string(name))

"""
    get_calspec_names([substring::AbstractString])
Returns a list of the names of the available spectral standards that can be download from CALSPEC and used as a standard in the Vega magnitude system.

If the optional `substring::AbstractString` argument is provided, then the list of names is filtered to only include those that contain the provided substring.

```jldoctest
julia> using PhotometricFilters: Vega, get_calspec_names

julia> names = get_calspec_names();

julia> names isa Vector{String}
true

julia> Vega(names[1]) isa Vega
true

julia> vega_standards = get_calspec_names("alpha_lyr");

julia> all(map(x -> occursin("alpha_lyr", x), vega_standards))
true
```
"""
function get_calspec_names()
    response = HTTP.get("https://ssb.stsci.edu/cdbs/calspec/")
    xml = read(IOBuffer(response.body), LazyNode)

    calspec_names = String[]

    for node in xml
        if tag(node) == "a"
            filename = value(children(node)[1])
            calspec_name, ext = splitext(filename)
            ext == ".fits" || continue
            push!(calspec_names, calspec_name)
        end
    end

    return calspec_names
end
get_calspec_names(f::AbstractString) = filter(x -> occursin(f, x), get_calspec_names())