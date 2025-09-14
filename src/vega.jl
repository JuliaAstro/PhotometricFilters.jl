import HDF5
import FITSIO
import FITSFiles

const calspec_url = "https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec"

const alpha_lyr_stis_011 = DataDep(
    "alpha_lyr_stis_011",
    "Reference spectrum for Vega (α Lyrae) sourced from STScI's [CALSPEC database]($calspec_url).",
    "https://rawcdn.githack.com/mfouesneau/pyphot/943b2a5e35c88b4419d5e93a613e44f22994eae3/pyphot/libs/vega.hd5",
    "b178445de4d079829b2806b13d203f8e197618d773989badbefc16c87a13990c"
)
# const VEGA_DATADEP = DataDep(
#     "vega",
#     "Vega spectrum for calibrating filters",
#     "https://rawcdn.githack.com/mfouesneau/pyphot/943b2a5e35c88b4419d5e93a613e44f22994eae3/pyphot/libs/vega.hd5",
#     "b178445de4d079829b2806b13d203f8e197618d773989badbefc16c87a13990c"
# )

"""
    Vega(name::String = "alpha_lyr_stis_011")
Loads the reference spectrum with filename `name` and returns it as an instance of `Vega`.

If the provided `name` is the full path to the existing file on disk, the spectrum is loaded from that file. Otherwise, it is downloaded from the [CALSPEC database]($calspec_url) of standard stars maintained by STScI and cached for future use. Specifically, we draw from the extended catalog [here](https://ssb.stsci.edu/cdbs/calspec/). Standard Vega spectra start with `"alpha_lyr"`. `name`s like `"alpha_lyr_stis_XXX"` are based on composites of calibrated stellar models and HST/STIS data, while `"alpha_lyr_mod_XXX"` are based on stellar models only.

Sometimes Vega is not used as the standard star for photometric systems even when the system follows the Vega magnitude convention. For example, in the near-IR Sirius is often used as the standard reference spectrum rather than Vega. This is the case for the definition of the [JWST zeropoints](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#gsc.tab=0), which presently use the `"sirius_stis_005.fits"` CALSPEC spectrum as their standard. To load a different standard, simply provide the corresponding `name` to `Vega`. For example, to load the Sirius spectrum, use `Vega("sirius_stis_005")`.
"""
function Vega(name::String = "alpha_lyr_stis_011")
    name = splitext(name)[1] * ".fits" # ensure .fits extension
    if isfile(name)
        try
            return FITSIO.FITS(name, "r") do f
                wavelength = read(f[2], "WAVELENGTH") .* u"angstrom"
                flux = read(f[2], "FLUX") .* u"erg/s/cm^2/angstrom"
                Vega(wavelength, flux, name)
            end
        catch e
            println("Failed to load Vega spectrum from local file $name with error")
            # showerror(stdout, e)  # Displays the error message
            throw(e)
        end
    end

    # Load from CALSPEC database or local cache
    fname = joinpath(vega_cache, name)
    if !isfile(fname)
        try
            download("https://ssb.stsci.edu/cdbs/calspec/" * name, fname)
        catch e
            println("Failed to download Vega spectrum $name from CALSPEC database with error")
            throw(e)
            # showerror(stdout, e)  # Displays the error message
        end
    end
    return Vega(fname)
end
Vega(name::AbstractString) = Vega(string(name))

# function Vega()
#     local wavelength, flux
#     HDF5.h5open(datadep"vega/vega.hd5") do fh
#         node = fh["spectrum"]
#         data = read(node)
#         wlunit = parse_unit(read(HDF5.attributes(node)["WAVELENGTH_UNIT"]))
#         wavelength = map(i -> uconvert(wave_unit, i.WAVELENGTH * wlunit), data)
#         fluxunit = u"erg/s/cm^2/angstrom"
#         flux = map(i -> i.FLUX * fluxunit, data)
#     end
#     return Vega(wavelength, flux, "vega.hdf5")
# end
