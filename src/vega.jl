import HDF5

const VEGA_DATADEP = DataDep(
    "vega",
    "Vega spectrum for calibrating filters",
    "https://rawcdn.githack.com/mfouesneau/pyphot/943b2a5e35c88b4419d5e93a613e44f22994eae3/pyphot/libs/vega.hd5",
    "b178445de4d079829b2806b13d203f8e197618d773989badbefc16c87a13990c"
)

function Vega()
    local wavelength, flux
    HDF5.h5open(datadep"vega/vega.hd5") do fh
        node = fh["spectrum"]
        data = read(node)
        wlunit = parse_unit(read(HDF5.attributes(node)["WAVELENGTH_UNIT"]))
        wavelength = map(i -> uconvert(wave_unit, i.WAVELENGTH * wlunit), data)
        fluxunit = u"erg/s/cm^2/angstrom"
        flux = map(i -> i.FLUX * fluxunit, data)
    end
    return wavelength, flux
end
