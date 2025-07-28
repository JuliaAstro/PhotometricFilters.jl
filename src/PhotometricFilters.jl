module PhotometricFilters

using DataDeps: register, DataDep, @datadep_str

export PhotometricFilter,
       get_filter,
       query_filters,
       wavelength, 
       throughput,
       detector_type,
       filtername,
       apply, 
       apply!,
       mean_flux_density,
       F_nu,
       F_lambda,
       AB_flux_zeropoint,
       AB_mag_zeropoint,
       AB_Jy_zeropoint,
       ST_flux_zeropoint,
       ST_mag_zeropoint,
       ST_Jy_zeropoint,
       Vega_flux_zeropoint,
       Vega_mag_zeropoint,
       Vega_Jy_zeropoint

include("core.jl") # types and base utilities (like size, length)
include("library.jl") # the library of filters
include("plots.jl") # plotting recipes
include("vega.jl") # vega calibration spectrum
include("svo.jl")  # Query SVO service for filter curves

function __init__()
    register(PYPHOT_DATADEP)
    register(VEGA_DATADEP)
end

end
