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
       AB_zeropoint_flux,
       AB_zeropoint_mag,
       AB_zeropoint_Jy,
       ST_zeropoint_flux,
       ST_zeropoint_mag,
       ST_zeropoint_Jy,
       Vega_zeropoint_flux,
       Vega_zeropoint_mag,
       Vega_zeropoint_Jy

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
