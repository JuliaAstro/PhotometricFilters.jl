module PhotometricFilters

using DataDeps: register, DataDep, @datadep_str
# This will be filled in inside `__init__()`
vega_cache = ""
filter_cache = ""
using Scratch: @get_scratch!

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
       Vega, AB, ST,
       zeropoint_Jy,
       zeropoint_flux,
       zeropoint_mag,
       magnitude

include("core.jl") # types and base utilities (like size, length)
include("library.jl") # the library of filters
include("plots.jl") # plotting recipes
include("vega.jl") # vega calibration spectrum
include("svo.jl")  # Query SVO service for filter curves

function __init__()
    register(PYPHOT_DATADEP)
    global vega_cache = @get_scratch!("vega_standards")
    global filter_cache = @get_scratch!("filter_cache")
end

end
