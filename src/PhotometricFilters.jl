module PhotometricFilters

using DataDeps
using DelimitedFiles

export PhotometricFilter,
       get_filter,
       wave, 
       throughput, 
       apply, 
       apply!

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
