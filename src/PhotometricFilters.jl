module PhotometricFilters

using DataDeps
using DelimitedFiles

export PhotometricFilter, 
       wave, 
       throughput, 
       apply, 
       apply!

include("core.jl") # types and base utilities (like size, length)
include("library.jl") # the library of filters
include("plots.jl") # plotting recipes
include("vega.jl") # vega calibration spectrum

function __init__()
    register(PYPHOT_DATADEP)
    register(VEGA_DATADEP)
end

end
