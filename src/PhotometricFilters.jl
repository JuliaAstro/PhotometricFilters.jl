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

function __init__()
    register(PYPHOT_DATADEP)
end

end
