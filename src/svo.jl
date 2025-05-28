import HTTP
import VOTables

const svo_url = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
const detector_type = (Energy(), Photon()) # SVO returns 0 or 1 for energy or photon

"""
    get_filter(filtername::AbstractString)

Query the online [SVO filter service](http://svo2.cab.inta-csic.es/theory/fps) for the photometric filter `filtername`, which must adhere to their naming scheme. Returns an instance of [`PhotometricFilter`](@ref PhotometricFilters.PhotometricFilter). 

# Examples
```jldoctest
julia> using PhotometricFilters: get_filter

julia> get_filter("2MASS/2MASS.J")
107-element PhotometricFilter{Float64}: 2MASS/2MASS.J
 min. wave.: 10820.0 Å
 max. wave.: 14060.0 Å
 effective wave.: 12285.654731403807 Å
 central wave.: 12410.5170694321 Å
 pivot wave.: 12358.089456559974 Å
 eff. width: 1624.3245065600008 Å
 fwhm: 2170.0 Å
```
"""
function get_filter(filtername::AbstractString)
    # PhotometricFilter constructor requires String, not AbstractString, so convert
    filtername = String(filtername)
    response = HTTP.get(svo_url; query = Dict("ID" => filtername))
    if response.status != 200 # If status is not normal,
        @info "HTTP request to SVO returned with status code $(response.status)."
    end
    table = VOTables.read(IOBuffer(response.body); unitful=true)
    # Need to get DetectorType parameter; VOTables does not currently do parameters
    # so we just have to index into the right place.
    # SVO returns 0 or 1 for energy or photon, so we add one to get index into
    # detector_type constant above
    dtype = parse(Int, split(String(response.body), "DetectorType")[2][10]) + 1
    return PhotometricFilter(table.Wavelength,
                             Vector(table.Transmission);
                             detector=detector_type[dtype],
                             name=filtername)
end
