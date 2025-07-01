import DataFrames: DataFrame
import HTTP
using OrderedCollections: OrderedDict
using Unitful: uparse, ustrip, NoUnits
import UnitfulAstro
using XML: Node, LazyNode, children, simple_value, attributes, tag, next
import VOTables

const svo_url = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
const detector_type = (Energy(), Photon()) # SVO returns 0 or 1 for energy or photon
const TYPE_VO_TO_JL = VOTables.TYPE_VO_TO_JL

"""
    get_filter(filtername::AbstractString, magsys::Symbol=:Vega)

Query the online [SVO filter service](http://svo2.cab.inta-csic.es/theory/fps) for data on a photometric filter.

# Arguments
 - `filtername::AbstractString`: The desired filter ID, in the correct SVO specification (e.g., `"2MASS/2MASS.J"`).
 - `magsys::Symbol`: Desired magnitude system for associated metadata (e.g., `"ZeroPoint"`). Can be any of `(:AB, :Vega, :ST)`. SVO uses Vega by default, so we mirror that choice here.

 # Returns
 A length-2 tuple, with elements
  1. a [`PhotometricFilter`](@ref) containing the transmission data of the filter.
  2. a dictionary containing additional metadata provided by SVO.

# Examples
```jldoctest
julia> using PhotometricFilters: get_filter

julia> filt = get_filter("2MASS/2MASS.J", :Vega);

julia> filt[1]
107-element PhotometricFilter{Float64}: 2MASS/2MASS.J
 min. wave.: 10806.470589792389 Å
 max. wave.: 14067.974683578484 Å
 effective wave.: 12285.654731403807 Å
 mean wave.: 12410.5170694321 Å
 central wave.: 12390.584132888223 Å
 pivot wave.: 12358.089456559974 Å
 eff. width: 1624.3245065600008 Å
 fwhm: 2149.1445403830403 Å

julia> filt[2] isa AbstractDict
true

julia> filt[2]["ZeroPoint"]
1594.0 Jy
```
"""
function get_filter(filtername::AbstractString, magsys::Symbol=:Vega)
    # PhotometricFilter constructor requires String, not AbstractString, so convert
    if !in(magsys, (:Vega, :AB, :ST))
        throw(ArgumentError("Valid `magsys` arguments to `get_filter` are `:Vega, :AB, :ST`, you provided $magsys."))
    end
    filtername = String(filtername)
    response = HTTP.get(svo_url; query = Dict("PhotCalID" => "$filtername/$magsys"))
    if response.status != 200 # If status is not normal,
        @info "HTTP request to SVO returned with status code $(response.status)."
    end
    
    # Parse metadata
    xml = read(IOBuffer(response.body), Node)
    info_resource = children(children(xml)[2])
    info = info_resource[1]

    if attributes(info)["value"] == "ERROR"
        errval = simple_value(children(info)[1])
        if errval == "Filter not found:"
            errval *= " $id"
        end
        error(errval)
    end

    resource = info_resource[2]
    param_nodes = children(children(resource)[1])
    d = OrderedDict{String,Any}()

    for n in param_nodes
        tag(n) == "PARAM" || continue

        key, val = _get_filter_param_node_keyval(n)
        d[key] = val

        if !isempty(children(n))
            d[key * "?"] = simple_value(children(n)[1])
        end
    end

    # Construct PhotometricFilter
    table = VOTables.read(IOBuffer(response.body); unitful=true)
    result = PhotometricFilter(table.Wavelength,
                               Vector(table.Transmission);
                               detector=detector_type[parse(Int, d["DetectorType"]) + 1],
                               name=filtername)
    return result, d
end

"""
    get_metadata()

Returns a table of the available parameters that can be used to query the SVO filter service
from the `FORMAT=metadata` VOTable they provide.

The table is a `DataFrame` from the [`DataFrames`](https://github.com/JuliaData/DataFrames.jl)
package with the following columns:
- `parameter`: parameter name that can be used for queries using [`query_filters`](@ref)
- `unit`: [`Unitful`](https://github.com/PainterQubits/Unitful.jl) unit of the parameter
- `datatype`: `Type` of the parameter
- `description`: description of the parameter
- `values`: vector of the possible values that the respective parameter can take on (e.g. for Instrument),
  or a vector of the minimum and maximum values that the parameter can assume (e.g. for WavelengthEff)

# Example
```jldoctest
julia> using DataFrames: DataFrame

julia> df = PhotometricFilters.get_metadata();

julia> df isa DataFrame
true

julia> facilities = df[findfirst(==("Facility"), df.parameter), :].values;

julia> facilities isa Vector{String}
true
```

This is not exported.
"""
function get_metadata()
    response = HTTP.get(svo_url; query = Dict("FORMAT" => "metadata"))
    xmldoc = read(IOBuffer(response.body), Node)
    param_nodes = children(children(children(xmldoc)[2])[1])
    table = [_get_metadata_param_node_table_row(n) for n in param_nodes if tag(n) == "PARAM" && !endswith(attributes(n)["name"], "_max")]
    return DataFrame(table)
end

function _get_metadata_param_node_table_row(n::Node)
    att = attributes(n)
    if endswith(att["name"], r"_min|_max")
        parameter = att["name"][7:end-4]
    else
        parameter = att["name"][7:end]
    end
    u = _get_unit(att["unit"])
    datatype = TYPE_VO_TO_JL[att["datatype"]]
    description = simple_value(children(n)[1])
    values = _get_values(children(n)[2], datatype)
    return (; parameter, unit=u, datatype, description, values)
end

function _get_unit(str::AbstractString)
    # Assume any integers are powers -> add power sign
    if !contains(str, "^")
        str = replace(str, r"(-?\d)" => x -> "^$x")
    end

    # Angstrom are not always written out
    if contains(str, "A") && !contains(str, "Angstrom")
        str = replace(str, "A" => "angstrom")
    end

    if str == "Angstrom"
        return u"angstrom"
    elseif str == ""
        return NoUnits
    else
        return uparse(str; unit_context=[Unitful, UnitfulAstro])
    end
end

function _get_values(n::Node, T)
    @assert tag(n) == "VALUES"
    c = children(n)
    if length(c) == 2 && tag(c[1]) == "MIN"
        return [parse(T, attributes(c[1])["value"]), parse(T, attributes(c[2])["value"])]
    elseif !isempty(c) && tag(c[1]) == "OPTION"
        return [attributes(valnodes)["value"] for valnodes in c]
    else
        return T[]
    end
end

_ustrip(_, val::Number) = val
_ustrip(u, val::Quantity) = ustrip(u, val)

"""
    query_filters(; queries...)

Queries the filters available from the SVO filter service with search parameters and returns a table of the filters found.

The available search parameters can be found with [`PhotometricFilters.get_metadata`](@ref). The following should be available
in general:
- `WavelengthRef`: Tuple of Numbers
- `WavelengthMean`: Tuple of Numbers
- `WavelengthEff`: Tuple of Numbers
- `WavelengthMin`: Tuple of Numbers
- `WavelengthMax`: Tuple of Numbers
- `WidthEff`: Tuple of Numbers
- `FWHM`: Tuple of Numbers
- `Instrument`: String
- `Facility`: String
- `PhotSystem`: String

The returned table is a `DataFrame` from the [`DataFrames`](https://github.com/JuliaData/DataFrames.jl)
package with all the columns of the response VOTable.

The filter information and transmission data can be obtained by calling [`get_filter`](@ref) with the ID obtained from the `filterID` column.

# Examples
```jldoctest
julia> using DataFrames: DataFrame

julia> df = query_filters(; Facility="SLOAN", WavelengthEff=(1000, 5000));

julia> df isa DataFrame
true

julia> id = df.filterID[3]
"SLOAN/SDSS.g"

julia> get_filter(id)[1] isa PhotometricFilter
true
```

```julia
# Other examples for querying
query_filters(; Facility="SLOAN") # all filters from a given facility
query_filters(; Instrument="BUSCA", WavelengthEff=(1000u"angstrom", 5000u"angstrom")) # Unitful wavelengths
```
"""
function query_filters(; kwargs...)
    isempty(kwargs) && error("At least one filter criterion must be passed to search for filters! The available parameters can be found with `PhotometricFilters.get_metadata()`.")

    # If parameters with non-Angstrom units show up in the future, the hard-coded part here has to be modified
    query = Dict(k[1] => k[2] isa AbstractString ? k[2] : join(_ustrip.(u"angstrom", k[2]), "/") for k in kwargs)
    response = HTTP.get(svo_url; query = query)

    err_str = try
        table = VOTables.read(IOBuffer(response.body))
        df = DataFrame(table)

        # Add units to fields like ZeroPoint
        n = read(IOBuffer(response.body), LazyNode)
        # Find fields
        while tag(n) != "FIELD"
            n = next(n)
        end

        while tag(n) == "FIELD"
            if haskey(attributes(n), "unit")
                field = attributes(n)["name"]
                u = _get_unit(attributes(n)["unit"])
                if hasproperty(df, field) && u != NoUnits
                    df[!, field] = df[!, field] .* u
                end
            end
            n = next(n)
        end

        return df

    catch e
        xml = read(IOBuffer(response.body), Node)
        n = children(children(xml)[2])[1]

        if attributes(n)["value"] == "ERROR"
            simple_value(children(n)[1])
        else
            throw(e)
        end
    end

    error(err_str)
end

function _get_filter_param_node_keyval(n::Node)
    att = attributes(n)
    key = att["name"]
    datatype = TYPE_VO_TO_JL[att["datatype"]]
    val = att["value"]

    if datatype <: Number
        val = parse(datatype, val)
    end

    if haskey(att, "unit")
        u =  _get_unit(att["unit"])
        val *= u
    end

    return key, val
end
