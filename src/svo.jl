import HTTP
import VOTables

const svo_url = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
const detector_type = (Energy(), Photon()) # SVO returns 0 or 1 for energy or photon

function get_filter(fname::AbstractString)
    # response = HTTP.get(svo_url, ["Accept" => ""], Dict("ID" => fname))
    response = HTTP.get(svo_url; query = Dict("ID" => fname))
    if response.status != 200 # If status is not normal,
        @info "HTTP request to SVO returned with status code $(response.status)."
    end
    table = VOTables.read(IOBuffer(response.body))
    # Need to get DetectorType parameter; VOTables does not currently do parameters
    # so we just have to index into the right place.
    # SVO returns 0 or 1 for energy or photon, so we add one to get index into
    # detector_type constant above
    dtype = parse(Int, split(String(response.body), "DetectorType")[2][10]) + 1
    # I think SVO always returns wavelengths in angstroms; the votable does contain
    # a unit "param" keyword but cannot currently read it
    return PhotometricFilter(table.Wavelength * u"angstrom",
                             Vector(table.Transmission);
                             detector=detector_type[dtype],
                             name=fname)
end
