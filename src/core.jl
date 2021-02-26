
abstract type AbstractFilter end

abstract type DetectorType end
struct Photon <: DetectorType end
struct Energy <: DetectorType end

function Base.parse(::Type{DetectorType}, s::AbstractString)
    tok = lowercase(s)
    if tok === "photon"
        return Photon()
    elseif tok === "energy"
        return Energy()
    else
        throw(ArgumentError("failed to parse DetectorType from $s"))
    end
end

struct PhotometricFilter{WLT,TT,DT<:DetectorType,ST<:Union{String,Nothing}} <: AbstractFilter
    wave::WLT
    throughput::TT
    detector::DT
    name::ST
end

Base.show(io::IO, f::PhotometricFilter) = print(io, f.name)

function Base.show(io::IO, ::MIME"text/plain", f::PhotometricFilter)
    print(io, "PhotometricFilter: $(f.name)\n wave: ", f.wave, "\n throughput: ", f.throughput)
end

wave(f::PhotometricFilter) = f.wave
throughput(f::PhotometricFilter) = f.throughput

Base.size(f::PhotometricFilter) = size(throughput(f))
Base.length(f::PhotometricFilter) = prod(size(f))

Base.iterate(f::PhotometricFilter, args...) = iterate(wave(f), args...), iterate(throughput(f), args...)
