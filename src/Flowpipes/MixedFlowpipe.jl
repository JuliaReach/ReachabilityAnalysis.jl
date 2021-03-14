# ============================================
# Flowpipes not contiguous in time
# ============================================

"""
    MixedFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same time, such that they are
not necessarily contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

This type does not assume that the flowpipes are contiguous in time.
"""
struct MixedFlowpipe{N, RT<:AbstractReachSet{N}, FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::VectorOfArray{RT, 2, Vector{FT}}
    ext::Dict{Symbol, Any}
end

function MixedFlowpipe(Fk::Vector{FT}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    ext = Dict{Symbol, Any}()
    return MixedFlowpipe{N, RT, FT}(voa, ext)
end

function MixedFlowpipe(Fk::Vector{FT}, ext::Dict{Symbol, Any}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    return MixedFlowpipe{N, RT, FT}(voa, ext)
end

# interface functions
array(fp::MixedFlowpipe) = fp.Fk
flowpipe(fp::MixedFlowpipe) = fp
setrep(::Type{MixedFlowpipe{N, RT, FT}}) where {N, RT, FT} = RT
numrsets(fp::MixedFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::MixedFlowpipe, I::Int...) = getindex(fp.Fk, I...)

@inline tspan(fp::MixedFlowpipe) = error("for mixed flowpipes you should specify the index, as in `tspan(fp, i)`")
@inline tspan(fp::MixedFlowpipe, i::Int) = tspan(fp[i])

function Base.similar(fp::MixedFlowpipe{N, RT, FT}) where {N, RT, FT}
    return MixedFlowpipe(Vector{FT}())
end

function (fp::MixedFlowpipe)(t::Number)
    error("not implemented yet")
end

function (fp::MixedFlowpipe)(dt::TimeInterval)
    error("not implemented yet")
end

function overapproximate(fp::MixedFlowpipe, args...)
    return MixedFlowpipe([overapproximate(Fi, args...) for Fi in fp], fp.ext)
end

function project(fp::MixedFlowpipe, args...)
    return [project(F, args...) for F in fp]
end

# LazySets interface
function LazySets.ρ(d::AbstractVector, fp::MixedFlowpipe)
    return maximum(ρ(d, F) for F in array(fp))
end

function LazySets.σ(d::AbstractVector, fp::MixedFlowpipe)
    error("not implemented")
end
