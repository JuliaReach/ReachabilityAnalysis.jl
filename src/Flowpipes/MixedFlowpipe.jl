# ============================================
# Flowpipes not contiguous in time
# ============================================

"""
    MixedFlowpipe{N,RT<:AbstractReachSet{N},FT<:AbstractFlowpipe} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same time, such that they are
not necessarily contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

This type does not assume that the flowpipes are contiguous in time.
"""
struct MixedFlowpipe{N,RT<:AbstractReachSet{N},FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::Vector{FT}
    ext::Dict{Symbol,Any}
end

function MixedFlowpipe(Fk::Vector{FT},
                       ext::Dict{Symbol,Any}=Dict{Symbol,Any}()) where {N,RT<:AbstractReachSet{N},FT<:Flowpipe{N,RT}}
    return MixedFlowpipe{N,RT,FT}(Fk, ext)
end

# interface functions
array(fp::MixedFlowpipe) = fp.Fk
flowpipe(fp::MixedFlowpipe) = fp
setrep(::Type{MixedFlowpipe{N,RT,FT}}) where {N,RT,FT} = RT
numrsets(fp::MixedFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::MixedFlowpipe, I::Int...) = getindex(fp.Fk, I...)

@inline function tspan(fp::MixedFlowpipe)
    return error("for mixed flowpipes you should specify the index, as in `tspan(fp, i)`")
end
@inline tspan(fp::MixedFlowpipe, i::Int) = tspan(fp[i])

function Base.similar(fp::MixedFlowpipe{N,RT,FT}) where {N,RT,FT}
    return MixedFlowpipe(Vector{FT}())
end

function (fp::MixedFlowpipe)(t::Number)
    return error("not implemented yet")
end

function (fp::MixedFlowpipe)(dt::TimeInterval)
    return error("not implemented yet")
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
    return error("not implemented")
end

function LazySets.linear_map(M, fp::MixedFlowpipe)
    out = [linear_map(M, F) for F in fp]
    return MixedFlowpipe(out, fp.ext)
end

function LazySets.affine_map(M, fp::MixedFlowpipe, b)
    out = [affine_map(M, F, b) for F in fp]
    return MixedFlowpipe(out, fp.ext)
end
