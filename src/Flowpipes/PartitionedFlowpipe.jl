# ============================================
# Hybrid flowpipe of possibly different types
# ============================================

"""
    PartitionedFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of possibly different types.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes
"""
struct PartitionedFlowpipe{T,S<:Tuple} <: AbstractFlowpipe # TODO: ask <:AbstractFlowpipe for each element in the tuple..?
    Fk::ArrayPartition{T,S}
    ext::Dict{Symbol,Any}
end
