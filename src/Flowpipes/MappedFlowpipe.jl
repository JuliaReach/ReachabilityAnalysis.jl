# =====================================
# Flowpipe composition with a lazy map
# =====================================

"""
    MappedFlowpipe{FT<:AbstractFlowpipe, ST} <: AbstractFlowpipe

### Fields

- `F`    -- flowpipe
- `func` -- function representing the map
"""
struct MappedFlowpipe{FT<:AbstractFlowpipe,ST} <: AbstractFlowpipe
    F::FT
    func::ST
end

"""
    Projection(fp::AbstractFlowpipe, vars::NTuple{D,Int}) where {D}

Return the lazy projection of a flowpipe.

### Input

### Output

### Notes

The projection is lazy, and consists of mapping each set
`X` in the flowpipe to `MX`, where `M` is the projection matrix associated with
the given variables `vars`.
"""
function LazySets.Projection(fp::AbstractFlowpipe, vars::NTuple{D,Int}) where {D}
    # TODO: assert that vars belong to the variables of the flowpipe
    M = projection_matrix(collect(vars), dim(F), Float64)
    func = @map(x -> M * x)
    return MappedFlowpipe(fp, func)
end

function overapproximate(fp::Flowpipe, args...)
    return Flowpipe(map(R -> overapproximate(R, args...), fp), fp.ext)
end

function overapproximate(fp::AbstractVector{<:AbstractReachSet}, args...)
    return Flowpipe(map(R -> overapproximate(R, args...), fp))
end
