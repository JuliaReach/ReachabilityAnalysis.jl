"""
    VREP{N, AM, S, D} <: AbstractContinuousPost

Implementation of a linear reachability solver using the vertex representation.
"""
struct VREP{N,AM,S,D} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    static::S
    dim::D
end

# convenience constructor using symbols
function VREP(; δ::N,
              dim::Union{Int,Missing}=missing,
              static::Bool=!ismissing(dim),
              backend=nothing,
              approx_model=Forward(; sih=:concrete, exp=:base, setops=:vrep, backend=backend)) where {N}
    n = ismissing(dim) ? missing : Val(dim)
    return VREP(δ, approx_model, Val(static), n)
end

step_size(alg::VREP) = alg.δ
numtype(::VREP{N}) where {N} = N

function setrep(::VREP{N,AM,Val{false},Missing}) where {N,AM}
    return VPolytope{N,Vector{N}}
end

function setrep(::VREP{N,AM,Val{false},D}) where {N,AM,D}
    return VPolygon{N,Vector{N}}
end

function setrep(alg::VREP{N,AM,Val{true},Missing}) where {N,AM}
    return error("the set representation of this algorithm requires the dimension to be specified, but it is $(alg.dim)")
end

function setrep(::VREP{N,AM,Val{true},Val{2}}) where {N,AM}
    return VPolygon{N,SVector{2,N}}
end

function setrep(::VREP{N,AM,Val{true},Val{n}}) where {N,AM,n}
    return VPolytope{N,SVector{n,N}}
end

function rsetrep(alg::VREP{N}) where {N}
    ST = setrep(alg)
    return ReachSet{N,ST}
end
