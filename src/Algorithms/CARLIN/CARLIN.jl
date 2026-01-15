using CarlemanLinearization: build_matrix, error_bound_specabs

"""
    CARLIN <: AbstractContinuousPost

Implementation of the reachability method using Carleman linearization from [ForetsS21](@citet).
"""
Base.@kwdef struct CARLIN <: AbstractContinuousPost
    N::Int = 2 # order of the algorithm
    compress::Bool = false # choose to use compressed Kronecker form
    δ::Float64 = 0.1 # step size for the linear reachability solver
    bloat::Bool = true # choose to include the error estimate in the result
    resets::Union{Int,Vector{Float64}} = 0 # choose the number of resets (equal spacing) or a vector specifying the reset times within the time interval [0, T]
end

include("kronecker.jl")
include("post.jl")
include("reach.jl")
