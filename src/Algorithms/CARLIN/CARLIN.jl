using CarlemanLinearization

"""
    CARLIN <: AbstractContinuousPost

Implementation of the reachability method using Carleman linearization from [1].

[1] Forets, Marcelo, and Christian Schilling. "Reachability of weakly nonlinear systems using Carleman linearization."
    International Conference on Reachability Problems. Springer, Cham, 2021.
"""
Base.@kwdef struct CARLIN <: AbstractContinuousPost
    N::Int=2 # order of the algorithm
    compressed::Bool=false # choose to use compressed Kronecker form
    δ::Float64=0.1 # step size for the linear reachability solver
    bloat::Bool=true # choose to include the error estimate in the result
    resets::Union{Int,Vector{Float64}}=0 # choose the number of resets (equal spacing) or a vector specifying the reset times within the time interval [0, T]
end

include("kronecker.jl")
include("post.jl")
include("reach.jl")
