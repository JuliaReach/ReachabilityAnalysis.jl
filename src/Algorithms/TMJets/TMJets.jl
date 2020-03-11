@with_kw struct TMJets{N} <: AbstractContinuousPost
    δ::N
    abs_tol::N
    orderQ::Int
    orderT::Int
    #setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
end

# include("post.jl")
# include("reach.jl")
