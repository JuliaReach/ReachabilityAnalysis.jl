@with_kw struct TMJets <: AbstractContinuousPost
    δ::Float64
    approx_model::AbstractApproximationModel=ForwardApproximation(sih_method="concrete",
                                                                  exp_method="base",
                                                                  set_operations="zonotope",
                                                                  phi2_method="base")
    max_order::Int=10
    #setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
end

step_size(alg::GLGM06) = alg.δ
approx_model(alg::GLGM06) = alg.approx_model
max_order(alg::GLGM06) = alg.max_order

include("post.jl")
include("reach.jl")
