@with_kw struct BFFPSV18 <: AbstractContinuousPost
    δ::Float64
    approx_model::AbstractApproximationModel=ForwardApproximation(sih_method="concrete",
                                                                  exp_method="base",
                                                                  set_operations="zonotope",
                                                                  phi2_method="base")
    max_order::Int=10
    #setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
end
