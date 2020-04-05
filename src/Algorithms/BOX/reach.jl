# ================
# Homogeneous case
# ================

function reach_homog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                          Ω0::Hyperrectangle{N, VNC, VNR},
                          Φ::AbstractMatrix,
                          NSTEPS::Integer,
                          δ::Float64,
                          X::Universe) where {N, VNC, VNR}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = overapproximate(Φ * set(F[k-1]), Hyperrectangle)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end
