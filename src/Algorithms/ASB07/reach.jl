function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::Float64,
                            max_order::Integer,
                            X::Universe) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    Ω0red = reduce_order(Ω0, max_order)
    F[1] = ReachSet(Ω0red, Δt)

    k = 2
    while k <= NSTEPS
        Rₖ = overapproximate(Φ * set(F[k-1]), Zonotope)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end
