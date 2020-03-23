function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::Float64,
                            max_order::Integer,
                            X::Universe) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)
    #println(diam_norm(Φ))
    k = 2
    @inbounds while k <= NSTEPS
        #println("k = $k")

        #println(norm(set(F[k-1])))

        # TODO: fix and use faster version _overapproximate
        Rₖ = overapproximate(Φ * set(F[k-1]), Zonotope)

        #println(norm(Rₖ))

        Rₖ = reduce_order(Rₖ, max_order)

        #println(norm(Rₖ))

        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end
