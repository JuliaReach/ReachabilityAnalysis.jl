# computes overapproximation of Φ * set(F[k-1]) with a zonotope
# this operations adds n generators, hence we use an order reduction
# function
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

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k <= NSTEPS - 1
        Zk = set(F[k])
        ck = Zk.center
        Gk = Zk.generators

        k = k + 1
        Zₖ₊₁ = _overapproximate_interval_linear_map(Φc, Φs, ck, Gk)
        Rₖ₊₁ = _reduce_order(Zₖ₊₁, max_order)
        Δt += δ
        F[k] = ReachSet(Rₖ₊₁, Δt)
    end
    return F
end
