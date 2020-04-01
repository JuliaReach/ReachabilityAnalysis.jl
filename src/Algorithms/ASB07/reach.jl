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
    println("max order = $max_order")
    k = 1
    @inbounds while k <= NSTEPS - 1
        Zk = set(F[k])
        ck = Zk.center
        Gk = Zk.generators
        Rₖ = _overapproximate_interval_linear_map(Φc, Φs, ck, Gk)
        println(typeof(Rₖ))
        println(order(Rₖ))
        Rₖ = _reduce_order(Rₖ, max_order)
        println(typeof(Rₖ))
        println(order(Rₖ))
        Δt += δ
        k += 1
        F[k] = ReachSet(Rₖ, Δt)
    end
    return F
end
