# computes overapproximation of Φ * set(F[k-1]) with a zonotope
# this operations adds n generators, hence we use an order reduction
# function
function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::Float64,
                            max_order::Integer,
                            X::Universe,
                            recursive::Val{true}) where {N, VN, MN}
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
        Zₖ₊₁ʳ = _reduce_order(Zₖ₊₁, max_order)
        Δt += δ
        F[k] = ReachSet(Zₖ₊₁ʳ, Δt)
    end
    return F
end


function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::Float64,
                            max_order::Integer,
                            X::Universe,
                            recursive::Val{false}) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    Z0 = Ω0
    c0 = Z0.center
    G0 = Z0.generators

    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        # split the interval matrix into center and radius
        Φc, Φs = _split(Φ_power_k)

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, c0, G0)
        Zₖʳ = _reduce_order(Zₖ, max_order)
        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end
