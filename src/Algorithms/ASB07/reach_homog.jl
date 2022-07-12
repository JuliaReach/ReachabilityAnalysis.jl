# computes overapproximation of Φ * set(F[k-1]) with a zonotope
# this operations adds n generators, hence we use an order reduction
# function
function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            X::Universe,
                            recursive::Val{true},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N, TN, VN, MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k <= NSTEPS - 1
        Zk = set(F[k])
        ck = Zk.center
        Gk = Zk.generators

        Zₖ₊₁ = _overapproximate_interval_linear_map(Φc, Φs, ck, Gk)
        Zₖ₊₁ʳ = reduce_order(Zₖ₊₁, max_order, reduction_method)

        k += 1
        Δt += δ
        F[k] = ReachSet(Zₖ₊₁ʳ, Δt)
    end
    return F
end

# non-recursive implementation; to get more accurate interval matrix powers Φ^k
# we use the IntervalMatrices.IntervalMatrixPower interface
function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            X::Universe,
                            recursive::Val{false},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N, TN, VN, MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)
    Z0 = Ω0
    c0 = Z0.center
    G0 = Z0.generators

    Φpow = IntervalMatrixPower(Φ) # lazy interval matrix power

    k = 2
    @inbounds while k <= NSTEPS
        Φ_power_k = get(Φpow)
        Φc, Φs = _split(Φ_power_k)

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, c0, G0)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)

        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
        increment!(Φpow)
        k += 1
    end
    return F
end

# case with an invariant
function reach_homog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                            Ω0::Zonotope{N, VN, MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            X::LazySet,
                            recursive::Val{true},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N, TN, VN, MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        Zₖ₋₁ = set(F[k-1])
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)
        _is_intersection_empty(X, Zₖʳ) && break

        k += 1
        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
    end
    if k < NSTEPS + 1
        resize!(F, k-1)
    end
    return F
end
