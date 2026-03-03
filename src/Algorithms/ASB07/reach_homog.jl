# computes overapproximation of Φ * set(F[k-1]) with a zonotope
# this operations adds n generators, hence we use an order reduction
# function
function reach_homog_ASB07!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                            Ω0::Zonotope{N,VN,MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            ::Universe,
                            recursive::Val{true},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Zₖʳ = Ω0
    @inbounds F[1] = ReachSet(Zₖʳ, Δt)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k < NSTEPS
        Zₖ₋₁ = Zₖʳ
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)

        k += 1
        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
    end
    return F
end

# non-recursive implementation; to get more accurate interval matrix powers Φ^k
# we use the IntervalMatrices.IntervalMatrixPower interface
function reach_homog_ASB07!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                            Ω0::Zonotope{N,VN,MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            ::Universe,
                            recursive::Val{false},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Z₀ = Ω0
    @inbounds F[1] = ReachSet(Z₀, Δt)

    c₀ = Z₀.center
    G₀ = Z₀.generators

    Φpow = IntervalMatrixPower(Φ) # lazy interval matrix power

    k = 1
    @inbounds while k < NSTEPS
        Φ_power_k = matrix(Φpow)
        Φc, Φs = _split(Φ_power_k)

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, c₀, G₀)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)

        Δt += δ
        k += 1
        F[k] = ReachSet(Zₖʳ, Δt)
        increment!(Φpow)
    end
    return F
end

# case with an invariant
function reach_homog_ASB07!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                            Ω0::Zonotope{N,VN,MN},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            max_order::Integer,
                            X::LazySet,
                            recursive::Val{true},
                            reduction_method::AbstractReductionMethod,
                            Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Zₖʳ = Ω0
    @inbounds F[1] = ReachSet(Zₖʳ, Δt)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k < NSTEPS
        Zₖ₋₁ = Zₖʳ
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)
        _isdisjoint(X, Zₖʳ) && break

        k += 1
        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
    end
    if k < NSTEPS
        resize!(F, k)
    end
    return F
end
