# case with input and without invariant
function reach_inhomog_ASB07!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                              Ω0::Zonotope{N,VN,MN},
                              Φ::AbstractMatrix,
                              NSTEPS::Integer,
                              δ::N,
                              max_order::Integer,
                              ::Universe,
                              U::Zonotope,
                              recursive::Val{true},
                              reduction_method::AbstractReductionMethod,
                              Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Zₖʳ = Ω0
    @inbounds F[1] = ReachSet(Zₖʳ, Δt)

    # input sequence
    Wk₊ = copy(U)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k < NSTEPS
        Zₖ₋₁ = Zₖʳ
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖ = minkowski_sum(Wk₊, Zₖ)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)

        Δt += δ
        k += 1
        F[k] = ReachSet(Zₖʳ, Δt)
        Wk₊ = _overapproximate_interval_linear_map(Φc, Φs, Wk₊.center, Wk₊.generators)
        Wk₊ = reduce_order(Wk₊, max_order, reduction_method)
    end
    return F
end

# case with input and with invariant
function reach_inhomog_ASB07!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                              Ω0::Zonotope{N,VN,MN},
                              Φ::AbstractMatrix,
                              NSTEPS::Integer,
                              δ::N,
                              max_order::Integer,
                              X::LazySet,
                              U::Zonotope,
                              recursive::Val{true},
                              reduction_method::AbstractReductionMethod,
                              Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Zₖʳ = Ω0
    @inbounds F[1] = ReachSet(Zₖʳ, Δt)

    # input sequence
    Wk₊ = copy(U)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 1
    @inbounds while k < NSTEPS
        Zₖ₋₁ = Zₖʳ
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖ = minkowski_sum(Wk₊, Zₖ)
        Zₖʳ = reduce_order(Zₖ, max_order, reduction_method)
        _isdisjoint(X, Zₖʳ) && break

        Δt += δ
        k += 1
        F[k] = ReachSet(Zₖʳ, Δt)

        Wk₊ = _overapproximate_interval_linear_map(Φc, Φs, Wk₊.center, Wk₊.generators)
        Wk₊ = reduce_order(Wk₊, max_order, reduction_method)
    end
    if k < NSTEPS
        resize!(F, k)
    end
    return F
end
