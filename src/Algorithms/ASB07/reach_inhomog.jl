# case with input and without invariant
function reach_inhomog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                              Ω0::Zonotope{N, VN, MN},
                              Φ::AbstractMatrix,
                              NSTEPS::Integer,
                              δ::N,
                              max_order::Integer,
                              X::Universe,
                              U::Zonotope,
                              recursive::Val{true},
                              reduction_method::AbstractReductionMethod,
                              time_shift::N) where {N, VN, MN}
    # initial reach set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # input sequence
    Wk₊ = copy(U)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        Zₖ₋₁ = set(F[k-1])
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖ = minkowski_sum(Wk₊, Zₖ)
        Zₖʳ = _reduce_order(Zₖ, max_order, reduction_method)

        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
        k += 1
        Wk₊ = _overapproximate_interval_linear_map(Φc, Φs, Wk₊.center, Wk₊.generators)
        Wk₊ = _reduce_order(Wk₊, max_order, reduction_method)
    end
    return F
end

# case with input and with invariant
function reach_inhomog_ASB07!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                              Ω0::Zonotope{N, VN, MN},
                              Φ::AbstractMatrix,
                              NSTEPS::Integer,
                              δ::N,
                              max_order::Integer,
                              X::LazySet,
                              U::Zonotope,
                              recursive::Val{true},
                              reduction_method::AbstractReductionMethod,
                              time_shift::N) where {N, VN, MN}
    # initial reach set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # input sequence
    Wk₊ = copy(U)

    # split the interval matrix into center and radius
    Φc, Φs = _split(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        Zₖ₋₁ = set(F[k-1])
        cₖ₋₁ = Zₖ₋₁.center
        Gₖ₋₁ = Zₖ₋₁.generators

        Zₖ = _overapproximate_interval_linear_map(Φc, Φs, cₖ₋₁, Gₖ₋₁)
        Zₖ = minkowski_sum(Wk₊, Zₖ)
        Zₖʳ = _reduce_order(Zₖ, max_order, reduction_method)
        _is_intersection_empty(X, Zₖʳ) && break

        Δt += δ
        F[k] = ReachSet(Zₖʳ, Δt)
        k += 1

        Wk₊ = _overapproximate_interval_linear_map(Φc, Φs, Wk₊.center, Wk₊.generators)
        Wk₊ = _reduce_order(Wk₊, max_order, reduction_method)
    end
    if k < NSTEPS + 1
        resize!(F, k-1)
    end
    return F
end
