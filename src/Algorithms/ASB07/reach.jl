function reach_ASB07!(R::Vector{<:ReachSet},
                      Ω0::LazySet,
                      U::Union{ConstantInput, Nothing},
                      Φ::AbstractMatrix,
                      N::Int,
                      δ::Float64,
                      max_order::Int)
    # initial reach set
    t0, t1 = zero(δ), δ
    R[1] = ReachSet(Ω0, t0, t1)

    k = 2
    while k <= N
        Rₖ = overapproximate(Φ * R[k-1].X, Zonotope)
        if U != nothing
            Rₖ = minkowski_sum(Rₖ, next_set(U))
        end
        Rₖ = reduce_order(Rₖ, max_order)  # reduce order

        # store reach set
        t0 = t1
        t1 += δ
        R[k] = ReachSet(Rₖ, t0, t1)

        k += 1
    end
    return R
end
