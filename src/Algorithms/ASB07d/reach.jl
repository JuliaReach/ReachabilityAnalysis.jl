function reach_ASB07_decomposed!(R::Vector{<:ReachSet},
                                 Ωhat0::Vector{<:LazySet},
                                 U::Union{ConstantInput, Nothing},
                                 ϕ::AbstractMatrix{IntervalArithmetic.Interval{NUM}},
                                 N::Int,
                                 δ::Float64,
                                 max_order::Int,
                                 n::Int,
                                 partition,
                                 blocks) where {NUM}
    # initial reach set
    t0, t1 = zero(δ), δ
    R[1] = ReachSet(CartesianProductArray(Ωhat0[blocks]), t0, t1)

    b = length(blocks)
    Rₖ_array = Vector{LazySet{NUM}}(undef, b)
    ϕpowerk = copy(ϕ)

    if U != nothing
        Whatk = Vector{Zonotope{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = linear_map(proj(bi, n), inputs)
        end
    end

    k = 2
    while k <= N
        for i in 1:b
            bi = partition[blocks[i]]
            Rₖ_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                block = IntervalMatrix(ϕpowerk[bi, bj])
                if !iszero(block)
                    Rₖ_bij = overapproximate(block * Ωhat0[j], Zonotope)
                    Rₖ_bi = minkowski_sum(Rₖ_bi, Rₖ_bij)
                end
            end
            if U != nothing
                Rₖ_bi = minkowski_sum(Rₖ_bi, Whatk[i])
            end
            Rₖ_bi = reduce_order(Rₖ_bi, max_order)  # reduce order
            Rₖ_array[i] = Rₖ_bi
        end
        Rₖ = CartesianProductArray(copy(Rₖ_array))

        # store reach set
        t0 = t1
        t1 += δ
        R[k] = ReachSet(Rₖ, t0, t1)

        ϕpowerk *= ϕ

        k += 1
    end
    return R
end
