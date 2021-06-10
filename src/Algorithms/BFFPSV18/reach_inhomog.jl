# =======================
# x' = Ax + U
# =======================

# Set representation: Interval
# Matrix operations: Dense
# Invariant: No
function reach_inhomog_BFFPSV18!(F, Xhat0::LazySet{N}, Φ::MT, NSTEPS, δ, X::Universe, U,
                                 ST::Type{<:Interval{N}}, vars,
                                 block_indices,
                                 row_blocks::AbstractVector{<:RBLKi},
                                 column_blocks::AbstractVector{<:CBLKj},
                                 Δt0,
                                 viewval::Val{true}) where {NM, MT<:AbstractMatrix{NM}, N, RBLKi, CBLKj}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # buffer for overapproximated Minkowski sum of each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))
    zero_int = Interval(zero(N), zero(N))

    # preallocate accumulated inputs and decompose it
    Whatk = Vector{ST}(undef, length(row_blocks))
    @inbounds for (i, bi) in enumerate(row_blocks)
        πX = Projection(U, bi)
        Whatk[i] = overapproximate(πX, ST)
    end

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            Xhatk[i] = zero_int
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                # concrete interval-interval Minkowski sum
                aux = Φpowerk[bi[1], bj[1]] * Xhat0.array[j].dat
                Xhatk[i] = Interval(Xhatk[i].dat + aux)
            end
            Xhatk[i] = Interval(Xhatk[i].dat + Whatk[i].dat)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        # update the input set
        for (i, bi) in enumerate(row_blocks)
            Whatk[i] = overapproximate(Whatk[i] ⊕ view(Φpowerk, bi, :) * U, ST)
        end

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end

# Set representation: Generic
# Matrix operations: Dense
# Invariant: No
function reach_inhomog_BFFPSV18!(F, Xhat0::LazySet{N}, Φ::MT, NSTEPS, δ, X::Universe, U,
                                 ST, vars, block_indices,
                                 row_blocks::AbstractVector{<:RBLKi},
                                 column_blocks::AbstractVector{<:CBLKj},
                                 Δt0,
                                 viewval::Val{true}) where {NM, MT<:AbstractMatrix{NM}, N, RBLKi, CBLKj}
``
    # initial reach-set
    Δt = (zero(N) .. δ) + Δt0

    # overapproximate with template
    #Xhat0 = CartesianProductArray([overapproximate(xi, PolarDirections(20)) for xi in Xhat0.array])

    aux0 = Xhat0.array[block_indices]
    aux0 = convert(Vector{LazySet{N}}, aux0)
    @inbounds F[1] = SparseReachSet(CartesianProductArray(aux0), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SubArray{N, 2, MT, Tuple{RBLKi, CBLKj}, false}
    ST = typeof(Xhat0.array[1])
    buffer = Vector{LazySet{N}}(undef, length(column_blocks)) # Vector{LazySets.LinearMap{N, ST, N, SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{LazySet{N}}(undef, length(row_blocks))

    # preallocate accumulated inputs and decompose it
    Whatk = Vector{LazySet{N}}(undef, length(row_blocks))
    @inbounds for (i, bi) in enumerate(row_blocks)
        πX = Projection(U, bi)
        Whatk[i] = πX
    end

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                buffer[j] = view(Φpowerk, bi, bj) * Xhat0.array[j]
            end
            Xhatk[i] = MinkowskiSumArray(copy(buffer)) ⊕ copy(Whatk[i])
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        # update the input set
        for (i, bi) in enumerate(row_blocks)
            Whatk[i] = Whatk[i] ⊕ view(Φpowerk, bi, :) * U
        end

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end

# Set representation: Generic
# Matrix operations: Sparse
# Invariant: No
function reach_inhomog_BFFPSV18!(F, Xhat0, Φ::MT, NSTEPS, δ::N,
                                 X::Universe, U,
                                 ST, vars,
                                 block_indices,
                                 row_blocks::AbstractVector{<:RBLKi},
                                 column_blocks::AbstractVector{<:CBLKj},
                                 Δt0,
                                 viewval::Val{true}) where {NM, IM, MT<:SparseMatrixCSC{NM, IM},
                                                            N, RBLKi, CBLKj}
    # store first element
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SubArray{N, 2, MT, Tuple{RBLKi, CBLKj}, false}
    buffer = Vector{LazySets.LinearMap{N, ST, N, SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    # preallocate accumulated inputs and decompose it
    Whatk = Vector{ST}(undef, length(row_blocks))
    @inbounds for (i, bi) in enumerate(row_blocks)
        πX = Projection(U, bi)
        Whatk[i] = overapproximate(πX, ST)
    end

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                Φpowerk_bi_bj = Φpowerk
                if !iszero(Φpowerk_bi_bj)
                    buffer[j] = view(Φpowerk, bi, bj) * Xhat0.array[j]
                end
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer) ⊕ Whatk[i], ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        # update the input set
        for (i, bi) in enumerate(row_blocks)
            Whatk[i] = overapproximate(Whatk[i] ⊕ view(Φpowerk, bi, :) * U, ST)
        end

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end

function reach_inhomog_BFFPSV18!(F, Xhat0, Φ::MT, NSTEPS, δ::N,
                                 X::Universe, U,
                                 ST, vars,
                                 block_indices,
                                 row_blocks::AbstractVector{<:RBLKi},
                                 column_blocks::AbstractVector{<:CBLKj},
                                 Δt0,
                                 viewval::Val{false}) where {NM, IM, MT<:SparseMatrixCSC{NM, IM}, N, RBLKi, CBLKj}

    # store first element
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SparseMatrixCSC{NM, IM}
    buffer = Vector{LazySets.LinearMap{N, ST, NM, SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    # preallocate accumulated inputs and decompose it
    Whatk = Vector{ST}(undef, length(row_blocks))
    @inbounds for (i, bi) in enumerate(row_blocks)
        πX = Projection(U, bi)
        Whatk[i] = overapproximate(πX, ST)
    end

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                Φpowerk_bi_bj = Φpowerk
                if !iszero(Φpowerk_bi_bj)
                    buffer[j] = Φpowerk[bi, bj] * Xhat0.array[j]
                end
            end

            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer) ⊕ Whatk[i], ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        # update the input set
        for (i, bi) in enumerate(row_blocks)
            Whatk[i] = overapproximate(Whatk[i] ⊕ Φpowerk[bi, :] * U, ST)
        end

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end





#######################################################################

# Set representation: Generic
# Matrix operations: Dense
# Invariant: LazySet
function reach_inhomog_BFFPSV18!(F, Xhat0::LazySet{N}, Φ::MT, NSTEPS, δ, X::LazySet, U,
                                 ST, vars, block_indices,
                                 row_blocks::AbstractVector{<:RBLKi},
                                 column_blocks::AbstractVector{<:CBLKj},
                                 Δt0,
                                 viewval::Val{true}) where {NM, MT<:AbstractMatrix{NM}, N, RBLKi, CBLKj}

    # initial reach-set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SubArray{N, 2, MT, Tuple{RBLKi, CBLKj}, false}
    buffer = Vector{LazySets.LinearMap{N, ST, N, SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    # preallocate accumulated inputs and decompose it
    Whatk = Vector{ST}(undef, length(row_blocks))
    @inbounds for (i, bi) in enumerate(row_blocks)
        πX = Projection(U, bi)
        Whatk[i] = overapproximate(πX, ST)
    end

    k = 2
    @inbounds while k <= NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                buffer[j] = view(Φpowerk, bi, bj) * Xhat0.array[j]
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer) ⊕ Whatk[i], ST)
        end
        Xk = CartesianProductArray(copy(Xhatk))
        _is_intersection_empty(X, Xk) && break
        Δt += δ
        F[k] = SparseReachSet(Xk, Δt, vars)

        # update the input set
        for (i, bi) in enumerate(row_blocks)
            Whatk[i] = overapproximate(Whatk[i] ⊕ view(Φpowerk, bi, :) * U, ST)
        end

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)

        k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k-1)
    end
    return F
end
