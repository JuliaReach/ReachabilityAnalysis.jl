# =======================
# x' = Ax
# =======================

# Set representation: Interval
# Matrix operations: Dense
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0::CartesianProductArray{N}, Φ::AbstractMatrix, NSTEPS, δ,
                               X::Universe, ST::Type{<:Interval{N}}, vars, block_indices,
                               row_blocks::AbstractVector, column_blocks::AbstractVector, Δt0,
                               viewval::Val{true}) where {N}
    return _reach_homog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, ST, vars, block_indices, row_blocks,
                                  column_blocks, Δt0, viewval)
end

# Set representation: Interval
# Matrix operations: Sparse  (TODO currently not exploited; only defined for disambiguation)
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0::CartesianProductArray{N}, Φ::SparseMatrixCSC, NSTEPS, δ,
                               X::Universe, ST::Type{<:Interval{N}}, vars, block_indices,
                               row_blocks::AbstractVector, column_blocks::AbstractVector, Δt0,
                               viewval::Val{true}) where {N}
    return _reach_homog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, ST, vars, block_indices, row_blocks,
                                  column_blocks, Δt0, viewval)
end

function _reach_homog_BFFPSV18!(F, Xhat0::CartesianProductArray{N}, Φ, NSTEPS, δ, ::Universe,
                                ST::Type{<:Interval{N}}, vars, block_indices,
                                row_blocks::AbstractVector, column_blocks::AbstractVector, Δt0,
                                viewval::Val{true}) where {N}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # buffer for overapproximated Minkowski sum of each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))
    zero_int = Interval(zero(N), zero(N))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            Xhatk[i] = zero_int
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                # concrete interval-interval Minkowski sum
                aux = Φpowerk[bi[1], bj[1]] * Xhat0.array[j].dat
                Xhatk[i] = Interval(Xhatk[i].dat + aux)
            end
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        Φpowerk, Φpowerk_cache = Φpowerk_cache, Φpowerk
    end
    return F
end

# Set representation: Generic
# Matrix operations: Dense
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0::CartesianProductArray{N}, Φ::MT, NSTEPS, δ,
                               ::Universe, ST, vars, block_indices,
                               row_blocks::AbstractVector{<:RBLKi},
                               column_blocks::AbstractVector{<:CBLKj}, Δt0,
                               viewval::Val{true}) where {NM,MT<:AbstractMatrix{NM},N,RBLKi,CBLKj}

    # initial reach-set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SubArray{N,2,MT,Tuple{RBLKi,CBLKj},false}
    buffer = Vector{LinearMap{N,ST,N,SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                buffer[j] = view(Φpowerk, bi, bj) * Xhat0.array[j]
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        Φpowerk, Φpowerk_cache = Φpowerk_cache, Φpowerk
    end
    return F
end

# Set representation: Generic
# Matrix operations: Sparse
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0::CartesianProductArray{N}, Φ::MT, NSTEPS, δ,
                               ::Universe, ST, vars, block_indices,
                               row_blocks::AbstractVector{<:RBLKi},
                               column_blocks::AbstractVector{<:CBLKj}, Δt0,
                               viewval::Val{true}) where {NM,MT<:SparseMatrixCSC{NM},N,RBLKi,CBLKj}

    # store first element
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    SMT = SubArray{N,2,MT,Tuple{RBLKi,CBLKj},false}
    buffer = Vector{LinearMap{N,ST,N,SMT}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                Φpowerk_bi_bj = Φpowerk
                if !iszero(Φpowerk_bi_bj)
                    buffer[j] = view(Φpowerk, bi, bj) * Xhat0.array[j]
                end
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        Φpowerk, Φpowerk_cache = Φpowerk_cache, Φpowerk
    end
    return F
end
