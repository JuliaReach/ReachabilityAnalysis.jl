# default continuous post

# discrete post
function post(alg::BFFPSV18{N,ST}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N,ST}
    @unpack δ, approx_model, vars, block_indices,
    row_blocks, column_blocks = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)
    Ω0 = initial_state(ivp)
    X = stateset(ivp)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

    # decompose the initial states into a cartesian product
    # TODO add option to do the lazy decomposition
    Xhat0 = decompose(Ω0, column_blocks, ST)

    # force using sparse type for the matrix exponential
    if alg.sparse
        Φ = SparseArrays.sparse(Φ)
    end

    # variables that are actually computed
    vars = reduce(vcat, alg.row_blocks)

    # preallocate output flowpipe
    CP = CartesianProductArray{N,ST}
    F = Vector{SparseReachSet{N,CP,length(vars)}}(undef, NSTEPS)

    # option to use array views
    viewval = Val(alg.view)

    if got_homogeneous
        reach_homog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, ST,
                              vars, block_indices,
                              row_blocks, column_blocks, Δt0, viewval)

    else
        U = inputset(ivp)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        reach_inhomog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, U, ST,
                                vars, block_indices,
                                row_blocks, column_blocks, Δt0, viewval)
    end
    return Flowpipe(F)
end
