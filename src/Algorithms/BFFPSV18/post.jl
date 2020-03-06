# TODO: dispatch on linear system only.. ?
function post(alg::BFFPSV18, ivp::IVP{<:AbstractContinuousSystem}, tspan, args...; kwargs...)

    # get time horizon
    T = _get_T(tspan)

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # discretize system
    δ = step_size(alg)
    ivp_discr = discretize(ivp_norm, δ, approx_model(alg))
    Ω0 = initial_state(ivp_discr)
    # TODO: convert or overapproximate to alg.ST
    #Ω0 = _convert_or_overapproximate(Zonotope, Ω0) # returns a Zonotope
    #Ω0 = Zonotope(center(Ω0), Matrix(genmat(Ω0))) # TEMP hard-code generators matrix
    #ZT = Zonotope{Float64, Vector{Float64}, Matrix{Float64}} # TODO: typeof(Ω0) # should be a concretely typed Zonotope

    # decompose into a cartesian product
    Xhat0 = _decompose(Ω0, partition(alg), setrep(alg))
    Φ = state_matrix(ivp_discr)
    N = eltype(Φ) # or num_type(alg)
    X = stateset(ivp_discr) # invariant

    # preallocate output flowpipe
    NSTEPS = round(Int, T / δ)
    SETREP = setrep(alg)
    ST = CartesianProductArray{N, SETREP}
    F = Vector{SparseReachSet{N, ST}}(undef, NSTEPS)

    if hasinput(ivp)
        error("not implemented yet")
    else
        reach_homog!(F, Xhat0, Φ, NSTEPS, δ, X, vars, block_indices,
                     row_blocks, column_blocks, NUM, ST)
    end

    return Flowpipe(F)
end
