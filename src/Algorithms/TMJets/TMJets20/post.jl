function post(alg::TMJets20{N}, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI,
              external::Bool=false,    # if `true`, use the external solver defined in TaylorModels.jl
              kwargs...) where {N}

    @unpack orderQ, orderT, abstol, maxsteps, adaptive, minabstol, disjointness = alg

    # initial time and final time
    t0 = tstart(timespan)
    T = tend(timespan)

    # vector field
    if islinear(ivp) || isaffine(ivp)
        f! = inplace_field!(ivp)
    else
        f! = VectorField(ivp)
    end
    parse_eqs = get(kwargs, :parse_eqs, false)

    n = statedim(ivp)
    ivp_norm = _normalize(ivp)
    X = stateset(ivp_norm)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x", numvars=n, order=2*orderQ)

    # initial set
    X0 = initial_state(ivp_norm)

    if external
        return _solve_external(f!, X0, t0, T, orderQ, orderT, abstol, maxsteps, Δt0; parse_eqs=parse_eqs, kwargs...)
    end

    X0tm = _initialize(X0, orderQ, orderT)

    # optionally absorb initial remainder
    shrink_wrapping = get(kwargs, :shrink_wrapping, true)
    if shrink_wrapping && isa(X0tm, TaylorModelReachSet)
        if !all(iszero, remainder(X0tm))
            X0tm = _shrink_wrapping(X0tm)
        end
    end

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, maxsteps)

    F, tv, xv, xTM1v, success, _t0 = validated_integ!(F, f!, X0tm, t0, T, orderQ, orderT,
                                                      abstol, maxsteps, X, disjointness, Δt0, adaptive; parse_eqs=parse_eqs)

    if success || !adaptive
        ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v, :actual_abs_tol=>[abstol])
        return Flowpipe(F, ext)
    end

    # save extra data, one vector per iteration
    tv_vec = Vector{typeof(tv)}()
    xv_vec = Vector{typeof(xv)}()
    xTM1v_vec = Vector{typeof(xTM1v)}()

    push!(tv_vec, tv)
    push!(xv_vec, xv)
    push!(xTM1v_vec, xTM1v)

    abstol_vec = [abstol]

    while !success
            # adapt the absolute tolerance
            if abstol > minabstol
                abstol = abstol / 10
                push!(abstol_vec, abstol)
            else
                @warn("Minimum absolute tolerance, $minabstol reached.")
                ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v)
                return Flowpipe(F, ext)
            end

            # new initial states
            if !isempty(F)
                # here we pass the box overapproximation of the final reach-set
                # X0_end_box = set(overapproximate(F[end], Hyperrectangle, orderQ=orderQ, orderT=orderT))
                # X0tm = _initialize(X0_end_box, orderQ, orderT)
                X0tm = xv_vec[end][end]
            end

            # new flowpipe
            Fk = Vector{TaylorModelReachSet{N}}()
            sizehint!(Fk, maxsteps)

            Fk, tv, xv, xTM1v, success, _t0 = validated_integ!(Fk, f!, X0tm, _t0, T, orderQ, orderT,
                                                               abstol, maxsteps, X, disjointness, Δt0, adaptive; parse_eqs=parse_eqs)

            # append the new flowpipe to the accumulated flowpipe and extra data
            append!(F, Fk)

            push!(tv_vec, tv)
            push!(xv_vec, xv)
            push!(xTM1v_vec, xTM1v)
    end

    ext = Dict{Symbol, Any}(:tv => tv_vec, :xv => xv_vec, :xTM1v => xTM1v_vec, :actual_abs_tol=>abstol_vec)
    return Flowpipe(F, ext)
end
