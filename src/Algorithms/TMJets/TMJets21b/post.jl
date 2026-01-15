function post(alg::TMJets21b{N}, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI,
              kwargs...) where {N}
    @unpack orderQ, orderT, abstol, maxsteps, absorb, adaptive, minabstol,
    validatesteps, ε, δ, absorb_steps, disjointness = alg

    # initial time and final time
    t0 = tstart(timespan)
    T = tend(timespan)

    # vector field
    if islinear(ivp) || isaffine(ivp) # TODO: refactor with inplace_field!
        f! = inplace_field!(ivp)
    else
        f! = VectorField(ivp)
    end

    # algorithm optional arguments
    parse_eqs = get(kwargs, :parse_eqs, false)
    params = get(kwargs, :params, nothing)

    n = statedim(ivp)
    ivp_norm = _normalize(ivp)
    X = stateset(ivp_norm)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x"; numvars=n, order=2 * orderQ)

    # initial set
    X0 = initial_state(ivp_norm)
    X0tm = _initialize(X0, orderQ, orderT)

    # optionally absorb initial remainder
    shrink_wrapping = get(kwargs, :shrink_wrapping, true)
    if shrink_wrapping && isa(X0tm, TaylorModelReachSet)
        if !all(iszero, remainder(X0tm))
            X0tm = _shrink_wrapping(X0tm)
        end
    end

    # call external solver
    TMSol = TM.validated_integ2(f!, X0tm, t0, T, orderQ, orderT,
                                abstol, params;
                                parse_eqs=parse_eqs,
                                maxsteps=maxsteps,
                                absorb=absorb,
                                adaptive=adaptive,
                                minabstol=minabstol,
                                validatesteps=validatesteps,
                                ε=ε, δ=ε, absorb_steps=absorb_steps)
    tv = TMSol.time
    xv = TMSol.fp
    xTM1v = TMSol.xTM

    # build flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, maxsteps)

    # loop over reach-sets (the first reach-set at the initial time point is ignored)
    @inbounds for i in 2:length(tv)

        # create Taylor model reach-set
        δt = TimeInterval(tv[i] + domain(xTM1v[1, i]))
        Ri = TaylorModelReachSet(xTM1v[:, i], δt + Δt0)

        # check intersection with invariant
        _isdisjoint(Ri, X, disjointness) && break

        push!(F, Ri)
    end
    ext = Dict{Symbol,Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v)
    return Flowpipe(F, ext)
end
