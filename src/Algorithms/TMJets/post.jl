function post(alg::TMJets{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...) where {N}

    @unpack max_steps, abs_tol, orderT, orderQ = alg

    # initial time and final tim
    t0 = tstart(tspan)
    T = tend(tspan)

    # vector field
    if islinear(ivp) || isaffine(ivp) # TODO: refactor with inplace_field!
        f! = inplace_field!(ivp)
    else
        f! = VectorField(ivp)
    end
    n = statedim(ivp)
    ivp_norm = _normalize(ivp)

    # TEMP should be handled by the normalization
    if ivp_norm isa ConstrainedBlackBoxContinuousSystem
        X = stateset(ivp_norm)
    else
        X = Universe(n)
    end

    # initial set
    X0 = initial_state(ivp_norm)
    box_x0 = box_approximation(X0)
    q0 = center(box_x0)
    δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x", numvars=length(q0), order=2*orderQ)

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, max_steps)

    F, tv, xv, xTM1v = validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT, abs_tol, max_steps, X)

    ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v) # keep Any or add the type param?
    return Flowpipe(F, ext)
end
