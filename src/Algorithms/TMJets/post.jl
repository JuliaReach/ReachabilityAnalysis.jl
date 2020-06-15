function post(alg::TMJets{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              time_shift::N=zero(N), kwargs...) where {N}

    @unpack max_steps, abs_tol, orderT, orderQ, disjointness, adaptive, min_abs_tol = alg

    # initial time and final time
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
    X = stateset(ivp_norm)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x", numvars=n, order=2*orderQ)

    # initial set
    X0 = initial_state(ivp_norm)
    box_x0 = box_approximation(X0)
    q0 = center(box_x0)
    δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, max_steps)

    F, tv, xv, xTM1v, success, _t0 = _validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT,
                                      abs_tol, max_steps, X, disjointness, time_shift, adaptive)

    if success || !adaptive
        ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v) # keep Any or add the type param?
        return Flowpipe(F, ext)
    end

    # save extra data, one vector per iteration
    #tv_vec = Vector{typeof(tv)}()
    #xv_vec = Vector{typeof(xv)}()
    #xTM1v_vec = Vector{typeof(xTM1v)}()

    #push!(tv_vec, tv)
    #push!(xv_vec, xv)
    #push!(xTM1v_vec, xTM1v)

    while !success
            # adapt the absolut tolerance
            if abs_tol > min_abs_tol
                abs_tol = abs_tol / 10
            else
                @warn("Minimum absolute tolerance, $min_abs_tol reached.")
                ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v) # keep Any or add the type param?
                return Flowpipe(F, ext)
            end

            # new initial states
            X0 = overapproximate(F[end], Zonotope) |> set
            box_x0 = box_approximation(X0)
            q0 = center(box_x0)
            δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

            # new flowpipe
            Fk = Vector{TaylorModelReachSet{N}}()
            sizehint!(Fk, max_steps)
            Fk, tv, xv, xTM1v, success, _t0 = _validated_integ!(Fk, f!, q0, δq0, _t0, T, orderQ, orderT,
                                                                abs_tol, max_steps, X, disjointness, time_shift, adaptive)

            # append the new flowpipe to the accumulated flowpipe and extra data
            append!(F, Fk)
            #push!(tv_vec, copy(tv))
            #push!(xv_vec, copy(xv))
            #push!(xTM1v_vec, copy(xTM1v))
    end
    #ext = Dict{Symbol, Any}(:tv => tv_vec, :xv => xv_vec, :xTM1v => xTM1v_vec) # keep Any or add the type param?
    #return Flowpipe(F, ext)
    return Flowpipe(F)
end
