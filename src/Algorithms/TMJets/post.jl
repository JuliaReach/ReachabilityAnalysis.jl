function post(alg::TMJets{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI,
              external::Bool=false,    # if `true`, use the external solver defined in TaylorModels.jl
              kwargs...) where {N}

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

    # FIXME refactor
    if external

        # pass box overapproximation of the initial states
        box_x0 = box_approximation(X0)
        q0 = center(box_x0)
        δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

        solver_name = haskey(kwargs, :solver_name) ? kwargs[:solver_name] : TM.validated_integ
        solver_kwargs = haskey(kwargs, :solver_kwargs) ? kwargs[:solver_kwargs] : Dict(:maxsteps=>max_steps)
        tv, xv, xTM1v = solver_name(f!, q0, δq0, t0, T, orderQ, orderT,
                                    abs_tol; solver_kwargs...)
        # build flowpipe
        F = Vector{TaylorModelReachSet{N}}()
        sizehint!(F, max_steps)
        for i in 2:length(tv)
            δt = TimeInterval(tv[i-1], tv[i])
            Ri = TaylorModelReachSet(xTM1v[:, i], δt + Δt0)
            push!(F, Ri)
        end
        ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v)
        return Flowpipe(F, ext)
    end

    # TODO generalize to using _convert_or_overapproximate(Zonotope, X0)
    #Ω0 = isa(X0, Zonotope) ? X0 : box_approximation(X0)

    # TEMP to use with NNA#mforets/TORA
    #X0tm = overapproximate(Ω0, TaylorModelReachSet)
    @assert isa(X0, Zonotope)
    if order(X0) == 1
        X0tm = overapproximate(X0, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)
    elseif order(X0) == 2
        X0tm = _overapproximate_structured(X0, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)
    else
        error("expected order 1 or 2")
    end

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, max_steps)

    F, tv, xv, xTM1v, success, _t0 = validated_integ!(F, f!, X0tm, t0, T, orderQ, orderT,
                                                     abs_tol, max_steps, X, disjointness, Δt0, adaptive)

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
            # adapt the absolute tolerance
            if abs_tol > min_abs_tol
                abs_tol = abs_tol / 10
            else
                @warn("Minimum absolute tolerance, $min_abs_tol reached.")
                ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v) # keep Any or add the type param?
                return Flowpipe(F, ext)
            end

            # new initial states
            if !isempty(F)
                # here we pass the last TM directly, without approximation
                X0s = F[end]
                X0tm = overapproximate(X0s, TaylorModelReachSet)
            end

            # new flowpipe
            Fk = Vector{TaylorModelReachSet{N}}()
            sizehint!(Fk, max_steps)

            Fk, tv, xv, xTM1v, success, _t0 = validated_integ!(Fk, f!, X0tm, _t0, T, orderQ, orderT,
                                                               abs_tol, max_steps, X, disjointness, Δt0, adaptive)

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
