function post(alg::TMJets{N}, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI,
              external::Bool=false,    # if `true`, use the external solver defined in TaylorModels.jl
              kwargs...) where {N}

    @unpack max_steps, abs_tol, orderT, orderQ, disjointness, adaptive, min_abs_tol = alg

    # initial time and final time
    t0 = tstart(timespan)
    T = tend(timespan)

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

    if external
        return _solve_external(f!, X0, t0, T, orderQ, orderT, abs_tol, max_steps, Δt0; kwargs...)
    end

    X0tm = _initialize(X0, orderQ, orderT)

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, max_steps)

    F, tv, xv, xTM1v, success, _t0 = validated_integ!(F, f!, X0tm, t0, T, orderQ, orderT,
                                                     abs_tol, max_steps, X, disjointness, Δt0, adaptive)

    if success || !adaptive
        ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v, :actual_abs_tol=>[abs_tol])
        return Flowpipe(F, ext)
    end

    # save extra data, one vector per iteration
    tv_vec = Vector{typeof(tv)}()
    xv_vec = Vector{typeof(xv)}()
    xTM1v_vec = Vector{typeof(xTM1v)}()

    push!(tv_vec, tv)
    push!(xv_vec, xv)
    push!(xTM1v_vec, xTM1v)

    abs_tol_vec = [abs_tol]

    while !success
            # adapt the absolute tolerance
            if abs_tol > min_abs_tol
                abs_tol = abs_tol / 10
                push!(abs_tol_vec, abs_tol)
            else
                @warn("Minimum absolute tolerance, $min_abs_tol reached.")
                ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v) # keep Any or add the type param?
                return Flowpipe(F, ext)
            end

            # new initial states
            if !isempty(F)
                # here we pass the box overapproximation of the final reach-set
                X0_end_box = overapproximate(F[end], Hyperrectangle)
                X0tm = _initialize(X0_end_box, orderQ, orderT)
            end

            # new flowpipe
            Fk = Vector{TaylorModelReachSet{N}}()
            sizehint!(Fk, max_steps)

            Fk, tv, xv, xTM1v, success, _t0 = validated_integ!(Fk, f!, X0tm, _t0, T, orderQ, orderT,
                                                               abs_tol, max_steps, X, disjointness, Δt0, adaptive)

            # append the new flowpipe to the accumulated flowpipe and extra data
            append!(F, Fk)

            push!(tv_vec, tv)
            push!(xv_vec, xv)
            push!(xTM1v_vec, xTM1v)
    end

    ext = Dict{Symbol, Any}(:tv => tv_vec, :xv => xv_vec, :xTM1v => xTM1v_vec, :actual_abs_tol=>abs_tol_vec)
    return Flowpipe(F, ext)
end

function _solve_external(f!, X0, t0, T, orderQ, orderT, abs_tol, max_steps, Δt0; kwargs...)
    N = eltype(X0)

    # box overapproximation of the initial states
    X0box = convert(IntervalBox, box_approximation(X0))

    # extract solver name and options
    solver_name = haskey(kwargs, :solver_name) ? kwargs[:solver_name] : TM.validated_integ
    solver_kwargs = haskey(kwargs, :solver_kwargs) ? kwargs[:solver_kwargs] : Dict(:maxsteps=>max_steps)

    # call external solver
    tv, xv, xTM1v = solver_name(f!, X0box, t0, T, orderQ, orderT, abs_tol; solver_kwargs...)

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

_initialize(X0::TaylorModelReachSet, orderQ, orderT) = set(X0)
_initialize(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT) where {T} = X0

_initialize(X0::AbstractHyperrectangle, orderQ, orderT) = convert(IntervalBox, box_approximation(X0))

# fallback
_initialize(X0::LazySet, orderQ, orderT) = _initialize(box_approximation(X0), orderQ, orderT)

function _initialize(X0::AbstractZonotope, orderQ, orderT)
    X = overapproximate(X0, TaylorModelReachSet, orderQ=orderQ, orderT=orderT)
    return set(X)
end
