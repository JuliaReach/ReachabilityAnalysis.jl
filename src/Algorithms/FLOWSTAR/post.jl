function post(alg::FLOWSTAR{ST,OT,PT,IT}, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroT, kwargs...) where {ST,OT,PT,IT}
    @unpack step_size, order, remainder_estimation, precondition, cutoff, precision, verbose, scheme = alg

    # initial time and final time
    t0, T = tstart(timespan), tend(timespan)
    @assert t0 == 0.0

    # extract model file
    if isa(ivp.s, BlackBoxContinuousSystem{String})
        # the model is provided in .flow form
        # consistency of t0, T, X0 with the model is not checked
        model = ivp.s.f
    else
        # extract model from Julia code (via Symbolics)
        model = _extract_model(ivp, step_size, T, order, remainder_estimation,
                               precondition, cutoff, precision, verbose, scheme)
    end

    # call Flow*; use Val(true) to obtain TaylorModel1 of TaylorN
    sol = _solution(model, nothing)
    flow = sol.flow

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{Float64,IA.Interval{Float64}}}()
    sizehint!(F, length(flow))

    counter = t0
    for Fi in flow
        dt = domain(first(Fi))
        δt = TimeInterval(counter + dt)
        counter += sup(dt)
        Ri = TaylorModelReachSet(Fi, δt + Δt0)
        push!(F, Ri)
    end
    ext = Dict{Symbol,Any}()
    return Flowpipe(F, ext)
end

# see ReachabilityAnalysisFlowstarSymbolicsExt.jl
function _extract_model(ivp, step_size, T, order, remainder_estimation,
                        precondition, cutoff, precision, verbose, scheme)
    mod = Base.get_extension(@__MODULE__, :ReachabilityAnalysisFlowstarSymbolicsExt)
    require(mod, [:Flowstar, :Symbolics]; fun_name="post")
    return error()
end

# see ReachabilityAnalysisFlowstarExt.jl
function _solution(model, dummy)
    mod = Base.get_extension(@__MODULE__, :ReachabilityAnalysisFlowstarExt)
    require(mod, :Flowstar; fun_name="post")
    return error()
end
