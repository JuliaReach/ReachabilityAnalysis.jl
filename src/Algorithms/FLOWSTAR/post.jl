using .Flowstar: FlowstarContinuousSolution

function post(alg::FLOWSTAR, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI, kwargs...)
    println("kwargs = $kwargs")
    @requires Flowstar

    @unpack remainder_estimation, precondition, order_min, order_max, cutoff, precision = alg

    # initial time and final time
    t0, T = tstart(timespan), tend(timespan)

    # vector field
    f! = (islinear(ivp) || isaffine(ivp)) ? inplace_field!(ivp) : VectorField(ivp)

    # initial set
    ivp_norm = _normalize(ivp)
    X0 = initial_state(ivp_norm)

    # extract model file if present
    model = if !haskey(kwargs, :model)
        # temporary
        throw(ArgumentError("the model file needs to be passed as a keyword argument"))
    else
        kwargs[:model]
    end

    # call Flow*
    sol = FlowstarContinuousSolution(model)
    flow = sol.flow

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{Float64, IA.Interval{Float64}}}()
    sizehint!(F, length(flow))

    counter = t0
    for Fi in flow
        dt = domain(first(Fi))
        δt = TimeInterval(counter + dt)
        counter += sup(dt)
        Ri = TaylorModelReachSet(Fi, δt + Δt0)
        push!(F, Ri)
    end
    ext = Dict{Symbol, Any}()
    return Flowpipe(F, ext)
end
