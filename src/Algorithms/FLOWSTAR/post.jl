function post(alg::FLOWSTAR{ST,OT,PT,IT}, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI, kwargs...) where {ST,OT,PT,IT}
    @required Flowstar

    @unpack step_size, order, remainder_estimation, precondition, cutoff, precision, verbose, scheme = alg

    # initial time and final time
    t0, T = tstart(timespan), tend(timespan)
    @assert t0 == 0.0

    # extract model file
    model = if isa(ivp.s, BlackBoxContinuousSystem{String})
        # the model is provided in .flow form
        # consistency of t0, T, X0 with the model is not checked
        ivp.s.f
    else
        s, n = system(ivp), statedim(ivp)

        # vector field in string form
        @required Symbolics
        fstr = convert(Vector{Tuple{String,String}}, s)

        # initial set
        X0 = initial_state(ivp)
        dom = isa(X0, IntervalBox) ? X0 : convert(IntervalBox, overapproximate(X0, Hyperrectangle))

        name = "flowstar" # TODO parse name of f

        # Plot first two variables if possible, otherwise plot the first variable vs time
        # TODO disable this option
        plot_states = n > 1 ? "$(fstr[1][1]),$(fstr[2][1])" : "t,$(fstr[1][1])" # seems mandatory?
        gen_plots = false # ignored?

        setting = FlowstarSetting(step_size, T, order, name, remainder_estimation, precondition,
                                  cutoff, precision, plot_states, verbose; gen_plots)

        # names of state variables
        states = join([fi[1] for fi in fstr], ",")

        # parameters; currently not supported
        params = nothing

        # equations of motion
        eom = join(["$(fi[1])' = $(fi[2])" for fi in fstr], "\n ")

        ContinuousReachModel(states, params, setting, scheme, eom, dom)
    end

    # call Flow*; use Val(true) to obtain TaylorModel1 of TaylorN
    sol = FlowstarContinuousSolution(model, Val(true))
    flow = sol.flow

    # preallocate output flowpipe
    F = Vector{TaylorModelReachSet{Float64,IA.Interval{Float64}}}()
    sizehint!(F, length(flow))

    counter = t0
    for Fi in flow
        dt = domain(first(Fi))
        δt = TimeInterval(counter + dt)
        counter += tend(dt)
        Ri = TaylorModelReachSet(Fi, δt + Δt0)
        push!(F, Ri)
    end
    ext = Dict{Symbol,Any}()
    return Flowpipe(F, ext)
end
