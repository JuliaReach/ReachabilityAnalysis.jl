using .DifferentialEquations
const DE = DifferentialEquations

LazySets.sample(X::IntervalArithmetic.Interval, d::Integer) = sample(convert(Interval, X), d)
LazySets.sample(X::IntervalArithmetic.IntervalBox, d::Integer) = sample(convert(Hyperrectangle, X), d)

# extend the solve API for initial-value problems
function _solve_ensemble(ivp::InitialValueProblem, args...;
                         trajectories=100,
                         trajectories_alg=DE.Tsit5(),
                         ensemble_alg=DE.EnsembleDistributed(),
                         inplace=true,
                         kwargs...)

    # get problem's vector field
    if inplace
        field = inplace_field!(ivp)
    else
        field = outofplace_field(ivp)
    end

    # get time span
    dt = _get_tspan(args...; kwargs...)
    tspan = (inf(dt), sup(dt))

    # sample initial states
    X0 = initial_state(ivp)
    X0_samples = sample(X0, trajectories)

    # formulate ensemble ODE problem
    ensemble_prob = ODEProblem(field, first(X0), tspan)
    _prob_func(prob, i, repeat) = remake(prob, u0 = X0_samples[i])

    ensemble_prob = EnsembleProblem(ensemble_prob, prob_func = _prob_func)
    return DE.solve(ensemble_prob, trajectories_alg, ensemble_alg; trajectories = trajectories)
end
