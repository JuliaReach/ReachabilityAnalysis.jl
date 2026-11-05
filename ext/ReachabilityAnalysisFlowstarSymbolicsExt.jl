module ReachabilityAnalysisFlowstarSymbolicsExt

import ReachabilityAnalysis.FLOWSTARModule
import Flowstar, Symbolics
using Flowstar: FlowstarSetting, ContinuousReachModel
import IntervalArithmetic as IA
using IntervalBoxes: IntervalBox
using LazySets: Hyperrectangle, LazySet, overapproximate
using MathematicalSystems: initial_state, statedim, system

function FLOWSTARModule._extract_model(ivp, step_size::Number, T, order, remainder_estimation,
                                       precondition, cutoff, precision, verbose, scheme)
    s, n = system(ivp), statedim(ivp)

    # vector field in string form
    fstr = convert(Vector{Tuple{String,String}}, s)

    # initial set
    X0 = initial_state(ivp)
    if X0 isa AbstractVector{<:IA.Interval}
        dom = X0
    elseif X0 isa IntervalBox
        dom = [X0[i] for i in 1:length(X0)]
    elseif X0 isa LazySet
        IB = convert(IntervalBox, overapproximate(X0, Hyperrectangle))
        dom = [IB[i] for i in 1:length(IB)]
    else
        throw(ArgumentError("invalid `X0` of type $(typeof(X0))"))
    end

    name = "flowstar" # TODO parse name of f

    # Plot first two variables if possible, otherwise plot the first variable vs time
    # TODO disable this option
    plot_states = n > 1 ? "$(fstr[1][1]),$(fstr[2][1])" : "t,$(fstr[1][1])" # seems mandatory?
    gen_plots = false # ignored?

    # names of state variables
    states = join([fi[1] for fi in fstr], ",")

    # parameters; currently not supported
    params = nothing

    # equations of motion
    eom = join(["$(fi[1])' = $(fi[2])" for fi in fstr], "\n ")

    setting = FlowstarSetting(step_size, T, order, name, remainder_estimation,
                              precondition, cutoff, precision, plot_states,
                              verbose; gen_plots)

    return ContinuousReachModel(states, params, setting, scheme, eom, dom)
end

end  # module
