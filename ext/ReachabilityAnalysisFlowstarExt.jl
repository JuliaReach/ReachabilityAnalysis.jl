module ReachabilityAnalysisFlowstarExt

import ReachabilityAnalysis.FLOWSTARModule
using ReachabilityAnalysis.FLOWSTARModule: FLOWSTAR
import Flowstar
using Flowstar: AbstractODEScheme, AbstractPreconditioner, AbstractTMOrder,
                 ContinuousReachModel, FixedTMOrder, FlowstarContinuousSolution,
                 FlowstarSetting, NonPolyODEScheme, QRPreconditioner
import IntervalArithmetic as IA

function FLOWSTARModule.FLOWSTAR(; δ::Union{Float64,Nothing}=nothing,
                                 step_size::Union{Float64,IA.Interval{Float64},Nothing}=nothing,
                                 order::AbstractTMOrder=FixedTMOrder(6),
                                 remainder_estimation=1e-4,
                                 precondition::AbstractPreconditioner=QRPreconditioner(),
                                 cutoff::Float64=1e-20,
                                 precision=53,
                                 verbose=false,
                                 scheme=NonPolyODEScheme())
    step_size = !isnothing(δ) ? δ :
                (!isnothing(step_size) ? step_size :
                 throw(ArgumentError("the step size should be specified")))
    T = FLOWSTAR{typeof(step_size),typeof(order),typeof(precondition),typeof(scheme)}
    return T(step_size, order, remainder_estimation, precondition, cutoff,
             precision, verbose, scheme)
end

function FLOWSTARModule._model(step_size::Number, T, order, name,
                               remainder_estimation, precondition, cutoff,
                               precision, plot_states, verbose, gen_plots,
                               states, params, scheme, eom, dom)
    setting = FlowstarSetting(step_size, T, order, name, remainder_estimation,
                              precondition, cutoff, precision, plot_states,
                              verbose; gen_plots)
    return ContinuousReachModel(states, params, setting, scheme, eom, dom)
end

# second argument is only used for making this method more specific
function FLOWSTARModule._solution(model, dummy::Nothing)
    return FlowstarContinuousSolution(model, Val(true))
end

end  # module
