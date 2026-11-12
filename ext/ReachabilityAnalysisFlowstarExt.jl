module ReachabilityAnalysisFlowstarExt

import ReachabilityAnalysis.FLOWSTARModule
using ReachabilityAnalysis.FLOWSTARModule: FLOWSTAR
import Flowstar
using Flowstar: AbstractPreconditioner, AbstractTMOrder, FixedTMOrder,
                FlowstarContinuousSolution, NonPolyODEScheme, QRPreconditioner
import IntervalArithmetic as IA

function FLOWSTARModule._FLOWSTAR(; δ::Union{Float64,Nothing}=nothing,
                                  step_size::Union{Float64,IA.Interval{Float64},Nothing}=nothing,
                                  remainder_estimation=1e-4,
                                  cutoff::Float64=1e-20,
                                  precision=53,
                                  verbose=false,
                                  order::AbstractTMOrder=FixedTMOrder(6),
                                  precondition::AbstractPreconditioner=QRPreconditioner(),
                                  scheme=NonPolyODEScheme())
    step_size = !isnothing(δ) ? δ :
                (!isnothing(step_size) ? step_size :
                 throw(ArgumentError("the step size should be specified")))
    T = FLOWSTAR{typeof(step_size),typeof(order),typeof(precondition),typeof(scheme)}
    return T(step_size, order, remainder_estimation, precondition, cutoff,
             precision, verbose, scheme)
end

# second argument is only used for making this method more specific
function FLOWSTARModule._solution(model, dummy::Nothing)
    return FlowstarContinuousSolution(model, Val(true))
end

end  # module
