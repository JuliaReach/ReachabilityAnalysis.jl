struct FLOWSTAR{ST,OT,PT,IT} <: AbstractContinuousPost
    step_size::ST
    order::OT
    remainder_estimation::Float64
    precondition::PT
    cutoff::Float64
    precision::Int64
    verbose::Bool
    scheme::IT
end

function FLOWSTAR(; kwargs...)
    mod = Base.get_extension(@__MODULE__, :ReachabilityAnalysisFlowstarExt)
    require(mod, :Flowstar; fun_name="FLOWSTAR")
    return _FLOWSTAR(; kwargs...)
end

# two-step function because of kwargs specificity
# see ReachabilityAnalysisFlowstarExt.jl
function _FLOWSTAR(dummy)
    mod = Base.get_extension(@__MODULE__, :ReachabilityAnalysisFlowstarExt)
    require(mod, :Flowstar; fun_name="FLOWSTAR")
    return error()
end

numtype(::FLOWSTAR) = Float64
rsetrep(::FLOWSTAR) = TaylorModelReachSet{Float64,IA.Interval{Float64}}
