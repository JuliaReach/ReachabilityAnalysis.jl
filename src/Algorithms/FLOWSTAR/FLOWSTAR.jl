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

function load_Flowstar_constructor()
    return quote
        function FLOWSTAR(; δ::Union{Float64,Nothing}=nothing,
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
            return T(step_size, order, remainder_estimation, precondition, cutoff, precision,
                     verbose, scheme)
        end
    end
end

numtype(::FLOWSTAR) = Float64
rsetrep(::FLOWSTAR) = TaylorModelReachSet{Float64,IA.Interval{Float64}}
