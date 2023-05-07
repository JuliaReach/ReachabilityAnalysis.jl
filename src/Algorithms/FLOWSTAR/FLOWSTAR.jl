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

numtype(::FLOWSTAR) = Float64
rsetrep(::FLOWSTAR) = TaylorModelReachSet{Float64,IA.Interval{Float64}}
