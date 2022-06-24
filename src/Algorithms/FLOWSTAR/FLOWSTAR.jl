Base.@kwdef struct FLOWSTAR <: AbstractContinuousPost
    remainder_estimation::Float64=1e-5
    precondition::String="QR"
    order_min::Int64=4
    order_max::Int64=10
    cutoff::Float64=1e-20
    precision::Int64=53
end

numtype(::FLOWSTAR) = Float64
rsetrep(::FLOWSTAR) = TaylorModelReachSet{Float64, IA.Interval{Float64}}
