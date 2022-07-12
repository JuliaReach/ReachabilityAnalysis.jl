using .Symbolics

# required to avoid conflicts with IntervalMatrices
using IntervalMatrices: infimum, supremum

include("../Continuous/symbolics.jl")

function Base.convert(::Type{Vector{Tuple{String, String}}}, s::BlackBoxContinuousSystem; params=[], t=nothing)
    @requires Symbolics

    n = statedim(s)
    f! = MathematicalSystems.mapping(s)

    # get right-hand side symbolic expression then turn it into a vector of strings
    Symbolics.@variables x[1:n]
    dx = Vector{Symbolics.Num}(undef, n)
    f!(dx, x, params, t)
    rhs = string.(dx)

    # make implicit times operators explicit, eg. 2x into 2*x
    rhs = rhs .|> Meta.parse .|> string
    
    # remove square brackets, eg. turning x[1] into x1 
    pat = ["x[$i]" => "x$i" for i in 1:n]
    rhs = [Symbolics.replace(fi, pat...) for fi in rhs]

    return [("x$i", rhs[i]) for i in 1:n]
end
