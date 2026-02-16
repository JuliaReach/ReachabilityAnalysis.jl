module XFZ18Module

using ..ReachabilityAnalysis: AbstractContinuousPost

using MathematicalSystems: PolynomialContinuousSystem
# using MultivariatePolynomials,
#       DynamicPolynomials,
#       JuMP,
#       PolyJuMP,
#       SumOfSquares,
#       MathOptInterfaceMosek,
#       SemialgebraicSets

import ..ReachabilityAnalysis: post, numtype

include("XFZ18.jl")
include("reach.jl")

export XFZ18

end  # module
