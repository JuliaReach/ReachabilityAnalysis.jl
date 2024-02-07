# =========================================
# Conservative time discretization methods
# =========================================

using Reexport
include("DiscretizationModule.jl")

include("ApplySetops.jl")

include("Overapproximate.jl")
@reexport using ..Overapproximate

# Approximation model in discrete time, i.e. without bloating
using ..DiscretizationModule # TODO Remove
using ..Exponentiation # TODO Remove
include("NoBloatingModule.jl")
@reexport using ..NoBloatingModule

# Forward approximation
include("ForwardModule.jl")
@reexport using ..ForwardModule

# Backward approximation
include("BackwardModule.jl")
@reexport using ..BackwardModule

# Intersect one step forward in time with one step backward
include("StepIntersect.jl")

# Discretize using the correction hull of the matrix exponential
include("CorrectionHullModule.jl")
@reexport using ..CorrectionHullModule

# First-order approximation from d/dt
include("SecondOrderddt.jl")

# First-order approximation with zonotope
include("FirstOrderZonotope.jl")

# First-order approximation
include("FirstOrderModule.jl")
@reexport using ..FirstOrderModule

# Forward-Backward discretization using continuous convex hull
include("ForwardBackwardModule.jl")
@reexport using ..ForwardBackwardModule
