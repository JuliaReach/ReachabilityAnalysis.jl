# =========================================
# Conservative time discretization methods
# =========================================

using Reexport

include("DiscretizationModule.jl")

include("ApplySetops.jl")

include("Overapproximate.jl")
@reexport using ..Overapproximate

# Approximation model in discrete time, i.e. without bloating
include("NoBloatingModule.jl")
@reexport using ..NoBloatingModule

# Forward approximation
include("ForwardModule.jl")
@reexport using ..ForwardModule

# Backward approximation
include("BackwardModule.jl")
@reexport using ..BackwardModule

# Intersect one step forward in time with one step backward
include("StepIntersectModule.jl")
@reexport using ..StepIntersectModule

# Discretize using the correction hull of the matrix exponential
include("CorrectionHullModule.jl")
@reexport using ..CorrectionHullModule

# First-order approximation from d/dt
include("SecondOrderddtModule.jl")
@reexport using ..SecondOrderddtModule

# First-order approximation with zonotope
include("FirstOrderZonotopeModule.jl")
@reexport using ..FirstOrderZonotopeModule

# First-order approximation
include("FirstOrderModule.jl")
@reexport using ..FirstOrderModule

# Forward-Backward discretization using continuous convex hull
include("ForwardBackwardModule.jl")
@reexport using ..ForwardBackwardModule
