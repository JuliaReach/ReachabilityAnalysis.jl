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
include("NoBloating.jl")

# Forward approximation
include("Forward.jl")

# Backward approximation
include("BackwardModule.jl")
@reexport using ..BackwardModule

# Intersect one step forward in time with one step backward
include("StepIntersect.jl")

# Discretize using the correction hull of the matrix exponential
include("CorrectionHull.jl")

# First-order approximation from d/dt
include("SecondOrderddt.jl")

# First-order approximation with zonotope
include("FirstOrderZonotope.jl")

# First-order approximation
include("FirstOrder.jl")

# Forward-Backward discretization using continuous convex hull
include("ForwardBackward.jl")
