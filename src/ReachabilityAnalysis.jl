module ReachabilityAnalysis

# ===========================================================
# Dependencies & User API
# ===========================================================
include("Initialization/init.jl")
include("Initialization/exports.jl")
include("Continuous/solve.jl")

# ===========================================================
# Structures to represent solutions of reachability problems
# ===========================================================
include("Flowpipes/reachsets.jl")
include("Flowpipes/setops.jl")
include("Flowpipes/flowpipes.jl")
include("Flowpipes/operators.jl")
include("Flowpipes/fields.jl")
include("Flowpipes/solutions.jl")
include("Flowpipes/recipes.jl")

# ===========================================================
# Shared functionality for continuous systems
# ===========================================================
include("Continuous/normalization.jl")
include("Continuous/exponentiation.jl")
include("Continuous/discretization.jl")

# ===========================================================
# Reachability solver algorithms
# ===========================================================
# Linear
#include("Algorithms/GLGM06/LGG09.jl")
include("Algorithms/GLGM06/GLGM06.jl")
include("Algorithms/INT/INT.jl")
include("Algorithms/BOX/BOX.jl")
include("Algorithms/BFFPSV18/BFFPSV18.jl")
include("Algorithms/ASB07/ASB07.jl")
#include("Algorithms/ASB07d/ASB07d.jl")
#include("Algorithms/A17/A17.jl")

# Nonlinear
include("Algorithms/TMJets/TMJets.jl")
#include("Algorithms/A13/A13.jl")
#include("Algorithms/KA19/KA19.jl")

# ===========================================================
# Discrete post-operators for hybrid systems
# ===========================================================
include("Hybrid/time_triggered.jl")
#include("DiscretePost/concrete.jl")
#include("DiscretePost/decomposed.jl")
#include("DiscretePost/lazy.jl")

end # module
