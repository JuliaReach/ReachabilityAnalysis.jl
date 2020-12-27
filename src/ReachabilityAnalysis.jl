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
include("Flowpipes/setops.jl")
include("Flowpipes/reachsets.jl")
include("Flowpipes/flowpipes.jl")
include("Flowpipes/arrayops.jl")
include("Flowpipes/clustering.jl")
include("Flowpipes/fields.jl")
include("Flowpipes/solutions.jl")

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
include("Algorithms/A20/A20.jl")
include("Algorithms/LGG09/LGG09.jl")
include("Algorithms/GLGM06/GLGM06.jl")
include("Algorithms/INT/INT.jl")
include("Algorithms/BOX/BOX.jl")
include("Algorithms/BFFPSV18/BFFPSV18.jl")
include("Algorithms/ASB07/ASB07.jl")
include("Algorithms/ORBIT/ORBIT.jl")

# Nonlinear
include("Algorithms/TMJets/TMJets.jl")
include("Algorithms/QINT/QINT.jl")

# ===========================================================
# Discrete post-operators for hybrid systems
# ===========================================================
include("Hybrid/constructors.jl")
include("Hybrid/waiting_list.jl")
include("Hybrid/transitions.jl")
include("Hybrid/time_triggered.jl")
include("Hybrid/solve.jl")

# visualization
include("Flowpipes/recipes.jl")

end # module
