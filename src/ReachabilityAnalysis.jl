module ReachabilityAnalysis

# ===========================================================
# Dependencies & User API
# ===========================================================
include("Initialization/init.jl")
include("Initialization/exports.jl")
include("Continuous/solve.jl")

# ===========================================================
# Discretization
# ===========================================================
include("Discretization/Exponentiation.jl")
include("Discretization/discretization.jl")

# ===========================================================
# Structures to represent solutions of reachability problems
# ===========================================================
include("Flowpipes/setops.jl")
include("ReachSets/reachsets.jl")
include("Flowpipes/flowpipes.jl")
include("Flowpipes/clustering.jl")
include("Flowpipes/solutions.jl")

# ===========================================================
# Shared functionality for continuous systems
# ===========================================================
include("Continuous/fields.jl")
include("Continuous/normalization.jl")
include("Continuous/homogenization.jl")
include("Continuous/linearization.jl")
include("Continuous/setops.jl")

# ===========================================================
# Reachability solver algorithms
# ===========================================================

# Linear
include("Algorithms/A20/A20Module.jl")
using ..A20Module: A20

include("Algorithms/ASB07/ASB07Module.jl")
using ..ASB07Module: ASB07

include("Algorithms/BFFPSV18/BFFPSV18Module.jl")
using ..BFFPSV18Module: BFFPSV18

include("Algorithms/BOX/BOXModule.jl")
using ..BOXModule: BOX

include("Algorithms/GLGM06/GLGM06Module.jl")
using ..GLGM06Module: GLGM06

include("Algorithms/INT/INTModule.jl")
using ..INTModule: INT

include("Algorithms/LGG09/LGG09Module.jl")
using ..LGG09Module: LGG09

include("Algorithms/ORBIT/ORBITModule.jl")
using ..ORBITModule: ORBIT

include("Algorithms/VREP/VREPModule.jl")
using ..VREPModule: VREP

include("Algorithms/linear_post.jl")

# Linear parametric
include("Algorithms/HLBS25/HLBS25.jl")

# Nonlinear
include("Algorithms/CARLIN/CARLINModule.jl")
using ..CARLINModule: CARLIN, CanonicalQuadraticForm

include("Algorithms/FLOWSTAR/FLOWSTARModule.jl")
using ..FLOWSTARModule: FLOWSTAR

include("Algorithms/TMJets/TMJets21a/TMJets21a.jl")

include("Algorithms/TMJets/TMJets21b/TMJets21b.jl")

include("Algorithms/TMJets/common.jl")

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
