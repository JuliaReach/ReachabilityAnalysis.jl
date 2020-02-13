ckearmodule ReachabilityAnalysis

# ==================
# Load dependencies
# ==================
include("init.jl")

# ===========================================================
# Structures to represent solutions of reachability problems
# ===========================================================
include("reachsets.jl")
include("flowpipes.jl")
include("operators.jl")
include("solutions.jl")

# ================================================
# Pre-processing functions for continuous systems
# ================================================

include("normalization.jl")
include("discretization.jl")

# ===============================
# Reachability solver algorithms
# ===============================

# Continuous post-operators for linear systems
#include("ContinuousPost/GLGM06/LGG09.jl")
include("ContinuousPost/GLGM06/GLGM06.jl")
include("ContinuousPost/BFFPSV18/BFFPSV18.jl")
#include("ContinuousPost/ASB07/ASB07.jl")
#include("ContinuousPost/ASB07d/ASB07d.jl")
#include("ContinuousPost/A17/A17.jl")

# Continuous post-operators for non-linear systems
#include("ContinuousPost/TMJets/TMJets.jl")
#include("ContinuousPost/A13/A13.jl")
#include("ContinuousPost/KA19/KA19.jl")

# ===========================================
# Discrete post-operators for hybrid systems
# ===========================================

#include("DiscretePost/concrete.jl")
#include("DiscretePost/decomposed.jl")
#include("DiscretePost/lazy.jl")

# =========
# User API
# =========

#include("logging.jl")
include("solve.jl")

end # module
