"""
Module for defining and checking properties.
"""
module Properties

using LazySets

import LazySets.dim
import LazySets.Approximations: project

# ==============================
# Property struct and evaluation
# ==============================
include("Property.jl")
export Property,
       dim,
       check,
       project

include("Conjunction.jl")
include("Disjunction.jl")
export Conjunction,
       Disjunction

include("BadStatesProperty.jl")
export BadStatesProperty

include("SafeStatesProperty.jl")
export SafeStatesProperty

end # module
