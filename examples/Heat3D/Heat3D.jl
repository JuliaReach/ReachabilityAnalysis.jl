# # Heat PDE
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/Heat3D.ipynb)
#
#md # !!! note "Overview"
#md #     System type: Affine system\
#md #     State dimension: 48\
#md #     Application domain: Mechanical Engineering
#
# ## Model description
#
# The system is described by the linear differential equations:
#
# ```math
#   \begin{array}{lcl}
#   \dot{x}(t) &=& Ax(t) + Bu(t),\qquad u(t) \in \mathcal{U} \\
#   y(t) &=& C x(t)
#   \end{array}
# ```
#
# There are two versions of this benchmark:
#
# - *(time-varying inputs):* The inputs can change arbitrarily over time: $\forall t: u(t)\in \mathcal{U}$.
#
# - *(constant inputs):* The inputs are uncertain only in their initial value, and
#    constant over time: ``u(0)\in \mathcal{U}``, ``\dot u (t)= 0``.

using ReachabilityAnalysis, SparseArrays, JLD2

examples_dir = normpath(@__DIR__, "..", "..", "..", "examples")
HEAT01_path = joinpath(examples_dir, "Heat3D", "HEAT01.jld2")
