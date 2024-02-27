# # Heat PDE

#md # !!! note "Overview"
#md #     System type: Affine continuous system\
#md #     State dimension: ``n^3`` for ``n ∈ \{5, 10, 20\}``\
#md #     Application domain: Mechanical engineering

# ## Model description

# The system is described by affine differential equations:
#
# ```math
# \begin{aligned}
#     \dot{x}(t) &= Ax(t) + Bu(t),\qquad u(t) ∈ \mathcal{U} \\
#     y(t) &= C x(t)
# \end{aligned}
# ```
#
# There are two versions of this benchmark:
#
# - **Time-varying inputs**: The inputs can change arbitrarily over time:
#    ``\forall t: u(t) ∈ \mathcal{U}``.
#
# - **Constant inputs**: The inputs are only uncertain in the initial value, and
#    constant over time: ``u(0) ∈ \mathcal{U}``, ``\dot{u}(t) = 0.``

using ReachabilityAnalysis, SparseArrays, JLD2, ReachabilityBase.CurrentPath  #!jl

path = @current_path("Heat3D", "HEAT01.jld2")
