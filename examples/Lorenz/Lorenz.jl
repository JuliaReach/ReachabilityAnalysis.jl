# # Lorenz equations
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/Lorenz.ipynb)
#
#md # !!! note "Overview"
#md #     System type: Continuous blackbox system\
#md #     State dimension: 3\
#md #     Application domain: Atmospheric convection
#
# ## Model description
#
# The model is a system of three ordinary differential equations now known as
# the Lorenz equations:
#
# ```math
#   \begin{array}{lcl}
#   \dfrac{dx}{dt} & = & \sigma (y-x), \\ \\
#   \dfrac{dy}{dt} & = & x(\rho - z) - y, \\ \\
#   \dfrac{dz}{dt} & = & xy-\beta z
#   \end{array}
# ```
#
# The equations relate the properties of a two-dimensional fluid layer uniformly
# warmed from below and cooled from above. In particular, the equations describe
# the rate of change of three quantities with respect to time:``x`` is
# proportional to the rate of convection, ``y`` to the horizontal temperature
# variation, and ``z`` to the vertical temperature variation. The constants
# ``\sigma`` , ``\rho`` , and ``\beta``  are system parameters proportional to
# the Prandtl number, Rayleigh number, and certain physical dimensions of the
# layer itself.

using ReachabilityAnalysis, Plots

@taylorize function lorenz!(dx, x, params, t)
    local σ = 10.0
    local β = 8.0 / 3.0
    local ρ = 28.0
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
    return dx
end

# ## Reachability settings
#
# The initial values considered are ``X_0 \in [0.9, 1.1] \times [0, 0] \times [0, 0]``,
# for a time span of `10`.

X0 = Hyperrectangle(low=[0.9, 0.0, 0.0], high=[1.1, 0.0, 0.0])
prob = @ivp(x' = lorenz!(x), dim=3, x(0) ∈ X0);

# ## Results

# We compute the flowpipe using the TMJets algorithm with ``n_T=10`` and ``n_Q=2``.

alg = TMJets(abs_tol=1e-15, orderT=10, orderQ=2, max_steps=50_000);

sol = solve(prob, T=10.0, alg=alg)

# Below we plot the flowpipe projected on the `(1, 3)` plane.

solz = overapproximate(sol, Zonotope);

plot(solz, vars=(1, 3))
