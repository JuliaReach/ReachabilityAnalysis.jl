# =================================================================
# Bouncing-ball model
#
# We model the bouncing ball as a hybrid automaton with one location and a self-loop.
# See for example [[LG09]](@ref) pp. 79-83.
#

function bouncing_ball(; X0 = Hyperrectangle(low=[10.0, 0.0], high=[10.2, 0.0]))

    # "falling" mode with invariant x >= 0
    invariant = HalfSpace([-1.0, 0.0], 0.0)
    flow = @system(z' = [0.0 1.0; 0.0 0.0] * z + [0.0, -9.81], z ∈ invariant)

    # guard x ≤ 0 && v ≤ 0
    guard = HPolyhedron([HalfSpace([1.0, 0.0], 0.0), HalfSpace([0.0, 1.0], 0.0)])

    # reset map v⁺ := -cv
    assignment = ConstrainedLinearMap([1.0 0.0; 0.0 -0.75], guard)

    # initial-value problem
    H = HybridSystem(flow, assignment)
    prob = @ivp(H, z(0) ∈ X0)
    tspan = (0.0, 5.0)

    return prob, tspan
end
