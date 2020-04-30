# =====================================
# Projectile
#
# system type: affine continuous system
# state dimension: 4
# =====================================

function projectile(; X=Universe(4))
    m = 1.0
    g = 9.81

    # variables: x, y, xdot, ydot
    A = [0 0 1 0;
         0 0 0 1;
         0 0 0 0;
         0 0 0 0.]
    X0 = Singleton([0., 1000.0, 10., 40.0])
    b = [0., 0., 0., -g]
    prob = @ivp(x'= A*x + b, x(0) ∈ X0, x ∈ X)
    tspan = (0.0, 20.0)
    return prob, tspan
end
