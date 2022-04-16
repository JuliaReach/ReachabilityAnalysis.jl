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
    @ivp(x'= A*x + b, x(0) ∈ X0, x ∈ X)
    tspan = (0.0, 20.0)
    return prob, tspan
end

@taylorize function _projectile!(dx, x, p, t)
    local m = 1.0
    local g = 9.81

    # variables: x, y, xdot, ydot
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = zero(x[3])
    dx[4] = -g * one(x[4])
    return dx
end

function projectile_taylor(; X=Universe(4))
    X0 = Singleton([0., 1000.0, 10., 40.0])
    prob = @ivp(x'= _projectile!(x), dim:4, x(0) ∈ X0, x ∈ X)
    tspan = (0.0, 20.0)
    return prob, tspan
end
