# ==============================================================================
# Five-dimensional system
#
# Example from [Althoff10; Section 3.2.3](@citet)
#
#
# system type: continuous LTI system
# state dimension: 5
# input dimension: 1
#
# This five-dimensional system is taken from [LeGuernic09; Example 4.1](@citet).
# ==============================================================================

function linear5D_homog()
    A = Float64[-1 -4 0 0 0;
                4 -1 0 0 0;
                0 0 -3 1 0;
                0 0 -1 -3 0;
                0 0 0 0 -2]
    X0 = Hyperrectangle(; low=fill(0.9, 5), high=fill(1.1, 5))
    prob = @ivp(x' = Ax, x(0) ∈ X0)
    tspan = (0.0, 5.0)
    return prob, tspan
end

function linear5D()
    # system matrix
    D = [-1.0 -4.0 0.0 0.0 0.0;
         4.0 -1.0 0.0 0.0 0.0;
         0.0 0.0 -3.0 1.0 0.0;
         0.0 0.0 -1.0 -3.0 0.0;
         0.0 0.0 0.0 0.0 -2.0]
    P = [0.6 -0.1 0.1 0.7 -0.2;
         -0.5 0.7 -0.1 -0.8 0.0;
         0.9 -0.5 0.3 -0.6 0.1;
         0.5 -0.7 0.5 0.6 0.3;
         0.8 0.7 0.6 -0.3 0.2]
    A = P * D * inv(P)
    X = Universe(5) # state domain
    U = Ball2(zeros(5), 0.01) # input domain

    X0 = BallInf([1.0, 0.0, 0.0, 0.0, 0.0], 0.1)
    prob = @ivp(x' = Ax + u, x ∈ X, u ∈ U, x(0) ∈ X0) # continuous LTI system
    tspan = (0.0, 5.0)
    return prob, tspan
end
