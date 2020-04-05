# ==============================================================================
# Motor
#
# Part of the SLICOT benchmarks
#
# system type: continuous LTI system
# state dimension: 8
# input dimension: 2
# ==============================================================================
using SparseArrays

function state_matrix_motor()
    # system matrix
    I = [1, 2, 2, 3, 3, 3, 3, 4, 5, 6, 6, 7, 7, 7, 7, 8]
    J = [2, 3, 2, 1, 2, 3, 4, 1, 6, 7, 6, 5, 6, 7, 8, 5]
    V = [1, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0, 1.0,
         8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0]
    A = sparse(I, J, V)
end

function initial_state_motor()
    # initial set:
    # xᵢ ∈ [0.002, 0.0025] for i = 1
    # xᵢ ∈ [0.001, 0.0015] for i = 5
    # xᵢ = 0 otherwise
    X0 = Hyperrectangle(low=[0.002, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0],
                        high=[0.0025, 0.0, 0.0, 0.0, 0.0015, 0.0, 0.0, 0.0])
end

function motor_homg()
    A = state_matrix_motor()
    X0 = initial_state_motor()

    prob = @ivp(x' = Ax + Bu, x(0) ∈ X0) # continuous LTI system
    tspan = (0.0, 20.0)
    return prob, tspan
end

function motor()
    A = state_matrix_motor()
    B = sparse([4, 8], [1, 2], [-1.0, -1.0]) # input matrix
    X = Universe(8) # state domain: unconstrained
    U = Hyperrectangle([0.23, 0.3], [0.07, 0.1]) # input domain
    X0 = initial_state_motor()

    prob = @ivp(x' = Ax + Bu, x ∈ X, u ∈ U) # continuous LTI system
    tspan = (0.0, 20.0)
    return prob, tspan
end
