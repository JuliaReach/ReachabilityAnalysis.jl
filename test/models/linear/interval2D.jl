# ============================================================================
# [AlthoffSB07; Example 1](@citet) TODO: Agregar referencia completa
# ============================================================================

# initial set
X0 = BallInf([1.0, 1.0], 0.1)

# linear ODE: x' = Ax
A = IntervalMatrix([interval(-1.05, -0.95) interval(-4.05, -3.95);
                    interval(3.95, 4.05) interval(-1.05, -0.95)])

# IVP(LCS(A), X0) TODO: remove
interval2D_linear = @ivp x' = A * x, x(0) ∈ X0

# affine ODE: x' = Ax + Bu
B = IntervalMatrix(hcat([interval(1); interval(1)])) # TODO why hcat?
U = Interval(-0.05, 0.05)

# IVP(CLCCS(A, B, nothing, U), X0) TODO: remove
interval2D_affine = @ivp x' = A * x + Bu, x(0) ∈ X0, u ∈ U
