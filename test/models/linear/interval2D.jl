# ============================================================================
# [AlthoffSB07; Example 1](@citet) TODO: Agregar referencia completa
# ============================================================================

# initial set
X0 = BallInf([1.0, 1.0], 0.1)

# linear ODE: x' = Ax
A = IntervalMatrix([-1.0±0.05 -4.0±0.05;
                    4.0±0.05 -1.0±0.05])

# IVP(LCS(A), X0) TODO: remove
interval2D_linear = @ivp x' = Ax, x(0) ∈ X0

# affine ODE: x' = Ax + Bu
B = IntervalMatrix(hcat([1.0 ± 0.0; 1.0 ± 0.0])) # why hcat?
U = Interval(-0.05, 0.05)

# IVP(CLCCS(A, B, nothing, U), X0) TODO: remove
interval2D_affine = @ivp x' = Ax + Bu, x(0) ∈ X0, u ∈ U
