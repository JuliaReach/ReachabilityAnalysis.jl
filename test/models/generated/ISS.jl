using ReachabilityAnalysis, JLD2, ReachabilityBase.CurrentPath
using ReachabilityAnalysis: add_dimension

path = @current_path("ISS", "ISS.jld2")

@load path C
const C3 = C[3, :]  # variable y₃
const C3_ext = vcat(C3, fill(0.0, 3));

function ISS_common(n)
    X0 = BallInf(zeros(n), 0.0001)
    U = Hyperrectangle(; low=[0.0, 0.8, 0.9], high=[0.1, 1.0, 1.0])
    return X0, U
end

function ISSF01()
    @load path A B
    X0, U = ISS_common(size(A, 1))
    return @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

function ISSC01()
    @load path A B
    A_ext = add_dimension(A, 3)
    A_ext[1:270, 271:273] = B
    X0, U = ISS_common(size(A, 1))
    return @ivp(x' = A_ext * x, x(0) ∈ X0 * U)
end;

dirs = CustomDirections([C3, -C3])
prob_ISSF01 = ISSF01()
sol_ISSF01 = solve(prob_ISSF01; T=20.0,
                   alg=LGG09(; δ=6e-4, template=dirs, sparse=true, cache=false));

dirs = CustomDirections([C3_ext, -C3_ext])
prob_ISSC01 = ISSC01()
sol_ISSC01 = solve(prob_ISSC01; T=20.0,
                   alg=LGG09(; δ=0.01, template=dirs, sparse=true, cache=false));
