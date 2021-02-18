using Reachability.ReachSets: tune_δ

#===== Projectile model =====
We test the line plot!(x->x, x->-24*x+375, 0., 50., line=2., color="red", linestyle=:solid, legend=:none)
=#
A = [0. 0.5 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.7 ; 0. 0. 0. 0.]
X0 = Singleton([0.,5.,100.,0])
U = ConstantInput(Singleton([0.,0.,0.,-9.81]))
S = IVP(CLCDS(A, Matrix(1.0I, 4, 4), nothing, U), X0)
time_horizon = 20.0
prec = 1e-4
initial_δ = 0.5
property = SafeStatesProperty(LinearConstraint([24., 0., 1, 0], 375.))
algorithm(N, δ) =
    solve(S,
          Options(:mode => "check", :plot_vars => [1, 3], :T => time_horizon,
                  :property => property),
          op=BFFPSV18(:δ => δ, :vars => [1, 3], :partition=>[1:2, 3:4])
         ).satisfied

tune_δ(algorithm, time_horizon, prec, initial_δ)
