using ReachabilityAnalysis, Symbolics

@variables x y z
positive_orthant = HPolyhedron([x >= 0, y >= 0, z >= 0], [x, y, z]);

function prod_dest_verif(sol; T=100.0, target=10.0)
    # obtain the final reach set and project onto the state variables x, y, z
    Xh = box_approximation(sol(T))
    X = project(Xh; vars=(1, 2, 3))

    # check that all variables are nonnegative
    nonnegative = X ⊆ positive_orthant

    # compute the volume
    vol = volume(X)

    # check that the target belongs to the Minkowski sum of the projections
    # onto each coordinate
    B = convert(IntervalBox, X) # get a product-of-intervals representation
    contains_target = target ∈ sum(B)

    return nonnegative && contains_target, vol
end;

@taylorize function prod_dest_I!(du, u, p, t)
    local a = 0.3
    x, y, z = u

    xyx = (x * y) / (1 + x)
    ay = a * y

    du[1] = -xyx
    du[2] = xyx - ay
    du[3] = ay
    return du
end

X0 = Hyperrectangle(; low=[9.5, 0.01, 0.01], high=[10, 0.01, 0.01])
prob = @ivp(x' = prod_dest_I!(x), dim:3, x(0) ∈ X0);

sol = solve(prob; T=100.0, alg=TMJets(; abstol=1e-12, orderT=6, orderQ=1));

property, vol = prod_dest_verif(sol)
@assert property "the property should be proven"

@taylorize function prod_dest_IP!(du, u, p, t)
    x, y, z, a = u[1], u[2], u[3], u[4]

    xyx = (x * y) / (1 + x)
    ay = a * y

    du[1] = -xyx
    du[2] = xyx - ay
    du[3] = ay
    du[4] = zero(a)
    return du
end

X0 = Hyperrectangle(; low=[9.98, 0.01, 0.01, 0.296], high=[9.98, 0.01, 0.01, 0.304])
prob = @ivp(x' = prod_dest_IP!(x), dim:4, x(0) ∈ X0);

sol = solve(prob; T=100.0, alg=TMJets(; abstol=9e-13, orderT=6, orderQ=1));

property, vol = prod_dest_verif(sol)
@assert property "the property should be proven"

X0 = Hyperrectangle(; low=[9.5, 0.01, 0.01, 0.296], high=[10, 0.01, 0.01, 0.304])
prob = @ivp(x' = prod_dest_IP!(x), dim:4, x(0) ∈ X0);

sol = solve(prob; T=100.0, alg=TMJets(; abstol=1e-12, orderT=6, orderQ=1));

property, vol = prod_dest_verif(sol)
@assert property "the property should be proven"
