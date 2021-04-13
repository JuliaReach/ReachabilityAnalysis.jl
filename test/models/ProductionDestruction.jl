using ReachabilityAnalysis, Symbolics, Plots

@variables x y z
const positive_orthant = HPolyhedron([x >= 0, y >= 0, z >= 0], [x, y, z])

function prod_dest_verif(sol; T=100.0, target=10.0)

    # convert to a zonotopic flowpipe
    solz = overapproximate(sol, Zonotope)

    # project the final reach-set onto the space variables x, y, z
    X = project(solz(T), vars=(1, 2, 3))

    # check that all variables are nonnegative
    nonnegative = X ⊆ positive_orthant

    # compute the volume of the last reach-set
    H = overapproximate(X, Hyperrectangle)
    vol = volume(H)

    # check that that target belongs to the minkowski sum of the reach-sets projected in each coordinate
    B = convert(IntervalBox, H) # get the product-of-intervals representation
    contains_target = target ∈ sum(B)

    return nonnegative && contains_target, vol
end

@taylorize function prod_dest_I!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    du[1] = - (x * y) / (1 + x)
    du[2] = (x * y) / (1 + x) - a * y
    du[3] = a * y
    return du
end

X0 = (9.5 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01)
prob = @ivp(x'= prod_dest_I!(x), dim:3, x(0) ∈ X0)

solI = solve(prob, T=100.0, alg=TMJets(abs_tol=1e-11, orderT=7, orderQ=1));

property, vol = prod_dest_verif(solI)

plot(solI, vars=(0, 3), linecolor=:orange, color=:orange, alpha=0.3, lab="I")

@taylorize function prod_dest_IP!(du, u, params, t)
    x, y, z, a = u[1], u[2], u[3], u[4]

    du[1] = - (x * y) / (1 + x)
    du[2] = (x * y) / (1 + x) - a * y
    du[3] = a * y
    du[4] = zero(x)
    return du
end

X0 = (9.98 .. 9.98) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.296 .. 0.304)
prob = @ivp(x'= prod_dest_IP!(x), dim:4, x(0) ∈ X0)

solP = solve(prob, T=100.0, alg=TMJets(abs_tol=1e-12, orderT=7, orderQ=1));

property, vol = prod_dest_verif(solP)

plot(solP, vars=(0, 3), linecolor=:blue, color=:blue, alpha=0.3, lab="P")

X0 = (9.5 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.296 .. 0.304)
prob = @ivp(x'= prod_dest_IP!(x), dim:4, x(0) ∈ X0)

solIP = solve(prob, T=100.0, alg=TMJets(abs_tol=1e-11, orderT=7, orderQ=1));

property, vol = prod_dest_verif(solIP)

plot(solIP, vars=(0, 3), linecolor=:red, color=:red, alpha=0.3, lab="I & P")

@taylorize function prod_dest_I_optimized!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    return du
end

@taylorize function prod_dest_IP_optimized!(du, u, params, t)
    x, y, z, a = u[1], u[2], u[3], u[4]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    du[4] = zero(x)
    return du
end

