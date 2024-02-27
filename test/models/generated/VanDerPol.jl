@taylorize function vanderpol!(du, u, p, t)
    x, y = u
    local μ = 1.0
    du[1] = y
    du[2] = (μ * y) * (1 - x^2) - x
    return du
end

function vanderpol()
    X0 = Hyperrectangle(; low=[1.25, 2.35], high=[1.55, 2.45])
    prob = @ivp(x' = vanderpol!(x), dim:2, x(0) ∈ X0)
    tspan = (0.0, 5.0)
    return prob, tspan
end

X0 = Hyperrectangle(; low=[1.25, 2.35], high=[1.55, 2.45])
prob = @ivp(x' = vanderpol!(x), dim:2, x(0) ∈ X0);

sol = solve(prob; T=7.0, alg=TMJets(; abstol=1e-12));

solz = overapproximate(sol, Zonotope);

y_bound = ρ([0.0, 1.0], solz)
@assert y_bound < 2.75 "the property should be proven"
y_bound

line = LineSegment([1, 2.0], [2.0, 2.5])

function cross_section(line::LineSegment, F, idx)
    p = VPolygon()
    for i in idx
        x = intersection(line, set(F[i]))
        if !isempty(x)
            p = convex_hull(p, x)
        end
    end
    vl = vertices_list(p)
    @assert length(vl) == 2
    return LineSegment(vl[1], vl[2])
end;

ifirst = cross_section(line, solz, 1:13)
ilast = cross_section(line, solz, [388]);

lfirst = norm(ifirst.q - ifirst.p)

llast = norm(ilast.q - ilast.p)

@assert ilast ⊆ ifirst "the cross section should get smaller"
