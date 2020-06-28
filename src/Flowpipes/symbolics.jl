# ============================================
# Functionality that requires ModelingToolkit
# ============================================
function load_modeling_toolkit_hyperplane()
return quote

"""
    Hyperplane(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)

Return the hyperplane given by a symbolic expression.

### Input

- `expr` -- symbolic expression that describes a hyperplane
- `vars` -- (optional, default: `get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may change; see `Notes` below fo details)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned half-space

### Output

A `Hyperplane`.

### Examples

```julia
julia> using ModelingToolkit

julia> vars = @variables x y
(x, y)

julia> Hyperplane(x - y == 2)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0], 2.0)

julia> Hyperplane(x == y)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0], -0.0)

julia> vars = @variables x[1:4]
(Operation[x₁, x₂, x₃, x₄],)

julia> Hyperplane(x[1] == x[2], x)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0, 0.0, 0.0], -0.0)
```

### Algorithm

It is assumed that the expression is of the form
`EXPR0: α*x1 + ⋯ + α*xn + γ CMP β*x1 + ⋯ + β*xn + δ`,
where `CMP` is `==`.
This expression is transformed, by rearrangement and substitution, into the
canonical form `EXPR1 : a1 * x1 + ⋯ + an * xn == b`. The method used to identify
the coefficients is to take derivatives with respect to the ambient variables `vars`.
Therefore, the order in which the variables appear in `vars` affects the final result.
Finally, the returned set is the hyperplane with normal vector `[a1, …, an]` and
displacement `b`.
"""
function Hyperplane(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)
    (expr.op == ==) || throw(ArgumentError("expected an expression of the form `ax == b`, got $expr"))

    # simplify to the form a*x + β == 0
    a, b = expr.args
    sexpr = simplify(a - b)

    # compute the linear coefficients by taking first order derivatives
    coeffs = [N(α.value) for α in gradient(sexpr, collect(vars))]

    # get the constant term by expression substitution
    dvars = Dict(to_symbolic(vi) => zero(N) for vi in vars)
    β = -N(ModelingToolkit.SymbolicUtils.substitute(to_symbolic(sexpr), dvars, fold=true))

    return Hyperplane(coeffs, β)
end

end end  # quote / load_modeling_toolkit_hyperplane()

function load_modeling_toolkit_halfspace()
return quote

"""
    HalfSpace(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)

Return the half-space given by a symbolic expression.

### Input

- `expr` -- symbolic expression that describes a half-space
- `vars` -- (optional, default: `get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may change; see `Notes` below fo details)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned half-space

### Output

A `HalfSpace`.

### Examples

```julia
julia> using ModelingToolkit

julia> vars = @variables x y
(x, y)

julia> HalfSpace(x - y <= 2)
HalfSpace{Float64,Array{Float64,1}}([1.0, -1.0], 2.0)

julia> HalfSpace(x >= y)
HalfSpace{Float64,Array{Float64,1}}([1.0, -1.0], -0.0)

julia> vars = @variables x[1:4]
(Operation[x₁, x₂, x₃, x₄],)

julia> HalfSpace(x[1] >= x[2], x)
HalfSpace{Float64,Array{Float64,1}}([-1.0, 1.0, 0.0, 0.0], -0.0)
```

Be careful with using the default `vars` value, becaus it may introduce a wrong
order.

```julia
julia> vars = @variables x y
(x, y)

julia> HalfSpace(2x ≥ 5y - 1) # wrong
HalfSpace{Float64,Array{Float64,1}}([5.0, -2.0], 1.0)

julia> HalfSpace(2x ≥ 5y - 1, vars) # correct
HalfSpace{Float64,Array{Float64,1}}([-2.0, 5.0], 1.0)
```

### Algorithm

It is assumed that the expression is of the form
`EXPR0: α*x1 + ⋯ + α*xn + γ CMP β*x1 + ⋯ + β*xn + δ`,
where `CMP` is one among `<`, `<=`, `≤`, `>`, `>=` or `≥`.
This expression is transformed, by rearrangement and substitution, into the
canonical form `EXPR1 : a1 * x1 + ⋯ + an * xn ≤ b`. The method used to identify
the coefficients is to take derivatives with respect to the ambient variables `vars`.
Therefore, the order in which the variables appear in `vars` affects the final result.
Note in particular that strict inequalities are relaxed as being smaller-or-equal.
Finally, the returned set is the half-space with normal vector `[a1, …, an]` and
displacement `b`.
"""
function HalfSpace(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)

    # find sense and normalize
    if expr.op == <
        a, b = expr.args
        sexpr = simplify(a - b)

    elseif expr.op == >
        a, b = expr.args
        sexpr = simplify(b - a)

    elseif (expr.op == |) && (expr.args[1].op == <)
        a, b = expr.args[1].args
        sexpr = simplify(a - b)

    elseif (expr.op == |) && (expr.args[2].op == <)
        a, b = expr.args[2].args
        sexpr = simplify(a - b)

    elseif (expr.op == |) && (expr.args[1].op == >)
        a, b = expr.args[1].args
        sexpr = simplify(b - a)

    elseif (expr.op == |) && (expr.args[2].op == >)
        a, b = expr.args[2].args
        sexpr = simplify(b - a)

    else
        throw(ArgumentError("expected an expression describing a half-space, got $expr"))
    end

    # compute the linear coefficients by taking first order derivatives
    coeffs = [N(α.value) for α in gradient(sexpr, collect(vars))]

    # get the constant term by expression substitution
    dvars = Dict(to_symbolic(vi) => zero(N) for vi in vars)
    β = -N(ModelingToolkit.SymbolicUtils.substitute(to_symbolic(sexpr), dvars, fold=true))

    return HalfSpace(coeffs, β)
end

end end  # quote / load_modeling_toolkit_halfspace()
