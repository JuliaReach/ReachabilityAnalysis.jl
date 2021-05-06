# ===========================================================================
# Conservative linearization and range enclosures
# ===========================================================================

# ------------------------------------------
# Enclosing the range of 1/2 * v' * Hi * v
# ------------------------------------------

abstract type AbstractRangeEnclosureMethod end

struct IntervalArithmeticEnclosure <: AbstractRangeEnclosureMethod end

# range enclosure of the Hessian using IA
function _enclose_hessian(Hsym::Vector{<:AbstractMatrix}, var, X0::IntervalBox, c, alg::AbstractRangeEnclosureMethod=IntervalArithmeticEnclosure())
    d = var - c
    ex = [1/2 * d' * Hi * d for Hi in Hsym]
    dict = [vi => ci for (vi, ci) in zip(var, X0.v)]
    y = Symbolics.substitute.(ex, Ref(dict))
    return IntervalBox(Symbolics.value.(y))
end

# fallback for generic LazySet
function _enclose_hessian(Hsym::Vector{<:AbstractMatrix}, var, X0::LazySet, c, alg::AbstractRangeEnclosureMethod=IntervalArithmeticEnclosure())
    B0 = convert(IntervalBox, box_approximation(X0))
    return _enclose_hessian(Hsym, var, B0, c, alg)
end

# ------------------------------------------
# Conservative linearization
# ------------------------------------------

# _new_variables(n) = [Num(SymbolicUtils.Sym{Real}(:x, i)) for i in 1:n]

function _linearize(ivp::InitialValueProblem, var::Vector{Num},
                    alg::AbstractRangeEnclosureMethod=IntervalArithmeticEnclosure())

    n = statedim(ivp)
    #var = @variables x[1:n]
    #var = reduce(vcat, var)

    J = _jacobian_function(ivp, var)
    #H = _hessian_function(ivp, var)
    Hsym = _hessian_expression(ivp, var)

    # linearization center
    X0 = initial_state(ivp)
    c = center(X0)

    # evaluations
    #=
    ff = eval(build_function(f, var)[1]) # TODO generalize with VectorField?
    fc = ff(c)
    =#
    f = ivp.s.f
    cc = Dict(vi => ci for (vi, ci) in zip(var, c))
    fc = Symbolics.value.(Symbolics.substitute.(f, Ref(cc)))

    Jc = J(c)
    A = Jc
    b = fc - Jc*c
#=
    # input set
    B = Matrix(1.0I, n, n)
    UH = _enclose_hessian(Hsym, var, X0, c, alg)
    UH = convert(Hyperrectangle, UH)
    ivplin = @ivp(x' = A*x + b + B*u, x(0) ∈ X0, x ∈ Universe(n), u ∈ UH)

    return ivplin
    =#
end
