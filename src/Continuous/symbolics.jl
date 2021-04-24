using .Symbolics: jacobian, hessian

# ===========================================================================
# Functionalities for Jacobian matrices
# ===========================================================================

# given an ODE right-hand-side, compute the Jacobian function symbolically,
# then compile a Julia function for that symbolic exprsesion
function _jacobian_function(f::AbstractVector{Num}, var::Vector{Num})

    # compute the symbolic jacobian
    Jsym = jacobian(f, var)

    # use the out-of-place version
    Jexpr = build_function(Jsym, var)[1]

    # return the associated Julia function
    return eval(Jexpr)
end

function _jacobian_function(sys::BBCS{Vector{Num}}, var)
    field = sys.f
    return _jacobian_function(sys.f, var)
end

# fallback for nonlinear systems: assumes in-place form f!(dx, x, p, t),
# *ignoring* p and t args and uses `var` to trace the associated symbolic function
function _jacobian_function(sys::BBCS, var)
    n = length(var) # dimension
    expr = Vector{Num}(undef, n)
    field! = sys.f
    field!(expr, var, [], [])
    return _jacobian_function(expr, var)
end

# fallback implementation: extract the vector field and use the symbolic variables
# `var` to trace the associated symbolic vector field
function _jacobian_function(sys::AbstractContinuousSystem, var)
    F = field(VectorField(sys))
    expr = F(var)
    return _jacobian_function(expr, var)
end

function _jacobian_function(ivp::InitialValueProblem, var)
    return _jacobian_function(ivp.s, var)
end

# ===========================================================================
# Functionalities for Hessian matrices
# ===========================================================================

# given an ODE right-hand-side, compute the Hessian functions symbolically,
# then compile a Julia function for each symbolic exprsesion

# ===========================================================================
# Conservative linearization
# ===========================================================================
