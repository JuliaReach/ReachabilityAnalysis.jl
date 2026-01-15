using .Symbolics: jacobian, hessian, build_function

# ===========================================================================
# Functionalities for Jacobian matrices
# ===========================================================================

# ----------------------------------------
# Symbolic expressions for the Jacobian
# ----------------------------------------

function _jacobian_expression(f::AbstractVector{Num}, var::Vector{Num})
    Jsym = jacobian(f, var)
    return Jsym
end

function _jacobian_expression(sys::BBCS{Vector{Num}}, var)
    field = sys.f
    return _jacobian_expression(field, var)
end

# fallback for nonlinear systems: assumes in-place form f!(dx, x, p, t),
# *ignoring* p and t args and uses `var` to trace the associated symbolic function
function _jacobian_expression(sys::BBCS, var)
    n = length(var) # dimension
    expr = Vector{Num}(undef, n)
    field! = sys.f
    field!(expr, var, [], [])
    return _jacobian_expression(expr, var)
end

# fallback implementation: extract the vector field and use the symbolic variables
# `var` to trace the associated symbolic vector field
function _jacobian_expression(sys::AbstractContinuousSystem, var)
    F = field(VectorField(sys))
    expr = F(var)
    return _jacobian_expression(expr, var)
end

function _jacobian_expression(ivp::IVP, var)
    return _jacobian_expression(ivp.s, var)
end

# ------------------------------------------
# Julia functions to evaluate the Jacobian
# ------------------------------------------

# given an ODE right-hand-side, compute the Jacobian function symbolically,
# then compile a Julia function for that symbolic expression
function _jacobian_function(f::AbstractVector{Num}, var::Vector{Num})

    # compute the symbolic jacobian
    Jsym = _jacobian_expression(f, var)

    # use the out-of-place version
    Jexpr = build_function(Jsym, var)[1]

    # return the associated Julia function
    return eval(Jexpr)
end

function _jacobian_function(sys::BBCS{Vector{Num}}, var)
    field = sys.f
    return _jacobian_function(field, var)
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

function _jacobian_function(ivp::IVP, var)
    return _jacobian_function(ivp.s, var)
end

# ===========================================================================
# Functionalities for Hessian matrices
# ===========================================================================

# ----------------------------------------
# Symbolic expressions for the Hessian
# ----------------------------------------

# return a vector of Hessian symbolic expressions, one for each component function
function _hessian_expression(f::AbstractVector{Num}, var::Vector{Num})
    return [_hessian_expression(fi, var) for fi in f]
end

function _hessian_expression(f::AbstractVector{Num}, i::Int, var::Vector{Num})
    return _hessian_expression(f[i], var)
end

function _hessian_expression(fi::Num, var::Vector{Num})
    Hsym = hessian(fi, var)
    return Hsym
end

function _hessian_expression(sys::BBCS{Vector{Num}}, var)
    field = sys.f
    return _hessian_expression(field, var)
end
function _hessian_expression(sys::BBCS, var)
    n = length(var) # dimension
    expr = Vector{Num}(undef, n)
    field! = sys.f
    field!(expr, var, [], [])
    return _hessian_expression(expr, var)
end

function _hessian_expression(sys::AbstractContinuousSystem, var)
    F = field(VectorField(sys))
    expr = F(var)
    return _hessian_expression(expr, var)
end

function _hessian_expression(ivp::IVP, var)
    return _hessian_expression(ivp.s, var)
end

# -----------------------------------------
# Julia functions to evaluate the Hessian
# -----------------------------------------

# return a vector of Hessian functions, one for each component function
function _hessian_function(f::AbstractVector{Num}, var::Vector{Num})
    return [_hessian_function(fi, var) for fi in f]
end

function _hessian_function(f::AbstractVector{Num}, i::Int, var::Vector{Num})
    return _hessian_function(f[i], var)
end

function _hessian_function(fi::Num, var::Vector{Num})

    # compute the symbolic jacobian
    Hsym = _hessian_expression(fi, var)

    # use the out-of-place version
    Hexpr = build_function(Hsym, var)[1]

    # return the associated Julia function
    return eval(Hexpr)
end

function _hessian_function(sys::BBCS{Vector{Num}}, var)
    field = sys.f
    return _hessian_function(field, var)
end

function _hessian_function(sys::BBCS, var)
    n = length(var) # dimension
    expr = Vector{Num}(undef, n)
    field! = sys.f
    field!(expr, var, [], [])
    return _hessian_function(expr, var)
end

function _hessian_function(sys::AbstractContinuousSystem, var)
    F = field(VectorField(sys))
    expr = F(var)
    return _hessian_function(expr, var)
end

function _hessian_function(ivp::IVP, var)
    return _hessian_function(ivp.s, var)
end
