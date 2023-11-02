abstract type AbstractVectorField end

MathematicalSystems.system(sys::AbstractContinuousSystem) = sys

struct VectorField{T} <: AbstractVectorField
    field::T
end

field(V::VectorField) = V.field

# no-op (for nonlinear systems)
field(f::Function) = f

# function-like evaluation
@inline function (V::VectorField)(args...)
    return evaluate(V, args...)
end

@inline function evaluate(V::VectorField, args...)
    return V.field(args...)
end

@inline function evaluate(V::VectorField, x::Number, args...)
    return V.field([x], args...)
end

VectorField(sys::InitialValueProblem) = VectorField(system(sys))

function VectorField(sys::AbstractContinuousSystem)
    if islinear(sys)
        if inputdim(sys) == 0
            field = (x) -> state_matrix(sys) * x
        else
            field = (x, u) -> state_matrix(sys) * x + input_matrix(sys) * u
        end
    elseif isaffine(sys)
        if inputdim(sys) == 0
            field = (x) -> state_matrix(sys) * x + affine_term(sys)
        else
            field = (x, u) -> state_matrix(sys) * x + affine_term(sys) + input_matrix(sys) * u
        end
    elseif sys isa BBCS || sys isa CBBCS || sys isa CBBCCS
        return sys.f
    else
        error("the vector field for a system of type $sys is not implemented yet")
    end

    return VectorField(field)
end

function vector_field(sys::AbstractSystem, args...)
    return evaluate(VectorField(sys), args...)
end

# (x, p, t) -> du
function outofplace_field(ivp::InitialValueProblem)
    vf = VectorField(ivp)

    # function closure over the initial-value problem
    f = function f_outofplace(x, p, t)
        return vf(x)
    end
    return f
end

# f!(dx, x, p, t)
function inplace_field!(ivp::InitialValueProblem)
    vf = VectorField(ivp)

    # function closure over the initial-value problem
    f! = function f_inplace!(dx, x, p, t)
        return dx .= vf(x)
    end
    return f!
end

# for "black-box" MathematicalSystems, return the function field directly

function inplace_field!(ivp::IVP{<:NonlinearSystem})
    return ivp.s.f # TODO getter?
end
