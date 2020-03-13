# See MathematicalSystems#160
MathematicalSystems.system(sys::AbstractSystem) = sys
#MathematicalSystems.system(sys::InitialValueProblem) = sys.s

struct VectorField{T}
    field::T
end

# function-like evaluation
@inline function (V::VectorField)(args...)
    evaluate(V, args...)
end

function evaluate(V::VectorField, args...)
    return V.field(args...)
end

function VectorField(sys::AbstractSystem)
    sys = system(sys)
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
