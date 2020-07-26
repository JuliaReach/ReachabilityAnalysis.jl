# ============================================
# Functionality that requires ModelingToolkit
# ============================================

#=
# outsourced to LazySets.jl
#
function load_modeling_toolkit_polyhedron()
return quote

function LazySets.HPolyhedron(expr::Vector{<:Operation}, vars=get_variables(first(expr)); N::Type{<:Real}=Float64)
    return HPolyhedron([HalfSpace(ex, vars) for ex in expr])
end

end end  # quote / load_modeling_toolkit_polyhedron()
=#
