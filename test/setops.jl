@testset "Intersection strategies" begin
    # default constructors
    @test TemplateHullIntersection() == TemplateHullIntersection{Float64,Array{Float64,1},Missing}(missing)
    @test TemplateHullIntersection{Float64, Vector{Float64}}() == TemplateHullIntersection{Float64,Array{Float64,1},Missing}(missing)
    @test TemplateHullIntersection{Float64, SVector{2, Float64}}() == TemplateHullIntersection{Float64,SArray{Tuple{2},Float64,1,2},Missing}(missing)

    # constructor passing the template
    td = TemplateHullIntersection(BoxDirections(2))
    @test setrep(td) == HPolytope{Float64,LazySets.Arrays.SingleEntryVector{Float64}}
end
