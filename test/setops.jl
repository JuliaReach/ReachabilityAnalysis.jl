@testset "Intersection methods: TemplateHullIntersection" begin

    # default constructor without types
    @test TemplateHullIntersection() ==
          TemplateHullIntersection{Float64,Vector{Float64},Missing,Val{false}}(missing,
                                                                               Val{false}())

    # default constructor passing normal vector types
    @test TemplateHullIntersection{Float64,Vector{Float64}}() ==
          TemplateHullIntersection{Float64,Vector{Float64},Missing,Val{false}}(missing,
                                                                               Val{false}())
    @test TemplateHullIntersection{Float64,SVector{2,Float64}}() ==
          TemplateHullIntersection{Float64,SArray{Tuple{2},Float64,1,2},Missing,Val{false}}(missing,
                                                                                            Val{false}())

    # constructor with a given template
    td = TemplateHullIntersection(BoxDirections(2))
    @test setrep(td) == Union{HPolytope{Float64,LazySets.Arrays.SingleEntryVector{Float64}},
                              HPolytope{Float64,Vector{Float64}}}
end
