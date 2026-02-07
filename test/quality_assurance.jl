using ReachabilityAnalysis, Test
import Aqua, ExplicitImports

@testset "ExplicitImports tests" begin
    ignores = (:AbstractDirections, :AbstractReductionMethod, :DEFAULT_ALPHA,
               :DEFAULT_ASPECT_RATIO, :DEFAULT_COLOR, :DEFAULT_GRID,
               :DEFAULT_LABEL, :PLOT_POLAR_DIRECTIONS, :PLOT_PRECISION,
               :_check_constrained_dimensions, :_expmv, :_plot_singleton_list,
               :plot_recipe, :shrink_wrapping!, :solve, :_exp_remainder,
               :correction_hull, :input_correction, :AbstractODEScheme,
               :AbstractPreconditioner, :AbstractTMOrder, :GIR05, :Slice,
               :_DEF_MINABSTOL)
    @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(ReachabilityAnalysis;
                                                                          ignore=ignores))
    @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(ReachabilityAnalysis))
    ignores = (:GIR05, :Slice, :_DEF_MINABSTOL, :_default_sampler, :GLOBAL_RNG,
               :value)
    @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(ReachabilityAnalysis;
                                                                            ignore=ignores))
    @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(ReachabilityAnalysis))
    # due to reexporting multiple packages
    ignores = (:HybridSystems, # due to reexporting HybridSystems:
               :AbstractHybridSystem, :AutonomousSwitching, :GraphAutomaton,
               :add_transition!, :guard, :mode, :nmodes, :nstates,
               :out_transitions, :resetmap, :source, :states, :target,
               # due to reexporting IntervalMatrices:
               :IntervalMatrices, :(..), :IntervalBox, :IntervalMatrix,
               :IntervalMatrixPower, :emptyinterval, :increment!, :inf,
               :interval, :mid, :mince, :pow, :sup,
               :IntervalArithmetic,  # superficially exported from IntervalMatrices
               # due to reexporting ReachabilityBase:
               :ReachabilityBase,
               # due to reexporting LazySets:
               :LazySets, :AbstractHyperrectangle, :AbstractPolyhedron,
               :AbstractPolytope, :AbstractSingleton, :AbstractZonotope,
               :Approximations, :BallInf, :BoxDirections, :CartesianProduct,
               :CartesianProductArray, :ConvexHull, :ConvexHullArray,
               :CustomDirections, :EmptySet, :ExponentialMap, :HPolygon,
               :HPolyhedron, :HPolytope, :HalfSpace, :Hyperplane,
               :Hyperrectangle, :Intersection, :LazySet, :LinearConstraint,
               :MatrixZonotope, :MatrixZonotopeExp, :MinkowskiSumArray,
               :OctDirections, :Singleton, :SparsePolynomialZonotope,
               :UnionSetArray, :Universe, :VPolygon, :VPolytope, :ZeroSet,
               :Zonotope, :affine_map, :area, :center, :complement, :concretize,
               :constrained_dimensions, :constraints_list, :convex_hull,
               :decompose, :element, :genmat, :high, :isbounded, :isbounding,
               :low, :minkowski_sum, :ngens, :order, :radius_hyperrectangle,
               :reduce_order, :remove_redundant_constraints!,
               :remove_redundant_generators, :symmetric_interval_hull, :tohrep,
               :translate, :vertices_list, :volume, :×, :ρ, :σ, :⊕,
               :Arrays,  # superficially exported from LazySets
               # due to reexporting MathematicalSystems:
               :MathematicalSystems, Symbol("@ivp"), Symbol("@map"),
               Symbol("@system"), :AbstractContinuousSystem,
               :AbstractDiscreteSystem, :AbstractInput, :AbstractMap,
               :AbstractSystem, :AffineContinuousSystem, :AffineDiscreteSystem,
               :BlackBoxContinuousSystem, :ConstantInput,
               :ConstrainedAffineContinuousSystem,
               :ConstrainedAffineControlContinuousSystem,
               :ConstrainedAffineControlDiscreteSystem,
               :ConstrainedAffineDiscreteSystem, :ConstrainedAffineMap,
               :ConstrainedBlackBoxContinuousSystem,
               :ConstrainedBlackBoxControlContinuousSystem,
               :ConstrainedIdentityMap, :ConstrainedLinearContinuousSystem,
               :ConstrainedLinearControlContinuousSystem,
               :ConstrainedLinearControlDiscreteSystem,
               :ConstrainedLinearDiscreteSystem, :ConstrainedLinearMap,
               :ConstrainedResetMap, :IVP, :Id, :IdentityMap, :IdentityMultiple,
               :LinearContinuousSystem, :LinearDiscreteSystem,
               :LinearParametricContinuousSystem,
               :LinearParametricDiscreteSystem,
               :SecondOrderAffineContinuousSystem,
               :SecondOrderConstrainedAffineControlContinuousSystem,
               :SecondOrderConstrainedLinearControlContinuousSystem,
               :SecondOrderLinearContinuousSystem, :VaryingInput, :affine_term,
               :input_matrix, :inputdim, :inputset, :isaffine, :islinear,
               :mass_matrix, :state_matrix, :stateset, :stiffness_matrix,
               :viscosity_matrix,
               # due to reexporting TaylorIntegration:
               :TaylorIntegration, :Taylor1, :get_numvars, :get_variables,
               :set_variables)
    @test isnothing(ExplicitImports.check_no_implicit_imports(ReachabilityAnalysis; ignore=ignores))
    @test isnothing(ExplicitImports.check_no_self_qualified_accesses(ReachabilityAnalysis))
    ignores = (:AbstractODEScheme, Symbol("@polyvar"))
    @test isnothing(ExplicitImports.check_no_stale_explicit_imports(ReachabilityAnalysis;
                                                                    ignore=ignores))
end

@testset "Aqua tests" begin
    # note that the check `persistent_tasks` may take several minutes
    Aqua.test_all(ReachabilityAnalysis;
                  # the piracies should be resolved in the future
                  piracies=(broken=true,),
                  # persistent tasks started failing occasionally on CI runs in January 2026
                  persistent_tasks=false)
end
