using ReachabilityAnalysis, Test
import Aqua

@testset "Aqua tests" begin
    # note that the check `persistent_tasks` may take several minutes
    Aqua.test_all(ReachabilityAnalysis; ambiguities=false,
                  # the unbound args should be resolved in the future
                  unbound_args=(broken=true,),
                  # the piracies should be resolved in the future
                  piracies=(broken=true,))

    # do not warn about ambiguities in dependencies
    # the ambiguities should be resolved in the future
    Aqua.test_ambiguities(ReachabilityAnalysis; broken=true)
end
