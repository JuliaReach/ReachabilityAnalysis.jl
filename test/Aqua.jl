using ReachabilityAnalysis, Test
import Aqua

@testset "Aqua tests" begin
    # note that the check `persistent_tasks` may take several minutes
    Aqua.test_all(ReachabilityAnalysis;
                  # the piracies should be resolved in the future
                  piracies=(broken=true,),
                  # the test for persistent tasks started failing in 2026
                  persistent_tasks=(broken=true,))
end
