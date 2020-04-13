ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, ReachabilityAnalysis

DocMeta.setdocmeta!(ReachabilityAnalysis, :DocTestSetup,
                   :(using ReachabilityAnalysis); recursive=true)

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS")), # disable for local builds
    sitename = "ReachabilityAnalysis.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md", # TODO: keep?
        "Manual" => Any["Set representations" => "man/setrep.md",
                        "Linear ODEs" => "man/linear.md",
                        "Exploiting structure" => "man/linear_high_dim.md",
                        "Nonlinear ODEs" => "man/nonlinear.md",
                        "Hybrid systems" => "man/hybrid.md",
                        "Parametric reachability" => "man/parametric.md",
                        #"FAQ" => "man/faq.md",
                        "Benchmarks" => "man/benchmarks.md",
                        "Model library" => "man/library.md"],
                        # Other topics: Distributed computations. Multithreading.
        "Applications" => Any["Electromechanic break" => "man/applications/embrake.md",
                              "Quadrotor altitude control" => "man/applications/quadrotor.md",
                              "Transmision line" => "man/applications/transmission_line.md",
                              "Epidemic disease" => "man/applications/epidemic.md"],
                              # Other topics: car control, power systems stability.
        "Algorithms" => Any["LGG09" => "lib/algorithms/LGG09.md",
                            "GLGM06" => "lib/algorithms/GLGM06.md",
                            "INT" => "lib/algorithms/INT.md",
                            "BOX" => "lib/algorithms/BOX.md",
                            "BFFPSV18" => "lib/algorithms/BFFPSV18.md",
                            "TMJets" => "lib/algorithms/TMJets.md"],
                            #
        "API Reference" => Any["Reach-sets" => "lib/reachsets.md",
                               "Flowpipes" => "lib/flowpipes.md",
                               "Solution types" => "lib/solution_types.md",
                               "Discretization" => "lib/discretize.md"],
        "References" => "references.md",
        "About" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/JuliaReach/ReachabilityAnalysis.jl.git",
    push_preview=true,
)
