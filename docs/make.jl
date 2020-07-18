ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, ReachabilityAnalysis

DocMeta.setdocmeta!(ReachabilityAnalysis, :DocTestSetup,
                   :(using ReachabilityAnalysis); recursive=true)

# generate Literate documentation
include("generate.jl")

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),  # disable for local builds
                             collapselevel = 1,
                             assets = ["assets/juliareach-light.css"]),
    sitename = "ReachabilityAnalysis.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
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
        "Examples" => Any["Building" => "models/building.md",
                          "Van der Pol oscillator" => "models/vanderpol.md",
                          "Brusselator" => "models/brusselator.md",
                          "Platoon" => "models/platoon.md",
                          "Quadrotor" => "models/quadrotor.md",
                          "SEIR model" => "models/seir.md",
                          "Spacecraft" => "models/spacecraft.md",
                          "Transmision line" => "models/transmission_line.md",
                          "Platoon" => "models/platoon.md",
                          "Spacecraft" => "models/spacecraft.md",
                          "Transmision line" => "models/transmission_line.md",
                          "Laub-Loomis" => "models/laub_loomis.md"],
                          # "Electromechanic break" => "man/applications/embrake.md"
        "Algorithms" => Any["ASB07" => "lib/algorithms/ASB07.md",
                            "BFFPSV18" => "lib/algorithms/BFFPSV18.md",
                            "BOX" => "lib/algorithms/BOX.md",
                            "GLGM06" => "lib/algorithms/GLGM06.md",
                            "INT" => "lib/algorithms/INT.md",
                            "LGG09" => "lib/algorithms/LGG09.md",
                            "TMJets" => "lib/algorithms/TMJets.md"],
                            #
        "API Reference" => Any["Reach-sets" => "lib/reachsets.md",
                               "Flowpipes" => "lib/flowpipes.md",
                               "Solutions" => "lib/solutions.md",
                               "Discretization" => "lib/discretize.md"],
        "References" => "references.md",
        "About" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/JuliaReach/ReachabilityAnalysis.jl.git",
    push_preview = true,
)
