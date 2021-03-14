# flowpipe interface
include("AbstractFlowpipe.jl")

# concrete implementations
include("Flowpipe.jl")
include("ShiftedFlowpipe.jl")
include("MappedFlowpipe.jl")
include("HybridFlowpipe.jl")
include("MixedFlowpipe.jl")
include("PartitionedFlowpipe.jl")
