module DAQPBase

using DAQP_jll

include("types.jl")
include("constants.jl")

include("api.jl")

export quadprog
export linprog

export setup
export solve
export update
export settings

export isfeasible

include("daqp_julia.jl")

end
