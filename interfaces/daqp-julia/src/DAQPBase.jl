module DAQPBase

using DAQP_jll
# Prefer a locally-built libdaqp (e.g. placed there by a CMake build) over the
# DAQP_jll artifact.  This allows downstream packages that `Pkg.develop` this
# module to pick up the latest compiled library automatically.
let _local_lib = joinpath(@__DIR__, "..", "libdaqp." * Libc.Libdl.dlext)
    if isfile(_local_lib)
        @info "DAQPBase: using local libdaqp at $_local_lib"
        DAQP_jll.libdaqp = _local_lib
    end
end

using LinearAlgebra 

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
include("avi_julia.jl")

end
