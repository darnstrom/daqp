# Constraint types
@enum SelectionRule begin
    DANTZIG = 0
    BLAND
end


const ACTIVE = 1
const LOWER = 2
const IMMUTABLE= 4
const EQUALITY = ACTIVE+IMMUTABLE
const SOFT= 8
const BINARY= 16

# LDP update masks
const DAQP_UPDATE_Rinv = 1
const DAQP_UPDATE_M = 2
const DAQP_UPDATE_v = 4
const DAQP_UPDATE_d = 8
const DAQP_UPDATE_sense = 16
const DAQP_UPDATE_hierarchy = 32
const DAQP_UPDATE_unconstrained = 64
const DAQP_UPDATE_eliminate = 128

# Exit Flags
const  CONSTRAINED_POINT   = 3
const  SOFT_OPTIMAL   =  2
const  OPTIMAL        =  1
const  INFEASIBLE     = -1
const  CYCLE          = -2
const  UNBOUNDED      = -3
const  ITERLIMIT      = -4
const  NONCONVEX      = -5
const  OVERDETERMINED = -6
const  TIMELIMIT      = -7
const  UNSUPPORTED    = -8

const flag2status= Dict{Int,Symbol}(3 => :Constrained_Point,
                                    2 => :Soft_Optimal,
                                    1 => :Optimal,
                                   -1 => :Primal_Infeasible,
                                   -2 => :Cycling,
                                   -3 => :Unbounded,
                                   -4 => :Iteration_Limit,
                                   -5 => :Nonconvex_Problem,
                                   -6 => :Initial_Overdetermined,
                                   -7 => :Time_Limit,
                                   -8 => :Unsupported_Problem)
