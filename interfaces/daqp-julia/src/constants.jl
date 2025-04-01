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

# Exit Flags
const  SOFT_OPTIMAL   =  2
const  OPTIMAL        =  1
const  INFEASIBLE     = -1
const  CYCLE          = -2
const  UNBOUNDED      = -3
const  ITERLIMIT      = -4
const  NONCONVEX      = -5
const  OVERDETERMINED = -6

const flag2status= Dict{Int,Symbol}(2 => :Soft_Optimal,
                                    1 => :Optimal,
                                   -1 => :Primal_Infeasible,
                                   -2 => :Cycling,
                                   -3 => :Unbounded,
                                   -4 => :Iteration_Limit,
                                   -5 => :Nonconvex_Problem,
                                   -6 => :Initial_Overdetermined)
