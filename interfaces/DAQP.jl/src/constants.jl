const flag2status= Dict{Int,Symbol}(1 => :Optimal,
								   -1 => :Primal_Infeasible,
								   -2 => :Cycling,
								   -3 => :Unbounded,
								   -4 => :Iteration_Limit
								   )
