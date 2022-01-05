module DAQP

const lib = "libdaqp.so"
include("types.jl")
include("constants.jl")

function quadprog(H::Matrix{Float64},f::Vector{Float64}, 
	A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64},sense::Vector{Int64})
  # TODO: check validity of dimensions
  # Setup QP
  (mA,n) = size(A);
  m = length(bupper);
  ms = m-mA;
  qp = QP(n,m,mA,H,f,A,bupper,blower,sense);
  # Setup settings
  settings = DAQPSettings();

  # Setup output struct
  xstar = zeros(Float64,n); 
  result= Ref(DAQPResult(xstar));

  ccall((:daqp_quadprog, DAQP.lib), Nothing,
		(Ref{DAQP.DAQPResult},Ref{DAQP.QP},Ref{DAQP.DAQPSettings}), 
		result,Ref(qp),Ref(settings))
  
  info = (status = DAQP.flag2status[result[].exitflag],
		  solve_time = result[].solve_time,
		  setup_time = result[].setup_time,
		  inner_iter = result[].iter,
		  outer_iter = result[].outer_iter)
  return xstar,result[].fval,result[].exitflag,info
end

end
