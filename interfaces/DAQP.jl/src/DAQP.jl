module DAQP

const lib = "libdaqp.so"
include("types.jl")
include("constants.jl")

function quadprog(H::Matrix{Float64},f::Vector{Float64}, 
	A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64},sense::Vector{Cint})
  return quadprog(QPj(H,f,A,bupper,blower,sense))
end
function quadprog(qpj::QPj)
  # TODO: check validity of dimensions
  # Setup QP
  qp = QPc(qpj);

  # Setup output struct
  xstar = zeros(Float64,qp.n); 
  result= Ref(DAQPResult(xstar));

  ccall((:daqp_quadprog, DAQP.lib), Nothing,
		(Ref{DAQP.DAQPResult},Ref{DAQP.QPc},Ref{DAQP.DAQPSettings}), 
		result,Ref(qp),Ptr{DAQP.DAQPSettings}(C_NULL))
  
  info = (status = DAQP.flag2status[result[].exitflag],
		  solve_time = result[].solve_time,
		  setup_time = result[].setup_time,
		  inner_iter = result[].iter,
		  outer_iter = result[].outer_iter)
  return xstar,result[].fval,result[].exitflag,info
end

mutable struct Model 
  work::Ptr{DAQP.Workspace}
  qpj::QPj
  qpc::QPc
  function Model()
	# Setup initial model
	work = Libc.calloc(1,sizeof(DAQP.Workspace))
	daqp= new(Ptr{DAQP.Workspace}(work))
	# Add default settings
	#settings = Ptr{DAQPSettings}(Libc.calloc(1,sizeof(DAQP.DAQPSettings)))
	#ccall((:daqp_default_settings,DAQP.lib),Nothing,(Ref{DAQP.DAQPSettings},),settings)
	#offset_settings = fieldoffset(DAQP.Workspace,findfirst(fieldnames(DAQP.Workspace).==:settings))
	#unsafe_store!(Ptr{Ptr{DAQP.DAQPSettings}}(work+offset_settings),settings)
	# Add finalizer to free memory
	finalizer(DAQP.delete!, daqp)
	return daqp 
  end
end
function delete!(daqp::DAQP.Model)
  workspace = unsafe_load(daqp.work);
  ccall((:free_daqp_workspace,DAQP.lib),Nothing,(Ptr{DAQP.Workspace},),daqp.work)
  Libc.free(daqp.work);
end

function setup(daqp::DAQP.Model, qp::DAQP.QPj)
  daqp.qpj = qp
  daqp.qpc = DAQP.QPc(daqp.qpj);
  return ccall((:setup_daqp,DAQP.lib),Cint,(Ref{DAQP.QPc}, Ptr{DAQP.Workspace}),
			   Ref{DAQP.QPc}(daqp.qpc), daqp.work)

end

function setup(daqp::DAQP.Model, H::Matrix{Cdouble},f::Vector{Cdouble},A::Matrix{Cdouble},bupper::Vector{Cdouble},blower::Vector{Cdouble},sense::Vector{Cint})
  setup(daqp,QPj(H,f,A,bupper,blower,sense))
end

function solve(daqp::DAQP.Model)
  xstar = zeros(Float64,daqp.qpc.n); 
  result= Ref(DAQPResult(xstar));

  exitflag=ccall((:daqp_solve, DAQP.lib), Cint,
		(Ref{DAQP.DAQPResult},Ref{DAQP.Workspace}), 
		result,daqp.work)
  
  info = (status = DAQP.flag2status[result[].exitflag],
		  solve_time = result[].solve_time,
		  setup_time = result[].setup_time,
		  inner_iter = result[].iter,
		  outer_iter = result[].outer_iter)
  return xstar,result[].fval,result[].exitflag,info
end

function settings(daqp::DAQP.Model)
  workspace = unsafe_load(daqp.work);
  if(workspace.settings != C_NULL)
	return unsafe_load(workspace.settings)
  end
end
function settings(daqp::DAQP.Model,changes::Dict{Symbol,<:Any})
  workspace = unsafe_load(daqp.work);
  if(workspace.settings == C_NULL) return end
  settings = unsafe_load(workspace.settings)
  new = [haskey(changes,f) ? changes[f] : getfield(settings,f) 
		 for f in fieldnames(DAQP.DAQPSettings)];
  new_settings = DAQP.DAQPSettings(new...)
  unsafe_store!(workspace.settings,new_settings);
  return new_settings;
end


end
