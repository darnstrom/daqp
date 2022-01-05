struct QP 
  n::Cint
  m::Cint
  ms::Cint
  H::Ptr{Cdouble}
  f::Ptr{Cdouble}
  A::Ptr{Cdouble}
  bupper::Ptr{Cdouble}
  blower::Ptr{Cdouble}
  sense::Ptr{Cint}
end
function QP(n::Int64,m::Int64,ms::Int64,
	H::Matrix{Float64},f::Vector{Float64},
	A::Matrix{Float64},bupper::Vector{Float64}, blower::Vector{Float64},
	sense::Vector{Int64})
  return QP(n,m,ms,pointer(H),pointer(f), 
			pointer(A),pointer(bupper),pointer(blower),pointer(sense))
end

struct DAQPSettings
  primal_tol::Cdouble
  dual_tol::Cdouble
  zero_tol::Cdouble
  pivot_tol::Cdouble
  progress_tol::Cdouble

  cycle_tol::Cint
  iter_limit::Cint

  eps_prox::Cdouble
  eta_prox::Cdouble
  prox_iter_limit::Cint

  rho_soft::Cdouble
end
function DAQPSettings()
  settings = Ref{DAQP.DAQPSettings}()
  ccall((:daqp_default_settings, DAQP.lib), Nothing,(Ref{DAQP.DAQPSettings},), settings)
  return settings[]
end

struct DAQPResult
  x::Ptr{Cdouble}
  fval::Cdouble
  soft_slack::Cdouble

  exitflag::Cint
  iter::Cint
  outer_iter::Cint
  solve_time::Cdouble
  setup_time::Cdouble
end

function DAQPResult(x::Vector{Float64})
  return DAQPResult(pointer(x),0,0,0,0,0,0,0)
end
