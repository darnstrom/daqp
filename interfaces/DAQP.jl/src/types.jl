mutable struct QPj
  n::Cint
  m::Cint
  ms::Cint
  H::Matrix{Cdouble}
  f::Vector{Cdouble}
  A::Matrix{Cdouble}
  bupper::Vector{Cdouble}
  blower::Vector{Cdouble}
  sense::Vector{Cint}
end
function QPj() 
  return QPj(0,0,0,Matrix{Cdouble}(undef,0,0),Vector{Cdouble}(undef,0),Matrix{Cdouble}(undef,0,0), Vector{Cdouble}(undef,0), Vector{Cdouble}(undef,0), Vector{Cint}(undef,0)) 
end
function QPj(H::Matrix{Float64},f::Vector{Float64},
	A::Matrix{Float64},bupper::Vector{Float64}, blower::Vector{Float64},
	sense::Vector{Cint})
  # TODO: check consistency of dimensions
  (mA,n) = size(A);
  m = length(bupper);
  ms = m-mA;
  return QPj(n,m,ms,H,f,A',bupper,blower,sense) # Transpose A for col => row major
end

struct QPc 
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
function QPc(qpj::QPj)
  return QPc(qpj.n,qpj.m,qpj.ms,
			 pointer(qpj.H),pointer(qpj.f),
			 pointer(qpj.A),pointer(qpj.bupper),pointer(qpj.blower),pointer(qpj.sense))
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

struct Workspace
  qp::Ptr{QPc}
  n::Cint
  m::Cint
  ms::Cint
  M::Ptr{Cdouble}
  dupper::Ptr{Cdouble}
  dlower::Ptr{Cdouble}
  Rinv::Ptr{Cdouble}
  v::Ptr{Cdouble}
  sense::Ptr{Cint}

  x::Ptr{Cdouble}
  xold::Ptr{Cdouble}

  lam::Ptr{Cdouble}
  lam_star::Ptr{Cdouble}

  u::Ptr{Cdouble}
  fval::Cdouble
  fval_bound::Cdouble


  L::Ptr{Cdouble}
  D::Ptr{Cdouble}
  xldl::Ptr{Cdouble}
  zldl::Ptr{Cdouble}
  reuse_ind::Cint

  WS::Ptr{Cint}
  n_active::Cint

  iterations::Cint
  outer_iter::Cint
  inner_iter::Cint
  sing_ind::Cint

  soft_slack::Cdouble

  settings::Ptr{DAQPSettings}
end

