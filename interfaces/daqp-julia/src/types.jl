struct QPj
    n::Cint
    m::Cint
    ms::Cint
    H::Matrix{Cdouble}
    f::Vector{Cdouble}
    A::Matrix{Cdouble}
    bupper::Vector{Cdouble}
    blower::Vector{Cdouble}
    sense::Vector{Cint}

    break_points::Vector{Cint}
    nh::Cint
end
function QPj() 
    return QPj(0,0,0,Matrix{Cdouble}(undef,0,0),Vector{Cdouble}(undef,0),Matrix{Cdouble}(undef,0,0), Vector{Cdouble}(undef,0), Vector{Cdouble}(undef,0), Vector{Cint}(undef,0),Vector{Cint}(undef,0),0)
end
function QPj(H::Matrix{Float64},f::Vector{Float64},
        A::Matrix{Float64},bupper::Vector{Float64}, blower::Vector{Float64},
        sense::Vector{Cint};A_rowmaj=false,break_points=Cint[])
    # TODO: check consistency of dimensions
    if(A_rowmaj)
        (n,mA) = size(A);
    else
        (mA,n) = size(A);
    end
    blower = isempty(blower) ? fill(-1e30,length(bupper)) : blower
    sense = isempty(sense) ? zeros(Cint,length(bupper)) : sense
    m = length(bupper);
    ms = m-mA;
    if(!A_rowmaj)
        A = A' # Transpose A for col => row major
    end
    return QPj(n,m,ms,H,f,A,bupper,blower,sense,break_points,length(break_points))
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

    break_points::Ptr{Cint}
    nh::Cint
end
function QPc(qpj::QPj)
    H_ptr = isempty(qpj.H) ? C_NULL : pointer(qpj.H)
    f_ptr = isempty(qpj.f) ? C_NULL : pointer(qpj.f)
    return QPc(qpj.n,qpj.m,qpj.ms,
               H_ptr,f_ptr,
               pointer(qpj.A),pointer(qpj.bupper),pointer(qpj.blower),pointer(qpj.sense),
               pointer(qpj.break_points),qpj.nh)
end

struct DAQPSettings
    primal_tol::Cdouble
    dual_tol::Cdouble
    zero_tol::Cdouble
    pivot_tol::Cdouble
    progress_tol::Cdouble

    cycle_tol::Cint
    iter_limit::Cint
    fval_bound::Cdouble

    eps_prox::Cdouble
    eta_prox::Cdouble

    rho_soft::Cdouble

    rel_subopt::Cdouble
    abs_subopt::Cdouble

    sing_tol::Cdouble
    refactor_tol::Cdouble
end
function DAQPSettings()
    settings = Ref{DAQPBase.DAQPSettings}()
    ccall((:daqp_default_settings, libdaqp), Nothing,(Ref{DAQPBase.DAQPSettings},), settings)
    return settings[]
end

struct DAQPResult
    x::Ptr{Cdouble}
    lam::Ptr{Cdouble}
    fval::Cdouble
    soft_slack::Cdouble

    exitflag::Cint
    iter::Cint
    nodes::Cint
    solve_time::Cdouble
    setup_time::Cdouble
end

function DAQPResult(x::Vector{Float64},lam::Vector{Float64})
    return DAQPResult(pointer(x),pointer(lam),0,0,0,0,0,0,0)
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
    scaling::Ptr{Cdouble}
    RinvD::Ptr{Cdouble}

    x::Ptr{Cdouble}
    xold::Ptr{Cdouble}

    lam::Ptr{Cdouble}
    lam_star::Ptr{Cdouble}

    u::Ptr{Cdouble}
    fval::Cdouble

    L::Ptr{Cdouble}
    D::Ptr{Cdouble}
    xldl::Ptr{Cdouble}
    zldl::Ptr{Cdouble}
    reuse_ind::Cint

    WS::Ptr{Cint}
    n_active::Cint

    iterations::Cint
    sing_ind::Cint

    soft_slack::Cdouble

    settings::Ptr{DAQPSettings}

    bnb::Ptr{Cvoid}

    nh::Cint
    break_points::Ptr{Cint}
end
