using Random
## QP
function generate_test_QP(n,m,ms,nActive,kappa)
  # * Transform QP to LDP with T = L Q' (where L is diagonal and Q is orthogonal)
  # to get solution  u =T(x+f) (=> H = T'*T, f = T'*v)
  # Construct eigenvalues such that cond(H) = kappa
  eigens = zeros(n);
  eigens[1] = 1; eigens[2] = kappa;
  eigens[3:end] = 1 .+ (kappa-1)*rand(n-2);
  # Randomly generate an orthogonal matrix
  Q = qr(randn(n,n)).Q;
  T= diagm(sqrt.(eigens))*Q';
  Tinv = Q*diagm( 1 ./ sqrt.(eigens));
  H = T' * T;

  # * Generate solution to min ||u||_2^2 s.t. dlower <= M u <= dupper *
  M = [Tinv[1:ms,:]; randn(m-ms,n)]; # First row ms rows corresponds to simple bounds
  dupper = zeros(m);
  dlower= zeros(m);
  shuffle_inds= shuffle(1:m);
  nActive_upper = rand(0:nActive);
  nActive_lower= nActive-nActive_upper;
  ids_active_upper = shuffle_inds[1:nActive_upper];
  ids_active_lower = shuffle_inds[nActive_upper+1:nActive];
  ids_inactive = shuffle_inds[nActive+1:m];

  # * Construct active bounds such that lam>=0 *
  lam = rand(nActive); # make lam>=0 by construction
  Ma = [M[ids_active_upper,:];-M[ids_active_lower,:]];
  da=  -Ma*(Ma')*lam;
  dupper[ids_active_upper] = da[1:nActive_upper];
  dlower[ids_active_lower] = -da[nActive_upper+1:nActive];
  u = -Ma'*lam;

  # * Make the inactive constraints feasible *
  bounds_gap = 1; # Scaling factor for distance between bounds
  slack_gap= 1; # Scaling factor for distance between bounds and optimizer
  dupper[ids_active_lower] = dlower[ids_active_lower] .+ bounds_gap*(0.01 .+ rand(nActive_lower));
  dlower[ids_active_upper] = dupper[ids_active_upper] .- bounds_gap*(0.01 .+ rand(nActive_upper));

  dupper[ids_inactive] = M[ids_inactive,:]*u + slack_gap*(0.01 .+ rand(length(ids_inactive)));
  dlower[ids_inactive] = M[ids_inactive,:]*u - slack_gap*(0.01 .+ rand(length(ids_inactive)));


  # Compute solution x = T\(u-v) (which yields f = (T)'*v)
  v= randn(n); f = T'*v;
  x = (T)\(u-v);
  # Transform constraint to x-domain (=> A = M*R, b = d-M*v)
  A = M[ms+1:end,:]*T;  # Remove simple bounds since they are defined implicitly
  bupper = dupper-M*v; blower = dlower-M*v;
  sense=zeros(Cint,m)
  return x,H,f,A,bupper,blower,sense
end

## LP
function generate_test_LP(n,m,ms)
    A = [Matrix(I(n)[1:ms,:]);randn(m-ms,n)];
    bupper = zeros(m);
    blower = zeros(m);
    shuffle_inds= randperm(m);
    nActive_upper = rand(1:n+1) - 1;
    nActive_lower= n - nActive_upper;
    ids_active_upper = shuffle_inds[1:nActive_upper];
    ids_active_lower = shuffle_inds[nActive_upper+1:n];
    ids_inactive = shuffle_inds[n+1:m];

    λ = rand(n); # λ ≥ 0 (dual feasible)
    x = randn(n);

    Aa = [A[ids_active_upper,:];-A[ids_active_lower,:]];
    f = -Aa'*λ;

    ba = Aa*x;
    bupper[ids_active_upper] = ba[1:nActive_upper];
    blower[ids_active_lower] = -ba[nActive_upper+1:n];


    # * Make the inactive constraints feasible *
    bounds_gap = 1; # Scaling factor for distance between bounds
    slack_gap= 1; # Scaling factor for distance between bounds and optimizer
    bupper[ids_active_lower] = blower[ids_active_lower]+bounds_gap*(0.01 .+ rand(nActive_lower));
    blower[ids_active_upper] = bupper[ids_active_upper]-bounds_gap*(0.01 .+ rand(nActive_upper));

    bupper[ids_inactive] = A[ids_inactive,:]*x + slack_gap*(0.01 .+ rand(length(ids_inactive)));
    blower[ids_inactive] = A[ids_inactive,:]*x - slack_gap*(0.01 .+ rand(length(ids_inactive)));
    A = A[ms+1:end,:];  #simple bounds are implicitly defined, so remove them from A.
    sense = zeros(Int32,m)
    return x,f,A,bupper,blower,sense
end
## Equality constrained QP
"""
    generate_test_QP_eq(n,m,ms,nActive,neq,kappa)

Generate a QP with `neq` equality constraints and a known solution.

`neq` of the constraints that are active at the optimum of the QP from
`generate_test_QP` are turned into equalities by collapsing their bounds onto
the active value. The multipliers of these constraints are nonnegative at the
optimum, so removing the sign requirement on them leaves the KKT conditions
satisfied and `xref` remains the solution.

General constraints are preferred over simple bounds, so that the equality
elimination in the solver operates on the rows of `A`.
"""
function generate_test_QP_eq(n,m,ms,nActive,neq,kappa)
  neq <= nActive || error("generate_test_QP_eq requires neq <= nActive")
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nActive,kappa)

  # The first ms constraints are the simple bounds x[1:ms]; the rest are A*x
  r = [xref[1:ms]; A*xref]
  # Active constraints sit exactly on a bound, inactive ones are >= 0.01 away
  is_active(i) = min(abs(r[i]-bupper[i]), abs(r[i]-blower[i])) <= 1e-6*(1+abs(r[i]))
  # Prefer general constraints (i > ms) over simple bounds
  active = [i for i in ms+1:m if is_active(i)]
  append!(active, [i for i in 1:ms if is_active(i)])

  length(active) >= neq ||
    error("only found $(length(active)) active constraints, need $neq")

  for i in active[1:neq]
    # Collapse the bounds onto whichever one is active at xref
    b = abs(r[i]-bupper[i]) <= abs(r[i]-blower[i]) ? bupper[i] : blower[i]
    bupper[i] = b
    blower[i] = b
    sense[i] = Cint(DAQPBase.EQUALITY)
  end

  return xref,H,f,A,bupper,blower,sense
end

## MIQP
"""
    generate_test_MIQP(n,m,ms,nb)

Generate a mixed-integer QP whose first `nb` variables are binary (requires
`ms >= nb`, since the binary constraints are placed on simple bounds).

The origin is feasible, and `f` is skewed towards the binary variables so that
it is lucrative to leave it. A cardinality constraint caps how many binaries
can be set, which makes the relaxation fractional and forces branch and bound
to actually search instead of diving straight to an integer solution.

Requires `m-ms >= 1` for the cardinality constraint.
"""
function generate_test_MIQP(n,m,ms,nb)
  ms >= nb || error("generate_test_MIQP requires ms >= nb")
  m-ms >= 1 || error("generate_test_MIQP requires at least one general constraint")
  M = randn(n,n);
  H = M'*M + I;                                # PD and reasonably conditioned
  A = randn(m-ms,n);
  bupper = 20*rand(m); blower = -20*rand(m);   # ensure that the origin is feasible
  f = 100*randn(n); f[1:nb] .= -abs.(f[1:nb]); # make it lucrative to leave the origin
  bupper[1:nb] .= 1.0;
  blower[1:nb] .= 0.0;
  sense = zeros(Cint,m);
  sense[1:nb] .= Cint(DAQPBase.BINARY);

  # Cardinality constraint sum(x[1:nb]) <= nb/2: the binaries compete, so the
  # relaxation is fractional and the search has to branch and backtrack.
  A[1,:] .= 0.0;
  A[1,1:nb] .= 1.0;
  bupper[ms+1] = floor(nb/2);
  blower[ms+1] = -1e30;

  return H,f,A,bupper,blower,sense
end

## AVI
function generate_test_avi(n,m)
    A = randn(m,n);
    shuffle_inds= randperm(m);
    nAS= rand(1:n+1) - 1;
    AS = shuffle_inds[1:nAS];

    λ = zeros(m);λ[AS]=rand(nAS);
    x = randn(n)

    # Generate H that is PD but not symmetric
    M,N = rand(n,n),randn(n,n)
    sym = M'*M
    asym = N-N'
    H = sym/norm(sym) + asym/norm(asym)

    # f to satisfy stationarity
    f = -H*x-A[AS,:]'*λ[AS]

    # b to make AS active and IS inactive
    Ax = A*x
    b = Ax + 5*rand(m)
    b[AS] = Ax[AS]

    return x,H,f,A,b
end
## Source files
function get_local_sources(srcdir)
    # Get local source
    # <repo>/interfaces/daqp-julia/test -> <repo>
    daqp_dir = joinpath(dirname(@__FILE__), "..","..","..","..")
    cfiles = ["daqp.c","auxiliary.c","factorization.c", "bnb.c", "hierarchical.c"]
    hfiles = ["daqp.h","auxiliary.h","factorization.h", "bnb.h", "hierarchical.h","constants.h", "types.h"]
    for cf in cfiles
        cp(joinpath(daqp_dir,"src",cf), joinpath(srcdir,cf))
    end
    for hf in hfiles
        cp(joinpath(daqp_dir,"include",hf), joinpath(srcdir,hf))
    end
end
