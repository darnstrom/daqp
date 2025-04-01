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

