function [u,M,dupper,dlower,sense] = generate_test_LDP(n,m,ms,nActive,cond)


    H = eye(n);
    % * Generate solution to min ||u||_2^2 s.t. dlower <= M u <= dupper *
    Mg = randn(m-ms,n);
    [U,S,V] = svd(Mg);
    S(S~=0) = linspace(cond,1,min(m-ms,n));

    Mg = U*S*V';

    M = [H(1:ms,:);Mg]; % First row ms rows corresponds to simple bounds 
    dupper = zeros(m,1);
    dlower= zeros(m,1);
    shuffle_inds= randperm(m);
    nActive_upper = randperm(nActive+1,1)-1;
    nActive_lower= nActive-nActive_upper; 
    ids_active_upper = shuffle_inds(1:nActive_upper);
    ids_active_lower = shuffle_inds(nActive_upper+1:nActive);
    ids_inactive = shuffle_inds(nActive+1:m);

    % * Construct active bounds such that lam>=0 *
    lam = rand(nActive,1); % make lam>=0 by construction
    Ma = [M(ids_active_upper,:);-M(ids_active_lower,:)];
    da=  -Ma*(Ma')*lam;
    dupper(ids_active_upper) = da(1:nActive_upper);
    dlower(ids_active_lower) = -da(nActive_upper+1:nActive);
    u = -Ma'*lam;

    % * Make the inactive constraints feasible *
    bounds_gap = 1; % Scaling factor for distance between bounds
    slack_gap= 1; % Scaling factor for distance between bounds and optimizer 
    dupper(ids_active_lower) = dlower(ids_active_lower)+bounds_gap*(0.01+rand(nActive_lower,1));
    dlower(ids_active_upper) = dupper(ids_active_upper)-bounds_gap*(0.01+rand(nActive_upper,1));

    dupper(ids_inactive) = M(ids_inactive,:)*u + slack_gap*(0.01+rand(length(ids_inactive),1)); 
    dlower(ids_inactive) = M(ids_inactive,:)*u - slack_gap*(0.01+rand(length(ids_inactive),1)); 


    M = M(ms+1:end,:);  % Remove simple bounds since they are defined implicitly 
    sense=zeros(m,1,'int32');
end
