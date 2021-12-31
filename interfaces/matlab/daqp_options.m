function opts = daqp_options()
 opts.primal_tol = 1e-6; 
 opts.dual_tol = 1e-12; 
 opts.zero_tol = 1e-14; 
 opts.pivot_tol= 1e-2; 
 opts.progress_tol= 1e-6; 
 opts.cycle_tol = 10; 
 opts.iter_limit = 1000;
 opts.eps_prox = 0;
 opts.eta_prox = 1e-6;
 opts.prox_iter_limit=1000;
end
