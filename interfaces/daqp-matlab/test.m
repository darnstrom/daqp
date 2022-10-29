%%  Basic
H = [];
f = [];
A = [1  0; 0 1; 1 0; 0 1];
bupper = [1; inf; inf; 0];
blower = [-inf; -1; 0.5; -inf];
sense= int32([8;8;8;8]);
break_points = [2;4];

d = daqp();
d.setup(H,f,A,bupper,blower,sense,break_points);
[x,fval,exitflag, hier_info] = d.solve(); 
hier_info
x

%%  Larger example 
rng(3);
break_points = [5;13;20;25;30;40;43;50];
n = 25; 
m = break_points(end);
bu_cpy = zeros(m,1)
bl_cpy = zeros(m,1)

H = []; 
f = []; 
A = randn(m,n); 
scale = 10;
bupper = scale*randn(m,1);
blower = bupper-scale*rand(m,1);
sense= int32(8*ones(m,1));

d = daqp();
bu_cpy(:) = bupper;
bl_cpy(:) = blower;
d.setup(H,f,A,bu_cpy,bl_cpy,sense,break_points);
[x_hi,fval,exitflag, hier_info] = d.solve(); 
hier_info


%% Solve without hierarchy
d = daqp();
bu_cpy(:) = bupper;
bl_cpy(:) = blower;
H = eye(n);
f = zeros(n,1);
sense= int32(8*ones(m,1));
d.setup(H,f,A,bu_cpy,bl_cpy,sense);
d.settings('rho_soft',1e-10);
[x_ref,fval,exitflag, hier_ref_info] = d.solve(); 
hier_ref_info
%% Compute slacks
start = 1;
slacks_hier = zeros(length(break_points),1);
for i = 1:length(break_points)
    inds = start:break_points(i);
    slacks_hier(i) = min(min(bupper(inds)-A(inds,:)*x_hi),... 
        min(-(blower(inds)-A(inds,:)*x_hi)));
    start = break_points(i)+1;
end
start = 1;
ref_slacks_hier = zeros(length(break_points),1);
for i = 1:length(break_points)
    inds = start:break_points(i);
    ref_slacks_hier(i) = min(min(bupper(inds)-A(inds,:)*x_ref),... 
        min(-(blower(inds)-A(inds,:)*x_ref)));
    start = break_points(i)+1;
end
[slacks_hier,ref_slacks_hier]

