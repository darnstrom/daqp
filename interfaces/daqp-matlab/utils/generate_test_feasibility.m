function [A,bupper,blower,sense] = generate_test_feasibility(n,m,r)
    c = randn(n,1);
    A = randn(m,n);
    A = A./vecnorm(A,2,2);
    bupper = A*c+r;
    blower = -1e30*ones(m,1);
    sense=zeros(m,1,'int32');
end
