function [dx,dl,ds] = solve_kkt_system(A,rc,rb,rxs,X,S)

n = length(rc);
m = length(rb);

% Theta = X * S^{-1}
Theta = X / S;

% normal matrix equations
M = A * Theta * A';

% RHS
rhs = -rb + A*(Theta*(-rc) + S\rxs);

% solve system
dl = M \ rhs;

% recuperar ds y dx
ds = -rc - A'*dl;
dx = -(S \ (rxs + X*ds));