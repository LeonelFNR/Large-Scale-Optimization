function [dx,dl,ds] = solve_kkt_system(A,rc,rb,rxs,X,S)

n = length(rc);
m = length(rb);

% Theta = X * S^{-1}
Theta = X / S;

% matriz normal equations
M = A * Theta * A';

% RHS
rhs = -rb + A*(Theta*(-rc) + S\rxs);

% resolver sistema
dl = M \ rhs;

% recuperar ds y dx
ds = -rc - A'*dl;
dx = -(S \ (rxs + X*ds));