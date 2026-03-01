function x = PrimalAffineScaling(A,b,c)
%PRIMALAFFINESCALING : Takes as input the constraint matrix A that has n
%columns
%the independent term vector b and
%the; cost vector c and applies the PrimalAffineScaling algorithm

M = 1000; % arbitrary neihter big neither small value for big M algorithm
eps_components = -1e-12; % tolerance for delta X components
n = size(A,2);

% Compute infeasibilities
en = ones(n,1);
r = b - A*en;

% Extend A with additional column r
A = [A r];
n = n+1;

%extend cost vector
c = [c' M]';

% create x0 as 1...1 (n+1 ones)
x = ones(n,1);

%set D = I, y = (ADA')^(-1)ADc
D = eye(n); % Initialize D as the identity matrix
y = (A * D * A') \ (A * D * c); % Compute y using the formula

% set k = 0, eps = 10^(-6) and rho in [0.95, 0.9995]
k = 0;
eps = 1e-6;
rho = 0.995;

while abs(c'*x-b'*y) / (1+abs(c'*x)) > eps

    %compute z = c - A'*y
    z = c-A'*y;

    %compute deltax = -D*z
    dx = -D*z;
    if all(dx > eps_components)
        % STOP and print that problem is unbounded
        fprintf('The problem is unbounded.\n');
        return;
    end

    %compute alfa
    negIdx = dx < 0;                 % índices donde dx_i < 0
    alpha = rho * min(-x(negIdx) ./ dx(negIdx));

    % compute x_k+1 = x_k + alpha*dx
    x = x + alpha*dx;

    % Update k
    k = k + 1;

    % Compute X^k: X^k = diag(x_1^k, ..., x_n^k) and
    %Compute D = (X^k)^2 in one step
    D = diag(x.^2);

    % y for the next iteration
    y = (A * D * A') \ (A * D * c);
    
end

% if last component of x is not ero then STOP INFEASIBLE
% else STOP and the optimal solution is the found x, return it
if abs(x(end)) > eps
    fprintf('The problem is infeasible.\n');
else
    % return x without the artifical variable
    x = x(1:end-1);
    fprintf('Optimal solution found.\n');
end