function x = PrimalAffineScaling(A,b,c, scaling)
%PRIMALAFFINESCALING : Takes as input the constraint matrix A that has n
%columns
%the independent term vector b and
%the; cost vector c and applies the PrimalAffineScaling algorithm
%scaling is a boolean. it's default value is true. if false it means no
%scaling

if nargin < 4
    scaling = true;
end

% set k = 0, eps = 10^(-6) and rho in [0.95, 0.9995]
k = 0;
eps = 1e-8;
rho = 0.95;
obj_unbounded_tol = 1e24;
% tol_result = 1e-11;

%work on A so that we only work with full row rank matrix
[Q,R,E] = qr(A',0);          % QR sobre A'
tol = 1e-12;
r = sum(abs(diag(R)) > tol);
ind = sort(E(1:r));          % filas independientes de A
A = A(ind,:);
b = b(ind);


M = 1e6 * max(1, norm(c, inf)); % arbitrary neihter big neither small value for big M algorithm
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

% lhs = A * D * A';
% rhs = A * D * c;
% y = lhs \ rhs;
y = (A * D * A') \ (A * D * c); % Compute y using the formula


while abs(c'*x-b'*y) / (1+abs(c'*x)) > eps

    %compute z = c - A'*y
    z = c-A'*y;

    %compute deltax = -D*z
    dx = -D*z;
    if all(dx >= 0)
        display(dx);
        % STOP and print that problem is unbounded
        fprintf('The problem is unbounded.\n');
        return;
    end

    %compute alfa
    negIdx = dx < 0;                 % indices where dx_i < 0
    alpha = rho * min(-x(negIdx) ./ dx(negIdx));

    % compute x_k+1 = x_k + alpha*dx
    x = x + alpha*dx;

    % Update k
    k = k + 1;

    if scaling == true
        % Compute X^k: X^k = diag(x_1^k, ..., x_n^k) and
        %Compute D = (X^k)^2 in one step
        D = diag(x.^2);
    end

    y = A * D * A' \ A * D * c;
    
end

% for printing information later
obj_value = c'*x;

%if the obj_value is too high (on absolute value), deem the problem as
%unbounded
if abs(obj_value) > obj_unbounded_tol
    fprintf('Status: LIKELY UNBOUNDED\n');
    fprintf('Heuristic triggered: |c''x| too large.\n');
    x = NaN(size(x));
    return;
end


fprintf('\n===== Primal Affine Scaling Results =====\n');
fprintf('Iterations: %d\n', k);
fprintf('Objective value (extended problem): %.10f\n', obj_value);

% if last component of x is not zero then STOP INFEASIBLE
% else STOP and the optimal solution is the found x, return it
%display(x);
if abs(x(end)) > 1e-12
    fprintf('Status: INFEASIBLE\n');
    fprintf('%d', x(end));
else
    x = x(1:end-1);  % remove artificial variable
    fprintf('Status: OPTIMAL\n');
    fprintf('Optimal objective value: %.10f\n', c(1:end-1)'*x);
end