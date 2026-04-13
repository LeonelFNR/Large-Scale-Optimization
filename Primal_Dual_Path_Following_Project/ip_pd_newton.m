function [x,lambda,s,history] = ip_pd_newton(A,b,c,info)

if nargin < 4
    info = 1;
end

[m,n] = size(A);

tol = 1e-8;
rho = 0.99;
max_iter = 100;

x = ones(n,1)*10;
s = ones(n,1)*10;
lambda = zeros(m,1);

e = ones(n,1);

history = [];

if info > 0
    fprintf('\n-------------------------------------------------------------\n');
    fprintf('Iter |   f_primal    f_dual    ||rc||    ||rb||       mu\n');
    fprintf('-------------------------------------------------------------\n');
end

for k = 1:max_iter
    
    rc = A'*lambda + s - c;
    rb = A*x - b;
    mu = (x'*s)/n;
    
    f_primal = c'*x;
    f_dual   = b'*lambda;
    
    history = [history; [f_primal, f_dual, norm(rc), norm(rb), mu]];
    
    if info > 0
        fprintf('%3d  | %10.4e %10.4e %10.2e %10.2e %10.2e\n',...
            k-1,f_primal,f_dual,norm(rc),norm(rb),mu);
    end
    
    if norm(rc)/(1+norm(c)) < tol && norm(rb)/(1+norm(b)) < tol && mu < tol
        break;
    end
    
    X = diag(x);
    S = diag(s);
    
    rxs = X*S*e;
    
    [dx,dl,ds] = solve_kkt_system(A,rc,rb,rxs,X,S);
    
    alpha_p = min([1; rho*min(-x(dx<0)./dx(dx<0))]);
    alpha_d = min([1; rho*min(-s(ds<0)./ds(ds<0))]);
    
    if isempty(alpha_p); alpha_p = 1; end
    if isempty(alpha_d); alpha_d = 1; end
    
    x = x + alpha_p*dx;
    lambda = lambda + alpha_d*dl;
    s = s + alpha_d*ds;
end

if info > 0
    fprintf('-------------------------------------------------------------\n');
end