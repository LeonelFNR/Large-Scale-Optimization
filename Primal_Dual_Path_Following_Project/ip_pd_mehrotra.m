function [x,lambda,s,history] = ip_pd_mehrotra(A,b,c,info)

if nargin < 4
    info = 1; % show logs by default
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
    
    % residues
    rc = A'*lambda + s - c;
    rb = A*x - b;
    mu = (x'*s)/n;
    
    f_primal = c'*x;
    f_dual   = b'*lambda;
    
    % save history
    history = [history; [f_primal, f_dual, norm(rc), norm(rb), mu]];
    
    % print
    if info > 0
        fprintf('%3d  | %10.4e %10.4e %10.2e %10.2e %10.2e\n',...
            k-1,f_primal,f_dual,norm(rc),norm(rb),mu);
    end
    
    % stopping criteria
    if norm(rc)/(1+norm(c)) < tol && norm(rb)/(1+norm(b)) < tol && mu < tol
        break;
    end
    
    X = diag(x);
    S = diag(s);
    
    %% --------- PREDICTOR ---------
    rxs_aff = X*S*e;
    
    [dx_aff,dl_aff,ds_aff] = solve_kkt_system(A,rc,rb,rxs_aff,X,S);
    
    alpha_aff_p = min([1; min(-x(dx_aff<0)./dx_aff(dx_aff<0))]);
    alpha_aff_d = min([1; min(-s(ds_aff<0)./ds_aff(ds_aff<0))]);
    
    if isempty(alpha_aff_p); alpha_aff_p = 1; end
    if isempty(alpha_aff_d); alpha_aff_d = 1; end
    
    x_aff = x + alpha_aff_p*dx_aff;
    s_aff = s + alpha_aff_d*ds_aff;
    
    mu_aff = (x_aff'*s_aff)/n;
    
    sigma = (mu_aff/mu)^3;
    
    %% --------- CORRECTOR ---------
    rxs = X*S*e + dx_aff.*ds_aff - sigma*mu*e;
    
    [dx,dl,ds] = solve_kkt_system(A,rc,rb,rxs,X,S);
    
    %% --------- STEP LENGTH ---------
    alpha_p = min([1; rho*min(-x(dx<0)./dx(dx<0))]);
    alpha_d = min([1; rho*min(-s(ds<0)./ds(ds<0))]);
    
    if isempty(alpha_p); alpha_p = 1; end
    if isempty(alpha_d); alpha_d = 1; end
    
    %% --------- UPDATE ---------
    x = x + alpha_p*dx;
    lambda = lambda + alpha_d*dl;
    s = s + alpha_d*ds;
end

if info > 0
    fprintf('-------------------------------------------------------------\n');
end