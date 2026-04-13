function results = compare_with_linprog(A,b,c)

fprintf('\n================= COMPARISON =================\n');

n = length(c);

%% --------- 1. MEHROTRA ---------
tic;
[x_ipm,lambda_ipm,s_ipm,h_ipm] = ip_pd_mehrotra(A,b,c,0);
t_ipm = toc;

f_ipm = c'*x_ipm;
iter_ipm = size(h_ipm,1);

%% --------- 2. NEWTON ---------
tic;
[x_new,lambda_new,s_new,h_new] = ip_pd_newton(A,b,c,0);
t_new = toc;

f_new = c'*x_new;
iter_new = size(h_new,1);

%% --------- 3. LINPROG ---------
options = optimoptions('linprog','Display','none');

tic;
[x_lp,f_lp,exitflag,output] = linprog(c,[],[],A,b,zeros(n,1),[],options);
t_lp = toc;

iter_lp = output.iterations;

%% --------- PRINT RESULTS ---------

fprintf('\nMethod        f(x)         Iter     Time(s)\n');
fprintf('-------------------------------------------------\n');
fprintf('Mehrotra   %10.6f    %4d     %8.4f\n',f_ipm,iter_ipm,t_ipm);
fprintf('Newton     %10.6f    %4d     %8.4f\n',f_new,iter_new,t_new);
fprintf('linprog    %10.6f    %4d     %8.4f\n',f_lp,iter_lp,t_lp);

fprintf('-------------------------------------------------\n');

%% --------- ERRORS ---------

fprintf('\nErrors vs linprog:\n');

err_ipm = norm(x_ipm - x_lp);
err_new = norm(x_new - x_lp);

fprintf('||x_ipm - x_lp|| = %.2e\n',err_ipm);
fprintf('||x_new - x_lp|| = %.2e\n',err_new);

%% --------- SAVE RESULTS ---------

results.ipm.x = x_ipm;
results.ipm.f = f_ipm;
results.ipm.iter = iter_ipm;
results.ipm.time = t_ipm;

results.newton.x = x_new;
results.newton.f = f_new;
results.newton.iter = iter_new;
results.newton.time = t_new;

results.linprog.x = x_lp;
results.linprog.f = f_lp;
results.linprog.iter = iter_lp;
results.linprog.time = t_lp;

results.errors.ipm = err_ipm;
results.errors.newton = err_new;

fprintf('===============================================\n');