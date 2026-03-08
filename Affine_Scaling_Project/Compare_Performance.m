function Compare_Performance(A,b,c)
%COMPARE_PERFORMANCE Summary of this function goes here
% Takes as input the constraint matrix A that has n
%columns
%the independent term vector b and
%the; cost vector c and applies the PrimalAffineScaling algorithm
% Compares solutions and results between primal affine scaling and lineprog

fprintf('\n===== Problem characteristics =====\n');
fprintf('Number of variables: %d\n', size(A,2));
fprintf('Number of constraints: %d\n', size(A,1));
fprintf('Number of nonzeros in A: %d\n', nnz(A));

%Call for the method of primal affine scaling and return the x coords
x = PrimalAffineScaling(A,b,c);
fprintf('Coordinates of x: \n');

%disp(x);

%comparison linprog
options = optimoptions('linprog','Display','none');

[x_lp,fval_lp,exitflag,output] = linprog(c, ...
                                         [],[], ...
                                         A, b, ...
                                         zeros(size(A,2),1), [], ...
                                         options);

fprintf('\n===== linprog Results =====\n');
fprintf('Exit flag: %d\n', exitflag);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Objective value: %.10f\n', fval_lp);

% %print x
% disp('Optimal solution coordinates:');
% disp(x);
