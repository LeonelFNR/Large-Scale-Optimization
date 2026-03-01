% Practical application of primal affine scaling algorithm

% The input parameters to the code will be:
%  the cost vector c, 
% the constraints matrix A, 
% and the right-hand-side vector b.  
% The constraints matrix A will have to be expanded 
% with a new variable needed by the Big-M method.

%Obtain somehow inputs A, b, c
% ===== TEST LP =====
% min 3x1 + x2
% s.t.
% 2x1 + x2 - s1 = 2
% 3x1 + 4x2 + s2 = 12
% x >= 0

A = [2 1 -1 0;
     3 4  0 1];

b = [2;
     12];

c = [3;
     1;
     0;
     0];   % slacks tienen coste 0
%Call for the method of primal affine scaling and return the x coords
x = PrimalAffineScaling(A,b,c);

%print x
disp('Optimal solution coordinates:');
disp(x);
