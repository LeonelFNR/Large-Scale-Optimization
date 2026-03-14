% In this script we introduce several test cases
% and call the function Compare_Performance to 
% compare primal affine scaling vs. lineprog
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')

%% FEASIBLE 1 — Small
A = [2 1 -1 0;
     3 4  0 1];

b = [2;
     12];

c = [3;1;0;0];

Compare_Performance(A,b,c);

%% FEASIBLE 2 — Multiple active constraints
A = [1 1 1 0 0;
     2 1 0 1 0;
     1 3 0 0 1];

b = [6;
     7;
     8];

c = [3;2;0;0;0];

Compare_Performance(A,b,c);

%% FEASIBLE 3 — Poorly scaled problem
A = [1e-3  2    -1  0;
     5e2   1     0  1];

b = [1;
     1000];

c = [1;3;0;0];

Compare_Performance(A,b,c);


%% UNBOUNDED 1 — Simple ray
A = [1 -1 1];

b = [1];

c = [-1;0;0];   % minimizar -x1 ⇒ crecer x1 indefinidamente

Compare_Performance(A,b,c);

%% UNBOUNDED 2 — Higher dimension
A = [1 -1 0 0 0;
     0  0 1 1 0;
     0  0 0 0 1];

b = [0;
     2;
     1];

c = [-1;0;0;0;0];

Compare_Performance(A,b,c);

%% INFEASIBLE 1 — Contradictory constraints
A = eye(4);

b = [-1000;
     -1000;
     -1000;
     -1000];

c = [1;1;1;1];
Compare_Performance(A,b,c);


%% MULTIPLE OPTIMA
A = [1 1 1];

b = [2];

c = [1;1;0];  % x3 free ⇒ infinite optimal solutions

Compare_Performance(A,b,c);


%% RANDOM LARGE FEASIBLE
rng(2);
m = 10;
n = 20;

A = rand(m,n);
x_true = rand(n,1);
b = A*x_true;     % garantizamos factibilidad
c = rand(n,1);

Compare_Performance(A,b,c);


%% lp_25fv47.mat
load('lp_25fv47.mat');  % Load the linear programming problem data
Compare_Performance(Problem.A, Problem.b, Problem.aux.c, true);  % Compare performance on the loaded problem

%% lp_adlittle.mat
load('lp_adlittle.mat');  % Load the next linear programming problem data
Compare_Performance(Problem.A, Problem.b, Problem.aux.c);  % Compare performance on the loaded problem

%% lp_afiro.mat
load('lp_degen2.mat');  % Load the next linear programming problem data
Compare_Performance(Problem.A, Problem.b, Problem.aux.c);  % Compare performance on the loaded problem


%% 
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:singularMatrix')