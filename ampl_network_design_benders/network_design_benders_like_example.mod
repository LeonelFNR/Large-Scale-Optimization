set N;
set O within N;
set A within N cross N;       # fixed/existing arcs
set Ahat within N cross N;    # candidate arcs
set AA := A union Ahat;       # all arcs

param xc {N};
param yc {N};
param t {i in N, l in O};

param c {(i,j) in AA, l in O} :=
    100 + (xc[i]-xc[j])^2 + (2*yc[i]-yc[j])^2;

param f {(i,j) in Ahat} :=
    10*(abs(xc[i]-xc[j]) + abs(yc[i]-yc[j]));

param rho > 0;
param Niter;
param eps >= 0 default 1e-6;

# Current fixed value of the master binary variables, used by the subproblem
param yb {(i,j) in Ahat} binary default 0;

# Benders information accumulated iteration by iteration
param nCUT integer >= 0 default 0;
param u {i in N, l in O, k in 1..nCUT} default 0;
param restric {(i,j) in Ahat, l in O, k in 1..nCUT} >= 0 default 0;
param ybk {(i,j) in Ahat, k in 1..nCUT} binary default 0;

# Network formulation of flow balance:
# net_out = t[i,l], so t>0 at origin and t<0 at destinations.
node I {i in N, l in O}: net_out = t[i,l];
arc xl {(i,j) in AA, l in O} >= 0, from I[i,l], to I[j,l];

# Binary design variables
var y {(i,j) in Ahat} binary;

# ------------------------------------------------------------
# Original problem, without Benders
# ------------------------------------------------------------
minimize z:
    sum {(i,j) in Ahat} f[i,j]*y[i,j]
  + sum {l in O, (i,j) in AA} c[i,j,l]*xl[i,j,l];

subject to caps1 {(i,j) in Ahat, l in O}:
    xl[i,j,l] <= rho*y[i,j];

# ------------------------------------------------------------
# Subproblem, with y fixed to yb
# ------------------------------------------------------------
minimize zd:
    sum {l in O, (i,j) in AA} c[i,j,l]*xl[i,j,l];

subject to caps {(i,j) in Ahat, l in O}:
    xl[i,j,l] <= rho*yb[i,j];

# ------------------------------------------------------------
# Master problem
# zmp3 represents total cost: investment + approximated exploitation
# ------------------------------------------------------------
var zmp3 >= 0;

minimize ZMP3:
    zmp3;

subject to Bcut {k in 1..nCUT}:
    zmp3 >=
        sum {(i,j) in Ahat} f[i,j]*y[i,j]
      + sum {l in O, i in N} t[i,l]*u[i,l,k]
      - rho * sum {l in O, (i,j) in Ahat} restric[i,j,l,k]*y[i,j];
