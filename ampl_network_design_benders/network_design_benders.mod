set N;
set O within N;
set A within N cross N;
set Ahat within N cross N;
set AA := A union Ahat;

param xc {N};
param yc {N};
param t {N,O} default 0;
param rho > 0;
param Niter integer >= 1 default 100;
param eps >= 0 default 1e-6;

param c {(i,j) in AA, l in O} := 100 + (xc[i]-xc[j])^2 + (2*yc[i]-yc[j])^2;
param f {(i,j) in Ahat} := 10*(abs(xc[i]-xc[j]) + abs(yc[i]-yc[j]));

param yb {(i,j) in Ahat} binary default 0;
param nCUT integer >= 0 default 0;
param u {N,O,1..nCUT} default 0;
param tauMinus {Ahat,O,1..nCUT} >= 0 default 0;
param cut_const {1..nCUT} default 0;

node I {i in N, l in O}: net_out t[i,l];
arc x {(i,j) in AA, l in O} >= 0, from I[i,l], to I[j,l];

var y {(i,j) in Ahat} binary;
var zmp >= 0;

# Problema original, para resolver sin Benders
minimize Total_Cost:
    sum {(i,j) in Ahat} f[i,j]*y[i,j]
  + sum {l in O, (i,j) in AA} c[i,j,l]*x[i,j,l];

subject to CapOrig {(i,j) in Ahat, l in O}:
    x[i,j,l] <= rho*y[i,j];

# Subproblema primal con y fijado en yb
minimize Sub_Exploit_Cost:
    sum {l in O, (i,j) in AA} c[i,j,l]*x[i,j,l];

subject to CapSub {(i,j) in Ahat, l in O}:
    x[i,j,l] <= rho*yb[i,j];

# Master de Benders: zmp representa coste total inversión + explotación aproximado
minimize Master_Obj:
    zmp;

subject to Bcut {k in 1..nCUT}:
    zmp >= sum {(i,j) in Ahat} f[i,j]*y[i,j]
         + cut_const[k]
         - rho * sum {l in O, (i,j) in Ahat} tauMinus[i,j,l,k]*y[i,j];
