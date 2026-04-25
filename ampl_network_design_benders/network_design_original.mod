set N;
set O within N;
set A within N cross N; #arcos fijos
set Ahat within N cross N; #arcos a�adibles
set AA := A union Ahat; #todos los arcos

param xc {N};
param yc {N};
param t {N,O} default 0;
param rho > 0;

# Coste de explotación. Si tu enunciado usa 100*(dx)^2 en vez de 100+(dx)^2,
param c {(i,j) in AA, l in O} := 100 + (xc[i]-xc[j])^2 + (2*yc[i]-yc[j])^2;
param f {(i,j) in Ahat} := 10*(abs(xc[i]-xc[j]) + abs(yc[i]-yc[j]));

node I {i in N, l in O}: net_out t[i,l];
arc x {(i,j) in AA, l in O} >= 0, from I[i,l], to I[j,l];

var y {(i,j) in Ahat} binary;

minimize Total_Cost:
    sum {(i,j) in Ahat} f[i,j]*y[i,j]
  + sum {l in O, (i,j) in AA} c[i,j,l]*x[i,j,l];

subject to CapCand {(i,j) in Ahat, l in O}:
    x[i,j,l] <= rho*y[i,j];

minimize Investment_Cost:
    sum {(i,j) in Ahat} f[i,j]*y[i,j];

minimize Exploitation_Cost:
    sum {l in O, (i,j) in AA} c[i,j,l]*x[i,j,l];
