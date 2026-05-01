# Model modificat
set N;
set artN;    #conjunt unitari
set O within N;
set A within N cross N; #arcos fijos
#set Aart1{O}; 
param M;
param y{A} default M; 

param xc {N};
param yc {N};
param t {i in N,l in O};
param c {(i,j) in A} :=100 + (xc[i]-xc[j])^2 + (2*yc[i]-yc[j])^2;
param rho>0;
param nCUT;
param mu{(i,j) in A} >=0, default 0; # els multiplicadors
var xx{A};
node I {i in N, l in O}: net_out=t[i,l]; #si positivo ==> inyecci�n, si negativo extracci�n
node art{o in artN, l in O}: net_out=0;

arc xl {(i,j) in A, l in O}>=0: from I [i,l], to I [j,l];
arc art1{o in artN, l in O}>=0: from I [l,l], to art [o,l];
arc art2{o in artN, i in N, l in O} >=0: from art[o,l], to I[i,l];

subject to total_flow {(i,j) in A}: xx[i,j] = sum{l in O} xl[i,j,l];

minimize w: (sum {(i,j) in A} c[i,j]*xx[i,j])+
                sum{(i,j) in A} mu[i,j]*(xx[i,j] - y[i,j])+
                M*(sum{o in artN, l in O} art1[o,l])+
                M*(sum{o in artN, i in N, l in O} art2[o,i,l]);
subject to caps {(i,j) in A}: xx[i,j]<=y[i,j];

var z;
var mu0{A}>=0;
param YY{1..nCUT} binary default 1;
maximize Z: z;
param xxX{A,{1..nCUT}} default 0;
param art1X{o in artN,l in O, k in {1..nCUT}} default 0;
param art2X{o in artN, i in N, l in O, k in {1..nCUT}} default 0;
subject to cuts{k in {1..nCUT}}: z <= (sum {(i,j) in A} c[i,j]*xxX[i,j,k])+
                sum{(i,j) in A} mu0[i,j]*(xxX[i,j,k] - y[i,j])+
                M*(sum{o in artN, l in O} art1X[o,l,k])+
                M*(sum{o in artN, i in N, l in O} art2X[o,i,l,k])+
                M * YY[k];

# ------------------------------------------------------------
# GENERALIZED LINEAR PROGRAMMING (PART 7)
# ------------------------------------------------------------

var alpha {k in 1..nCUT} >= 0;

minimize GLP_Obj:
    sum {k in 1..nCUT} alpha[k] *
        sum {(i,j) in A} c[i,j] * xxX[i,j,k];

subject to Convexity:
    sum {k in 1..nCUT} alpha[k] = 1;

subject to Conv_Caps {(i,j) in A}:
    sum {k in 1..nCUT} alpha[k] * xxX[i,j,k] <= y[i,j];