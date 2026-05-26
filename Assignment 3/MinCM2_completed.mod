# ------------------------------------------------------------------
# Traffic Assignment Problem with linear volume-delay functions
# s_ij(x) = c_ij + cap_ij*x
# Beckmann objective: integral_0^v s_ij(x) dx = c_ij*v + 0.5*cap_ij*v^2
# ------------------------------------------------------------------

set nusos;
set centr within nusos;
set links within (nusos cross nusos);
set origens within centr;
set destins within centr;
set odpair within (origens cross destins);

set destperorig {i in origens} :=
    setof {(i1,j1) in odpair: i = i1} j1;

param g {odpair} > 0;

# Link coordinates, used to reproduce the costs from the previous assignments
param xc {nusos};
param yc {nusos};

# Constant term of the link travel time
param c {(i,j) in links} :=
    100 + (xc[i]-xc[j])^2 + (2*yc[i]-yc[j])^2;

# Slope of the linear volume-delay function
param dcoef default 0.005;
param cap {links} default dcoef;

# Current linear costs used in the all-or-nothing subproblem
param t0 {links} default 1;

param Tdreta {i in nusos, k in origens} :=
    if i in destperorig[k] then -g[k,i]
    else if i = k then sum {j in destperorig[k]} g[k,j]
    else 0;

node N {i in nusos, k in origens}: net_out = Tdreta[i,k];

arc v_k {(i,j) in links, k in origens} >= 0,
    from N[i,k], to N[j,k];

var v {(i,j) in links} >= 0;

subject to flux_total {(i,j) in links}:
    v[i,j] = sum {k in origens} v_k[i,j,k];

# All-or-nothing assignment with current link costs t0
minimize Vg:
    sum {(i,j) in links} t0[i,j]*v[i,j];

# Direct nonlinear equilibrium problem
minimize Vnl:
    sum {(i,j) in links} (c[i,j]*v[i,j] + 0.5*cap[i,j]*v[i,j]^2);

# ------------------------------------------------------------------
# Restricted Simplicial Decomposition master problem
# ------------------------------------------------------------------

param maxVERT integer >= 1 default 501;
param nVERT integer >= 0 default 0;
param activeVertex {1..maxVERT} binary default 0;
param vertexFlow {(i,j) in links, k in 1..maxVERT} default 0;

var alpha {k in 1..maxVERT} >= 0;
var xbar {(i,j) in links} >= 0;

subject to AlphaSum:
    sum {k in 1..maxVERT} alpha[k] = 1;

subject to AlphaActive {k in 1..maxVERT}:
    alpha[k] <= activeVertex[k];

subject to RSD_link {(i,j) in links}:
    xbar[i,j] = sum {k in 1..maxVERT} alpha[k]*vertexFlow[i,j,k];

minimize RSD_Obj:
    sum {(i,j) in links} (c[i,j]*xbar[i,j] + 0.5*cap[i,j]*xbar[i,j]^2);
