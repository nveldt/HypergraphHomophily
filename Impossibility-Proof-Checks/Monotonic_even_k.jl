using SparseArrays
using LinearAlgebra
using JuMP

include("proof_helper.jl")
gurobi_env = Gurobi.Env()
k = 5
r = round(Int64,(k+1)/2)
alph = 0.57
standard = false

if standard
    n = 200
    n1 = round(Int64,floor(n*alph))
    n2 = n - n1
    b1, b2 = n1_n2_baselines(k,n1,n2)
else
    b1 = arbitrary_baselines(k,alph)
    b2 = arbitrary_baselines(k,1-alph)
end

X,eps = monotonic_lp_jump(b1,b2,k)

Dk, Db, R, B = Bmat(b1,b2)

Xs = Xstanard(k,n1,n2)
Xa = Xasym(k,alph)

# Try another way
B = monotonic_LP_B(b1,b2,k)
Xopt, Yopt, eps = PrimalDualXY(gurobi_env,B,k)


if standard
    type = "Standard"
else
    type = "Asymptotic"
end
println("")

if eps> 0
    println("\tMonotonic homophily is possible!")
else
    println("\n\tMonotonic homophily is impossible.")
end


## Check alternate dual variables
x = Xa
x = X
D1 = sum(i*x[i+1] for i = 1:k)
D2 = sum(i*x[k-i+1] for i = 1:k)
r = round(Int64,(k+1)/2)
for t = r:k
    @assert((t)*x[t+1]/b1[t] - (t-1)*x[t]/b1[t-1] >= 0 )
    @assert((t)*x[k-t+1]/b2[t] - (t-1)*x[k-t+2]/b2[t-1] >= 0)
end

## Check more complicated baselines (from other hypergraph)

include("../src/helper.jl")
data = matread("../data/Polblogs_H.mat")

H = data["H"]; classes = vec(round.(Int64,data["labels"]))

## Get scores
r = 13                    # maximum k size we care about
Hw, Hm = relative_affinities(H,classes,r)
N = get_hyperedge_counts(H,classes)
Nmat = N[1:r,1:r+1]
# Get baselines
b1c = Hw[k,1:k]
b2c = Hm[k,1:k]

## Find solution to the LP

k = 11
N = Nmat[k,1:k+1]
N = N.+1
r = round(Int64,(k+1)/2)

b1 = zeros(k)
b2 = zeros(k)
D1 = sum(i*N[i+1] for i = 1:k)
D2 = sum(i*N[k-i+1] for i = 1:k)
for i = 1:k
    b1[i] = i*N[i+1]/D1
    b2[i] = i*N[k-i+1]/D2
end

B = monotonic_LP_B(b1,b2,k)
Xopt, Yopt, eps = PrimalDualXY(gurobi_env,B,k)


## Check dual variables

S = [r/b2[r] -(r-1)/b1[r-1];
    -(r-1)/b2[r-1] r/b1[r]]

S*[D1;D2]
