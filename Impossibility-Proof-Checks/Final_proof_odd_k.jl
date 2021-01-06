using SparseArrays
using LinearAlgebra
using JuMP

include("proof_helper.jl")
include("../src/hypergraph-affinity-functions.jl")
k = 3

r = round(Int64,(k+1)/2)
alph = 0.56
b1 = arbitrary_baselines(k,alph)
b2 = arbitrary_baselines(k,1-alph)

n = 100
n1 = round(Int64,floor(n*alph))
n2 = n - n1
b1, b2 = n1_n2_baselines(k,n1,n2)


Dk, Db, R, B = Bmat(b1,b2)
Xopt, Yopt, gam = PrimalDualXY(gurobi_env,B,k)
Ydual = OptimalDual(b1,b2)
Ybal, delta =  BalancedDual(k)
@assert(norm(Ydual-Yopt) < 1e-15)

## Now we start to check how things hold

V = R'*Db*Yopt
y = y_double(b1,b2)
t = k-1
c = 2

## I know how to recover V
for c = 1:2
    for t = k:-1:r
        r2t = t*sum(b2[i]*y[2,i] for i = r:k) + (k-t)*sum(b1[i]*y[1,i] for i = r:k)
        V2t = V[tt(2,t)]

        @show norm(r2t- V2t)
    end
end

## Now let's manipulate terms
t = k-1
c = 2
yvec[2,t]*t
y2k = Yopt[1]
V2t = V[tt(2,t)]

# ID (1)
s1 = sum(b2[i] for i = r:k) + sum(b1[i] for i = r:k) - b2[r]

C1 = y2k*t*sum((2-k/i)*b2[i] for i = r:k) + y2k*t*sum((2-k/i)*b1[i] for i = r:k)
C2 = y2k*t*sum((2-k/i)*(b2[i]+b1[i]) for i = r:k)


AA = 2*k*t/delta*sum(b2[i]/i for i = r:k)
BB = -2*t/delta*sum(b2[i] for i = r:k)
CC = t*y2k*sum((2-k/i)*b2[i] for i = r:k)
DD = 2*(k-t)/delta*sum(b1[i] for i = r:k)
FF = -k*y2k*sum((2-k/i)*b1[i] for i = r:k)
GG = t*y2k*sum((2-k/i)*b1[i] for i = r:k)
ans1 = AA+BB+CC+DD+FF+GG

## Checking baselines

i = k-1

k/i*b1[i]

(1-alph)^(k-i)*alph^(i-1)*binomial(k,i)



## Four piece
y = y2k
bet2 = sum(b2[i] for i = r:k)
bet1 = sum(b1[i] for i = r:k)
bet2i = sum(b2[i]/i for i = r:k)
bet1i = sum(b1[i]/i for i = r:k)
ans2 = bet2*( -2*t/delta + 2*y*t) + bet2i*(2*k*t/delta - t*y*k) + bet1*(2*(k-t)/delta + (t-k)*y*2) + bet1i*((k-t)*k*y)


## New approach GOOD (1) in notes

ls = (k-t)*sum(b1[i]*(2/delta - (2-k/i)*y) for i = r:k)
rs = (2/delta*(k - t) + (2*t-k)*y) - t*y
r2 = 2/delta*(k-t) + t*y - k*y
r3 = (2/delta - y)*(k-t)

## simplifies to...

ls = sum(b1[i]*(2/delta - (2-k/i)*y) for i = r:k)
rs = (2/delta - y)
r2 = 2/delta*( (1- sum(b2[i] for i = r:k))/( 1 - sum(b2[i]*(2-k/i) for i = r:k)))
l2 = sum(b1[i] for i = r:k)*(2/delta - 2*y) + sum(b1[i]*(k/i)*y for i = r:k)


yguess = 2/delta*(1- sum(b1[i] for i = r:k))/(1-sum(b1[i]*(2-k/i) for i = r:k))

norm(y-yguess)
