# Useful functions for checking proof details of the impossibility proof
using Gurobi
using JuMP
gurobi_env = Gurobi.Env()
"""
Builds the B matrix from the proof
"""
function Bmat(b1,b2)
    k = length(b1)
    @assert(isodd(k))
    r = round(Int64,(k+1)/2)
    b_half = [reverse(b2[r:end]); b1[r:end]]
    Db = diagm(b_half)
    Dk = zeros(k+1,k+1)
    for i = 1:r
        Dk[i,i] = k-i+1
        Dk[k-i+2,k-i+2] = k-i+1
    end
    row = collect(k:-1:0.0)
    r2 = collect(0.0:1:k)
    Top = repeat(row', r)
    Bottom = repeat(r2',r)
    R = [Top; Bottom]
    # @show size(Dk), size(Db), size(R), r, k

    B = Dk - Db*R

    return Dk, Db, R, B

end

"""
Computes primal variables corresponding to a complete hypergraph in the
case of using standard baseline scores for the proof.

These will often be optimal, but not if there is significant imbalance between
the two class sizes
"""
function Xstanard(k,n1,n2)
    n = n1+n2
    Xs = zeros(k+1)
    for i = 1:k+1
        t = i-1
        Xs[i] = binomial(n1,t)*binomial(n2,k-t)/binomial(n,k)
    end

    return Xs
end


"""
Optimal primal solution for but for asymptotic baselines
"""
function Xasym(k,alph)
    Xa = zeros(k+1)
    for i = 1:k+1
        t = i-1
        Xa[i] = alph^t*(1-alph)^(k-t)*binomial(k,t)
    end

    return Xa
end

"""
Use Gurobi to get the solution numerically
"""
function PrimalDualXY(gurobi_env,B,k)

    A = [B -ones(k+1,1)]
    p = [ones(1,k+1) 0]
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env),
        "OutputFlag" => 0, "FeasibilityTol"=>1e-9,"OptimalityTol"=>1e-9))
    @variable(m, x[1:k+2])
    @constraint(m, y, A*x .>= zeros(k+1,1))
    p = [ones(1,k+1) 0]
    @constraint(m, a, p*x .== 1)
    @constraint(m, z, x .>= zeros(k+2,1))

    @objective(m, Max, x[k+2])
    JuMP.optimize!(m)

    X = JuMP.value.(x)
    Y= JuMP.dual.(y)
    gam = X[end]

    return X[1:end-1], Y, gam
end

"""
Get the optimal dual variables for balanced classes.
"""
function BalancedDual(k)
    L = zeros(k+1)
    r = round(Int64,(k+1)/2)
    for i = 1:r
        L[i] = k/(k-i+1)
        L[k-i+2] = k/(k-i+1)
    end
    delta = sum(L)
    Ybal = L/delta
    return Ybal, delta
end

"""
Get optimal dual variables for any class size.
"""
function OptimalDual(b1,b2)
    k = length(b1)
    @assert(isodd(k))
    Ybal, delta = BalancedDual(k)
    r = round(Int64,(k+1)/2)
    s1 = 0
    s2 = 0
    for i = 1:r
        t = k-i+1
        s1 += b2[t]*(k/t - 1)
        s2 += (2-k/t)*b2[t]
    end

    y2k = 2/delta*s1/(1-s2)
    Ydual = zeros(2*r)
    Ydual[1] = y2k
    Ydual[k+1] = 2/delta - y2k

    for i = 2:r
        t = k-i+1
        Ydual[i] = 2/delta*(k/t-1) + (2-k/t)*y2k
        Ydual[k-i+2] = 2*k/(delta*t) - Ydual[i]
    end
    return Ydual
end


"""
For conceptual simplicity, store b and y values in a double-index
fashion
"""
function y_double(b1,b2)
    k = length(b1)
    b = zeros(2,k)
    y = zeros(2,k)
    @assert(isodd(k))

    Y = OptimalDual(b1,b2)
    r = round(Int64,(k+1)/2)

    for i = 1:r
        t = k-i+1
        y[2,t] = Y[i]
        y[1,t] = Y[k-i+2]
    end

    return y
end

# Convert from double indexing to linear
function tt(c,t)
    @assert(c == 1 || c == 2)

    if c == 2
        i = k-t+1
    else
        i = t+1
    end
    return i
end


# Now we have some code for the monotonic LP proof
function monotonic_lp_jump(b1,b2,k)
    """
    Maximize epsilon so that relative affinity scores satisfy
    h(t)/b(t) > epsilon + h(t-1)/b(t-1) for top r = (k+1)/2 scores
    """
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => 0))

    @variable(m, x[1:k+1])  # x[i+1] = proportion of hyperedges of type i
    @variable(m, y)         # variable we maximize to check beta value
    # @variable(m, d1)        # denominator for class 1 affinities
    # @variable(m, d2)        # denominator for class 2 affinities
    for i = 1:k+1
        @constraint(m, x[i] >= 0)
    end
    @constraint(m, 1 == sum(x[i] for i = 1:k+1))

    # Constraint for whether top J are a beta factor above baseline
    r = round(Int64,(k+1)/2)
    for t = r:k
        @constraint(m, (t)*x[t+1]/b1[t] - (t-1)*x[t]/b1[t-1]>= y)
        @constraint(m, (t)*x[k-t+1]/b2[t] - (t-1)*x[k-t+2]/b2[t-1] >= y)
    end

    @objective(m, Max, y)

    JuMP.optimize!(m)
    Y = JuMP.value.(y)
    X = JuMP.value.(x)

    return X,Y
end


function monotonic_LP_B(b1,b2,k)
    B = zeros(k+1,k+1)
    r = round(Int64,(k+1)/2)
    for i = 1:r
        t = k-i+1
        B[i,i] = t/(b2[t])
        B[i,i+1] = -(t-1)/(b2[t-1])
        B[k-i+2,k-i+2] = t/(b1[t])
        B[k-i+2,k-i+1] = -(t-1)/(b1[t-1])
    end
    return B
end



# Now we have some code for the monotonic LP proof
function monotonic_lp_jump_even_k(b1,b2,k)
    """
    Maximize epsilon so that relative affinity scores satisfy
    h(t)/b(t) > epsilon + h(t-1)/b(t-1) for top r = (k+1)/2 scores
    """
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => 0))

    @variable(m, x[1:k+1])  # x[i+1] = proportion of hyperedges of type i
    @variable(m, y)         # variable we maximize to check beta value
    for i = 1:k+1
        @constraint(m, x[i] >= 0)
    end
    @constraint(m, 1 == sum(x[i] for i = 1:k+1))

    # Constraint for whether top J are a beta factor above baseline
    half = round(Int64,k/2)
    for t = r:k
        @constraint(m, (t)*x[t+1]/b1[t] - (t-1)*x[t]/b1[t-1]>= y)
        @constraint(m, (t)*x[k-t+1]/b2[t] - (t-1)*x[k-t+2]/b2[t-1] >= y)
    end

    @objective(m, Max, y)

    JuMP.optimize!(m)
    Y = JuMP.value.(y)
    X = JuMP.value.(x)

    return X,Y
end
