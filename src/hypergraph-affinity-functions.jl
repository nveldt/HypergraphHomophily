using SparseArrays

function clique_expansion_homophily(N,k_range)
    # Compute homophily index in the clique expansion of the hypergraph
    In1 = 0
    In2 = 0
    Cross = 0
    for k = k_range
        for j = 1:k+1
            i = j-1         # number of women in this type hyperedge
            Hik = N[k,j]    # number of hyperedges of this type
            In1 += binomial(i,2)*Hik    # number of extra edges of this type
            In2 += binomial(k-i,2)*Hik
            Cross += Hik*i*(k-i)
        end
    end

    hw = 2*In1/(2*In1 + Cross)
    hm = 2*In2/(2*In2 + Cross)
    return hw, hm
end

# Used specifically for group picture measures. Computes the 
# proportion of men/women 
function clique_expansion_homophily_proportions(N,k_range,d1=3,d2=2)
    # Compute homophily index in the clique expansion of the hypergraph
    In1 = 0
    In2 = 0
    Cross = 0
    numw = 0
    numm = 0
    for k = k_range
        for j = 1:k+1
            i = j-1         # number of women in this type hyperedge
            Hik = N[k,j]    # number of hyperedges of this type
            In1 += binomial(i,2)*Hik    # number of extra edges of this type
            In2 += binomial(k-i,2)*Hik
            Cross += Hik*i*(k-i)

            # get the number of women and men in pictures to see if this is roughly even
            numw += Hik*i
            numm += Hik*(k-i)
        end
    end

    hw = round(2*In1/(2*In1 + Cross), digits = d1)
    hm = round(2*In2/(2*In2 + Cross), digits = d1)

    NN = numm+numw
    rw = round(numw/NN,digits = d2)
    rm = round(numm/NN,digits = d2)
    return hw, hm, rw, rm
end

# For this version, we use an unweighted clique expansion where two nodes are connected if they share in any hyperedge;
# the edge is not weighted based on the NUMBER of times two nodes are in a hyperedge together.
function clique_expansion_homophily_unweighted(A,classes)
    @assert(maximum(classes) == 1)
    @assert(minimum(classes) == 0)
    @assert(length(unique(classes)) == 2)
    n = size(A,1)
    num1 = 0
    denom1 = 0
    num0 = 0
    denom0 = 0
    I,J,V = findnz(triu(A))

    for t = 1:length(I)
        i = I[t]
        j = J[t]
        if i == j
            continue
        end
        li = classes[i]
        lj = classes[j]
        if li == lj
            if li == 1
                num1 += 2
                denom1 += 2
            else
                num0 += 2
                denom0 += 2
            end
        else
            denom1 += 1
            denom0 += 1
        end
    end

    return num1/denom1, num0/denom0
end

function n_baselines(k,n,alpha)
    n1 = round(Int64,floor(n*alpha))
    n2 = n - n1

    b = zeros(k)
    for t = 1:k
        b[t] = binomial(n1-1,t-1)*binomial(n2,k-t)/binomial(n-1,k-1)
    end
    return b

end

function n1_n2_baselines(k,n1,n2)
    n = n1+n2

    b1 = zeros(k)
    b2 = zeros(k)
    for t = 1:k
        b1[t] = binomial(n1-1,t-1)*binomial(n2,k-t)/binomial(n-1,k-1)
        b2[t] = binomial(n2-1,t-1)*binomial(n1,k-t)/binomial(n-1,k-1)
    end
    return b1, b2
end

function n_alt_baselines(k,n,alpha)
    n1 = round(Int64,floor(n*alpha))
    n2 = n - n1

    b = zeros(k)
    denom = 0
    for t = 1:k
        b[t] = binomial(n1,t)*binomial(n2,k-t)
    end
    return b/sum(b)

end

function arbitrary_baselines(k,alpha)
    """
    Returns asymptotic relative baseline scores for a single class that makes up
        a ratio alpha of the total population. Assumes k-uniform hypergraph.

        k = hyperedge size
        alpha = proportion of nodes in the class of interest
    """
    b = zeros(k)
    for j = 1:k
        b[j] = alpha^(j-1)*(1-alpha)^(k-j)*binomial(k-1,j-1)
    end
    return b
end

function alternative_baselines(k,alpha)
    """
    This is for the alternative definition of affinity score, that doesn't
        degree correct.

    Returns asymptotic relative baseline scores for a single class that makes up
        a ratio alpha of the total population. Assumes k-uniform hypergraph.

        k = hyperedge size
        alpha = proportion of nodes in the class of interest
    """
    b = zeros(k)
    denom = sum(alpha^(i)*(1-alpha)^(k-i)/(factorial(i)*factorial(k-i)) for i = 1:k)
    for j = 1:k
        b[j] = alpha^(j)*(1-alpha)^(k-j)/(denom*(factorial(j)*factorial(k-j)))
    end
    return b
end

function baselines(r,alpha)
    """
    Returns asymptotic relative baseline scores for a single class that makes up
        a ratio alpha of the total population. Does this for all hyperedge sizes
        from 1 to r.
    """
    B = zeros(r,r)
    for k = 1:r
        for t = 1:k
            B[k,t] = alpha^(t-1)*(1-alpha)^(k-t)*binomial(k-1,t-1)
        end
    end
    return B
end

function order_degree(H::SparseMatrixCSC{Float64,Int64})
    """Returns hyperedge sizes and degrees for a hypergraph."""
    order = vec(sum(H,dims=2))
    deg = vec(sum(H,dims=1))
    return order, deg
end

function get_hyperedge_counts(H::SparseMatrixCSC{Float64,Int64},classes::Vector{Int64})
    """
    Compute number of each type of hyperedge. Output is a matrix N, such that
        N[k,j] = number of hyperedge of size k that have j-1 nodes of class 1
    """
    order, degree = order_degree(H)
    r = round(Int64,maximum(order))

    # Ensure classes are 0 and 1
    cs = unique(classes)
    @assert(length(cs) == 2)
    @assert(minimum(cs) == 0 && maximum(cs) == 1)
    Hyp = incidence2elist(H)
    N = zeros(r,r+1)
    for edge in Hyp
        board = classes[edge]
        w = sum(board)
        k = length(board)
        if k > 0
            N[k,w+1] += 1
        end
    end
    return N
end

function relative_affinities(H::SparseMatrixCSC{Float64,Int64},classes::Vector{Int64},r::Int64=0)
    """
    Compute relative affinity scores for both classes in a hypergraph H with
        two classes.

        Aff1[k,t] = type-t relative affinity score for size k hyperedges, for class 1
        Aff0[k,t] = type-t relative affinity score for size k hyperedges, for class 0
    """
    N = get_hyperedge_counts(H,classes)
    if r == 0
        r = size(N,1)   # max hyperedge size
    end
    Aff1 = zeros(r,r)
    Aff0 = zeros(r,r)
    for k = 1:r
        c1_denom = sum(i*N[k,i+1] for i = 1:k)
        c0_denom = sum(i*N[k,k-i+1] for i = 1:k)
        for t = 1:k
            Aff1[k,t] = t*N[k,t+1]/c1_denom
            Aff0[k,t] = t*N[k,k-t+1]/c0_denom
        end
    end
    return Aff1, Aff0
end

function affinities_from_N(N)
    r = size(N,1)
    Aff1 = zeros(r,r)
    Aff0 = zeros(r,r)
    for k = 1:r
        c1_denom = sum(i*N[k,i+1] for i = 1:k)
        c0_denom = sum(i*N[k,k-i+1] for i = 1:k)
        for t = 1:k
            Aff1[k,t] = t*N[k,t+1]/c1_denom
            Aff0[k,t] = t*N[k,k-t+1]/c0_denom
        end
    end
    return Aff1, Aff0
end


function relative_affinities_k(H::SparseMatrixCSC{Float64,Int64},classes::Vector{Int64},k::Int64)
    """
    Compute relative affinity scores for both classes in a hypergraph H with
        two classes, only for hyperedges of size k.

        Aff1[t] = type-t relative affinity score for class 1
        Aff0[t] = type-t relative affinity score for class 0
    """
    N_all = get_hyperedge_counts(H,classes)
    N = N_all[k,:]
    Aff1 = zeros(k)
    Aff0 = zeros(k)
    c1_denom = sum(i*N[i+1] for i = 1:k)
    c0_denom = sum(i*N[k-i+1] for i = 1:k)
    for t = 1:k
        Aff1[t] = t*N[t+1]/c1_denom
        Aff0[t] = t*N[k-t+1]/c0_denom
    end
    return Aff1, Aff0
end

function incidence2elist(Hin::SparseMatrixCSC{Float64,Int64})
    """ Convert a hyperedge list to a hypergraph binary incidence matrix. """
    H = sparse(Hin')
    rp = H.rowval
    ci = H.colptr
    Hyperedges = Vector{Vector{Int64}}()
    n,m = size(H)
    for i = 1:m
        startedge = ci[i]
        endedge = ci[i+1]-1
        edge = rp[startedge:endedge]
        push!(Hyperedges,edge)
    end
    return Hyperedges
end



function elist2incidence(Hyperedges::Vector{Vector{Int64}}, N::Int64)
    """
    Take a list of hyperedges and turn it into a hyperedge incidence matrix
    H. N is the number of nodes in the hypergraph.
    H(e,u) = 1  iff node u is in hyperedge e.
    """
    U = Vector{Int64}()
    E = Vector{Int64}()
    M = length(Hyperedges)
    for enum = 1:length(Hyperedges)
        e = Hyperedges[enum]
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end

    H = sparse(E,U,ones(length(U)),M,N)
    return H
end


function read_hypergraph_data(dataname::String, maxsize::Int64=25)
    """
    Read data from .txt file to incidence matrix, and extract
    class labels.
    """
    classes = Int64[]
    open("../original-data/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(classes, parse(Int64, line))
        end
    end
    n = length(classes)

    # hyperedges
    EdgeList = Vector{Vector{Int64}}()
    open("../original-data/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            push!(EdgeList,edge)
        end
    end

    H = elist2incidence(EdgeList,n)

    return H,classes
end


function Hypergraph_to_Scores(H,classes,r)
    n = size(H,2)
    alpha = sum(classes)/n
    B1 = baselines(r,alpha)
    B2 = baselines(r,1-alpha)
    H1, H2 = relative_affinities(H,classes,r)

    N = get_hyperedge_counts(H,classes)
    N = N[1:r,1:r+1]
    R1 = H1./B1
    R2 = H2./B2

    return alpha, B1, B2, H1, H2, R1, R2, N
end

function MaHI(H,classes,r)
    """
    Majority homophily index: largest value of j such that the top j affinity
    scores are above baseline. 
    """
    alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,r)
    Sets1 = zeros(r)
    Sets2 = zeros(r)
    for k = 1:r
        last1 = k
        last2 = k
        c1searching = true
        c2searching = true

        if R1[k,k] < 1
            c1searching = false
            last1 = k+1
        end

        if R2[k,k] < 1
            c2searching = false
            last2 = k+1
        end

        for j = k:-1:1
            # @assert(R1[k,k] >= 1)
            # @assert(R2[k,k] >= 1)
            if R1[k,j] > 1 && c1searching
                last1 = j
            else
                c1searching = false
            end
            if R2[k,j] > 1 && c2searching
                last2 = j
            else
                c2searching = false
            end
            if R1[k,j] < 1 && R2[k,j] < 1
                break
            end
        end
        Sets1[k] = last1
        Sets2[k] = last2
    end

    ## Set-in plot
    S1 = collect(1:r) - Sets1 .+ 1
    S2 = collect(1:r) - Sets2 .+ 1

    return S1, S2
end


function MoHI(H,classes,r)
    """
    Monotonic homophily index.
    """
    alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,r)
    Sets1 = zeros(r)
    Sets2 = zeros(r)
    for k = 2:r
        last1 = k
        last2 = k
        c1searching = true
        c2searching = true

        if R1[k,k] < R1[k,k-1]
            c1searching = false
            last1 = k+1
        end

        if R2[k,k] < R2[k,k-1]
            c2searching = false
            last2 = k+1
        end

        for j = k:-1:2
            # @assert(R1[k,k] >= 1)
            # @assert(R2[k,k] >= 1)
            if R1[k,j] > R1[k,j-1] && c1searching
                last1 = j
            else
                c1searching = false
            end
            if R2[k,j] > R2[k,j-1] && c2searching
                last2 = j
            else
                c2searching = false
            end
            if R1[k,j] < R1[k,j-1] && R2[k,j] < R2[k,j-1]
                break
            end
        end
        Sets1[k] = last1
        Sets2[k] = last2
    end

    ## Set-in plot
    S1 = collect(1:r) - Sets1 .+ 1
    S2 = collect(1:r) - Sets2 .+ 1

    return S1, S2
end

function Bootstrap_Affinities(H,k,B,classes)
    n = size(H,2)
    R1 = zeros(B,k)
    R0 = zeros(B,k)

    A1 = zeros(B,k)
    A0 = zeros(B,k)

    order = vec(sum(H,dims = 2))
    Ek = findall(x->x==k,order)

    alpha = sum(classes)/n
    B1 = arbitrary_baselines(k,alpha)
    B0 = arbitrary_baselines(k,1-alpha)
    for i = 1:B
        sam = sample(Ek,length(Ek))
        Hnew = H[sam,:]
        Elist = incidence2elist(Hnew)
        Nk = zeros(k+1)
        for edge in Elist
            elabels = classes[edge]
            t = sum(elabels)
            Nk[t+1] += 1
        end
        c1_denom = sum(i*Nk[i+1] for i = 1:k)
        c0_denom = sum(i*Nk[k-i+1] for i = 1:k)
        for t = 1:k
            A1[i,t] = (t*Nk[t+1]/c1_denom)
            A0[i,t] = (t*Nk[k-t+1]/c0_denom)
            R1[i,t] = (t*Nk[t+1]/c1_denom)/B1[t]
            R0[i,t] = (t*Nk[k-t+1]/c0_denom)/B0[t]
        end
    end

    # Average of the random trials
    MR1 = vec(mean(R1,dims = 1))
    MR0 = vec(mean(R0,dims = 1))
    MA1 = vec(mean(A1,dims = 1))
    MA0 = vec(mean(A0,dims = 1))

    # Compute standard error
    SR1 = zeros(k)
    SR0 = zeros(k)
    SA1 = zeros(k)
    SA0 = zeros(k)
    for t = 1:k
        SR1[t] = StatsBase.std(R1[:,t])
        SR0[t] = StatsBase.std(R0[:,t])
        SA1[t] = StatsBase.std(A1[:,t])
        SA0[t] = StatsBase.std(A0[:,t])
    end

    return R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0

end

function Bootstrap_Affinities_NormBias(H,k,B,classes)
    n = size(H,2)
    R1 = zeros(B,k)
    R0 = zeros(B,k)

    A1 = zeros(B,k)
    A0 = zeros(B,k)

    G1 = zeros(B,k)
    G0 = zeros(B,k)

    order = vec(sum(H,dims = 2))
    Ek = findall(x->x==k,order)

    alpha = sum(classes)/n
    B1 = arbitrary_baselines(k,alpha)
    B0 = arbitrary_baselines(k,1-alpha)
    for i = 1:B
        sam = sample(Ek,length(Ek))
        Hnew = H[sam,:]
        Elist = incidence2elist(Hnew)
        Nk = zeros(k+1)
        for edge in Elist
            elabels = classes[edge]
            t = sum(elabels)
            Nk[t+1] += 1
        end
        c1_denom = sum(i*Nk[i+1] for i = 1:k)
        c0_denom = sum(i*Nk[k-i+1] for i = 1:k)
        for t = 1:k
            h1 = (t*Nk[t+1]/c1_denom)
            h0 = (t*Nk[k-t+1]/c0_denom)
            A1[i,t] = h1
            A0[i,t] = h0
            R1[i,t] = h1/B1[t]
            R0[i,t] = h0/B0[t]
            
            if h1 > B1[t]
                G1[i,t] = (h1 - B1[t])/(1 - B1[t])
            else
                G1[i,t] = (h1 - B1[t])/(B1[t])
            end

            if h0 > B0[t]
                G0[i,t] = (h0 - B0[t])/(1 - B0[t])
            else
                G0[i,t] = (h0 - B0[t])/(B0[t])
            end
        end
    end

    # Average of the random trials
    MR1 = vec(mean(R1,dims = 1))
    MR0 = vec(mean(R0,dims = 1))
    MA1 = vec(mean(A1,dims = 1))
    MA0 = vec(mean(A0,dims = 1))

    MG1 = vec(mean(G1,dims = 1))
    MG0 = vec(mean(G0,dims = 1))

    # Compute standard error
    SR1 = zeros(k)
    SR0 = zeros(k)
    SA1 = zeros(k)
    SA0 = zeros(k)
    SG1 = zeros(k)
    SG0 = zeros(k)
    for t = 1:k
        SR1[t] = StatsBase.std(R1[:,t])
        SR0[t] = StatsBase.std(R0[:,t])
        SA1[t] = StatsBase.std(A1[:,t])
        SA0[t] = StatsBase.std(A0[:,t])
        SG1[t] = StatsBase.std(G1[:,t])
        SG0[t] = StatsBase.std(G0[:,t])
    end

    return R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0

end

function Bootstrap_From_N_normbias(Nk_old,B,alpha)

    k = length(Nk_old)-1

    R1 = zeros(B,k)
    R0 = zeros(B,k)
    A1 = zeros(B,k)
    A0 = zeros(B,k)
    G1 = zeros(B,k)
    G0 = zeros(B,k)
    B1 = arbitrary_baselines(k,alpha)
    B0 = arbitrary_baselines(k,1-alpha)

    for i = 1:B
        Nk = zeros(k+1)
        # @show Nk_old
        for j = 1:round(Int64,sum(Nk_old))
            Nk[sample(1:k+1,Nk_old)] += 1
        end
        c1_denom = sum(i*Nk[i+1] for i = 1:k)
        c0_denom = sum(i*Nk[k-i+1] for i = 1:k)
        for t = 1:k
            h1 = (t*Nk[t+1]/c1_denom)
            h0 = (t*Nk[k-t+1]/c0_denom)
            A1[i,t] = h1
            A0[i,t] = h0
            R1[i,t] = h1/B1[t]
            R0[i,t] = h0/B0[t]
            
            if h1 > B1[t]
                G1[i,t] = (h1 - B1[t])/(1 - B1[t])
            else
                G1[i,t] = (h1 - B1[t])/(B1[t])
            end

            if h0 > B0[t]
                G0[i,t] = (h0 - B0[t])/(1 - B0[t])
            else
                G0[i,t] = (h0 - B0[t])/(B0[t])
            end
        end
    end

    # Average of the random trials
    MR1 = vec(mean(R1,dims = 1))
    MR0 = vec(mean(R0,dims = 1))
    MA1 = vec(mean(A1,dims = 1))
    MA0 = vec(mean(A0,dims = 1))

    MG1 = vec(mean(G1,dims = 1))
    MG0 = vec(mean(G0,dims = 1))

    # Compute standard error
    SR1 = zeros(k)
    SR0 = zeros(k)
    SA1 = zeros(k)
    SA0 = zeros(k)
    SG1 = zeros(k)
    SG0 = zeros(k)
    for t = 1:k
        SR1[t] = StatsBase.std(R1[:,t])
        SR0[t] = StatsBase.std(R0[:,t])
        SA1[t] = StatsBase.std(A1[:,t])
        SA0[t] = StatsBase.std(A0[:,t])
        SG1[t] = StatsBase.std(G1[:,t])
        SG0[t] = StatsBase.std(G0[:,t])
    end

    return R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0

end

function Bootstrap_From_N(Nk_old,B,alpha)

    k = length(Nk_old)-1

    R1 = zeros(B,k)
    R0 = zeros(B,k)
    A1 = zeros(B,k)
    A0 = zeros(B,k)
    B1 = arbitrary_baselines(k,alpha)
    B0 = arbitrary_baselines(k,1-alpha)

    for i = 1:B
        Nk = zeros(k+1)
        # @show Nk_old
        for j = 1:round(Int64,sum(Nk_old))
            Nk[sample(1:k+1,Nk_old)] += 1
        end
        c1_denom = sum(i*Nk[i+1] for i = 1:k)
        c0_denom = sum(i*Nk[k-i+1] for i = 1:k)
        for t = 1:k
            A1[i,t] = (t*Nk[t+1]/c1_denom)
            A0[i,t] = (t*Nk[k-t+1]/c0_denom)
            R1[i,t] = (t*Nk[t+1]/c1_denom)/B1[t]
            R0[i,t] = (t*Nk[k-t+1]/c0_denom)/B0[t]
        end
    end

    # Average of the random trials
    MR1 = vec(mean(R1,dims = 1))
    MR0 = vec(mean(R0,dims = 1))
    MA1 = vec(mean(A1,dims = 1))
    MA0 = vec(mean(A0,dims = 1))

    # Compute standard error
    SR1 = zeros(k)
    SR0 = zeros(k)
    SA1 = zeros(k)
    SA0 = zeros(k)
    for t = 1:k
        SR1[t] = StatsBase.std(R1[:,t])
        SR0[t] = StatsBase.std(R0[:,t])
        SA1[t] = StatsBase.std(A1[:,t])
        SA0[t] = StatsBase.std(A0[:,t])
    end

    return R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0

end


function normalized_bias(H1, H2, B1, B2)

    I1 = (H1 - B1) ./ (1 .- B1)
    I2 = (H2 - B2) ./ (1 .- B2)
    J1 = (H1 - B1) ./ (B1)
    J2 = (H2 - B2) ./ (B2)
    
    K,R = size(I1)
    G1 = copy(I1)
    G2 = copy(I2)
    for k = 1:K
        for r = 1:K
            if I1[k,r] < 0
                G1[k,r] = J1[k,r]
            end
            if I2[k,r] < 0
                G2[k,r] = J2[k,r]
            end
        end
    end

    return G1, G2
end