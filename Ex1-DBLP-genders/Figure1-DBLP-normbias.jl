using MAT
using Random
include("../src/hypergraph-affinity-functions.jl")
D = matread("DBLP-Author-Gender_H.mat")
H = D["H"]
classes = D["class"]
order = vec(sum(H,dims = 2))
r = 21
alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,r)
class1 = "female"
class2 = "male"
n = length(classes)

## homophily indices for clique expanded hypergraph, computed in a few different ways

# Relative class sizes
alp1 = sum(classes)/n
alp2 = 1-alp1

maxr = 4    # largest group size we consider in our hypergraph affinity plots

# Weighted clique projection results
gh1, gh2 = clique_expansion_homophily(N,2:r)
gh1_lim, gh2_lim = clique_expansion_homophily(N,2:maxr)

# Unweighted clique projection results
A = H'*H
k = findall(x->x<=maxr,order)
Hr = H[k,:]
gh1_un, gh2_un = clique_expansion_homophily_unweighted(A,classes)
Ar = Hr'*Hr
gh1_unlim, gh2_unlim = clique_expansion_homophily_unweighted(Ar,classes)

## Some plots
using Plots

lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
k = 2
color1 = :green
color2 = :blue
title = "$k-author papers"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline "
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond, markerstrokecolor = color2, label = class2)
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Ratios_$(k)_DBLP.pdf")


## Normalized Ratio
G1, G2 = normalized_bias(H1, H2, B1, B2)

using Plots

lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
k = 2
color1 = :green
color2 = :blue
title = "$k-author papers"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Normalized Bias"
p = Plots.plot([1; k],[0; 0],fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,G2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond, markerstrokecolor = color2, label = class2)
plot!(p,1:k,G1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Normbias_$(k)_DBLP.pdf")
savefig("FinalFigs/Normbias_$(k)_DBLP.svg")


## Compute the graph homophily index from a projection

## Compute pairs in the clique projection, where we project all at once
N = Nf
In1 = 0
In2 = 0
Cross = 0
r = 10
numw = 0
numm = 0
for k = 2:10
    global In1, In2, Cross, numw, numm
    for j = 1:k+1
        i = j-1         # number of women in this type hyperedge
        Hik = N[k,j]    # number of hyperedges of this type
        In1 += binomial(i,2)*Hik    # number of extra edges of this type
        In2 += binomial(k-i,2)*Hik
        Cross += Hik*i*(k-i)
        # println("$k, $i")

        # get the number of women and men in pictures to see if this is roughly even
        numw += Hik*i
        numm += Hik*(k-i)
    end
end

hw = round(2*In1/(2*In1 + Cross),digits = 3)
hm = round(2*In2/(2*In2 + Cross),digits = 3)
println("Clique expansion homophily scores: women = $hw, men = $hm")
NN = numm+numw
rw = numw/NN
rm = numm/NN
println("$rw $rm")