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

# homophily indices for clique expanded hypergraph

hw, hm = clique_expansion_homophily(N,10)

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