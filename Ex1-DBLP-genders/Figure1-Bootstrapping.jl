using MAT
using Random
include("../src/hypergraph-affinity-functions.jl")
D = matread("DBLP-Author-Gender_H.mat")
H = D["H"]
classes = D["class"]

n = size(H,2)
r = 4
alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,r)
class1 = "female"
class2 = "male"


## Bootstrap
alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,4)

k =  2
B = 100
T1 = R1[k,1:k]  # True affinity scores
T0 = R2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

# Bootstraps
R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0 = Bootstrap_Affinities(H,k,B,classes)

E1 = vec(maximum(R1,dims = 1))
E1m = vec(minimum(R1,dims = 1))
E0 = vec(maximum(R0,dims = 1))
E0m = vec(minimum(R0,dims = 1))


## Plot affinity scores
lw = 2.5
ms = 8
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
sc = 2
color1 = :green
color2 = :blue
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[MA0 MA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[MA1 MA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("Figures/Affinities_$(k)_DBLP_Bootstrap.pdf")

## Ratio scores

xlab = "Affinity type t"
ylab = "Affinity / Baseline"
p = plot()
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

p = plot!(1:k,[MR0 MR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[MR1 MR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
savefig("Figures/Ratios_$(k)_DBLP_Bootstrap.pdf")
