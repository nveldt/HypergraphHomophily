include("../src/hypergraph-affinity-functions.jl")
using StatsBase

dataset = "congress-bills"
r = 25
H, classes = read_hypergraph_data(dataset,r)
# 1 = democrat
# 2 = republican
m,n = size(H)
classes = classes .- 1
# 0 = democrat = 2 mod 0
# 1 = republican
class1 = "Republican"
class2 = "Democrat"

Hyp = incidence2elist(H)

## Bootstrap
alpha, B1, B2, H1, H2, tR1, tR2, N = Hypergraph_to_Scores(H,classes,15)

for k =  5:5:15
B = 100
TR1 = tR1[k,1:k]  # True affinity scores
TR0 = tR2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

# Bootstraps
R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0 = Bootstrap_Affinities(H,k,B,classes)

E1 = vec(maximum(R1,dims = 1))
E1m = vec(minimum(R1,dims = 1))
E0 = vec(maximum(R0,dims = 1))
E0m = vec(minimum(R0,dims = 1))

## Plot Ratio Scores
lw = 2.5
ms = 8
gfs = 15
tfs = 15
titlesize = 14
s1 = 400
s2 = 300

sc = 2
color1 = :red
color2 = :blue
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline "
p = plot()
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/5:k/5:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TR0 TR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TR1 TR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
savefig("Figures/Ratios_$(k)_Congress_Bootstrap.pdf")

## Plot just affinity scores
lw = 2.5
ms = 8
gfs = 15
tfs = 15
titlesize = 14
s1 = 400
s2 = 300

color1 = :red
color2 = :blue
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"

p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/5:k/5:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TA0 TA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TA1 TA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
savefig("Figures/Affinities_$(k)_Congress_Bootstrap.pdf")
end
