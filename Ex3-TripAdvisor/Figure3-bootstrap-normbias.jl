include("../src/hypergraph-affinity-functions.jl")
using StatsBase
using MAT
D = matread("TripAdvisor_NAM_Europe_H.mat")
H = D["H"]
classes = D["classes"]

class1 = "N. America"
class2 = "Europe"
n = size(H,2)
alpha = sum(classes)/n
# B1 = arbitrary_baselines(k,alpha)
# B0 = arbitrary_baselines(k,1-alpha)
## Bootstrap
alpha, B1, B2, H1, H2, tR1, tR2, N = Hypergraph_to_Scores(H,classes,9)
tG1, tG2 = normalized_bias(H1, H2, B1, B2)

for k =  3:3:9
B = 100
TR1 = tR1[k,1:k]  # True ratio scores
TR0 = tR2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

TG1 = tG1[k,1:k]  # True normalized bias scores
TG0 = tG2[k,1:k]

# Bootstraps
R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0  = Bootstrap_Affinities_NormBias(H,k,B,classes)

## Plot Ratio Scores
lw = 2.5
ms = 8
gfs = 15
tfs = 15
titlesize = 14
s1 = 400
s2 = 300
sc = 2
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TA0 TA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TA1 TA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Affinities_$(k)_TA_Bootstrap.pdf")

## Ratio scores

xlab = "Affinity type t"
ylab = "Affinity / Baseline"
p = plot()
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TR0 TR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TR1 TR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
savefig("FinalFigs/Ratios_$(k)_TA_Bootstrap.pdf")


xlab = "Affinity type t"
ylab = "Normalized Bias"
p = plot()
p = Plots.plot([1; k],[0; 0],fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TG0 TG0], fillrange=[MG0-sc*SG0 MG0+sc*SG0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TG1 TG1], fillrange=[MG1-sc*SG1 MG1+sc*SG1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)
savefig("FinalFigs/NormBias_$(k)_TA_Bootstrap.svg")

end
