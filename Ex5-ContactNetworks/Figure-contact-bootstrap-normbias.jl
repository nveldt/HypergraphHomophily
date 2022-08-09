include("../src/hypergraph-affinity-functions.jl")
using StatsBase
name = "primary"
name = "high"

dataset = "contact-$name-school-classes-gender"
r = 25
H, labels= read_hypergraph_data(dataset,r)

m = [2 3 6 7 9 11 14 15 17]
w = [1 4 5 8 10 12 13 16 18]

keep = findall(x->x<19,labels)
H = H[:,keep]
labels = labels[keep]
n = size(H,2)
classes = zeros(Int64,n)
women = findall(x->in(x,w),labels)
classes[women] .= 1


## Bootstrap
alpha, B1, B2, H1, H2, TR1, TR2, N = Hypergraph_to_Scores(H,classes,5)
tG1, tG2 = normalized_bias(H1, H2, B1, B2)

for k = 2:4
B = 100
T1 = TR1[k,1:k]  # True ratio scores
T0 = TR2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

TG1 = tG1[k,1:k]  # True normalized bias scores
TG0 = tG2[k,1:k]

# Bootstraps
R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0  = Bootstrap_Affinities_NormBias(H,k,B,classes)

E1 = vec(maximum(R1,dims = 1))
E1m = vec(minimum(R1,dims = 1))
E0 = vec(maximum(R0,dims = 1))
E0m = vec(minimum(R0,dims = 1))

# Plot Ratio Scores
lw = 2.5
ms = 8
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300

sc = 1
color1 = :green
color2 = :blue
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline "
p = plot()
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[T0 T0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[T1 T1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
savefig("FinalFigs/Ratios_$(k)_$(name)_Bootstrap.svg")

# Plot just affinity scores
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"

p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

# [M0-2*SE0 M0+2*SE]
# [E0m E0]
# [E1m E1]
p = plot!(1:k,[TA0 TA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TA1 TA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
    savefig("FinalFigs/Affinities_$(k)_$(name)_Bootstrap.svg")


xlab = "Affinity type t"
ylab = "Normalized Bias"
p = plot()
p = Plots.plot([1; k],[0; 0],fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

p = plot!(1:k,[TG0 TG0], fillrange=[MG0-sc*SG0 MG0+sc*SG0], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
p = plot!(1:k,[TG1 TG1], fillrange=[MG1-sc*SG1 MG1+sc*SG1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)
savefig("FinalFigs/NormBias_$(k)_$(name)_Bootstrap.svg")

end
