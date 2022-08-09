include("../src/hypergraph-affinity-functions.jl")
using MAT
using Plots
M = matread("Total_Ns.mat")
Ng = M["Ng"]
Nf = M["Nf"]
Nw = M["Nw"]

N = Nf
r = 10


## Bootstrap
for k = 2:4
B = 100
for j = 1:3
if  j == 1
    M = Nf
elseif j == 2
    M = Nw
else
    M = Ng
end

H1, H2 = affinities_from_N(M)
alpha = 1/2
B1 = baselines(r,alpha)
B2 = baselines(r,1-alpha)
tR1 = H1./B1
tR2 = H2./B2

tG1, tG2 = normalized_bias(H1, H2, B1, B2)


# Set of types of hyperedges
Nk = weights(M[k,1:k+1])

TR1 = tR1[k,1:k]
TR0 = tR2[k,1:k]

TG1 = tG1[k,1:k]  # True normalized bias scores
TG0 = tG2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0  = Bootstrap_From_N_normbias(Nk,B,alpha)


# Plot
sc = 1
lw = 2
ms = 9
s1 = 375
s2 = 300
gfs = 18
tfs = 18
titlesize = 18
color1 = :green
color2 = :blue
title = ""
leg = false
if j == 3
    xlab = "Affinity type t"
else
    xlab = ""
end
if k == 2
    ylab = "Affinity / Baseline"
else
    ylab = ""
end

if k == 2
    yl = [.2,1.7]
elseif k ==3
    yl = [.2,2]
elseif k == 4
    yl = [.1,2.5]
end

p = Plots.plot([1; k],[1; 1],yscale = :identity,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs,ylim = yl)
plot!(p,[TR0 TR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color2, label = class1)
plot!(p,1:k,[TR1 TR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color1, label = class2)

savefig("FinalFigs/Ratios_$(k)_$(j)_sc_$sc.pdf")

if k == 2
    ylab = "Normalized Bias"
else
    ylab = ""
end
p = Plots.plot([1; k],[0; 0],yscale = :identity,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
plot!(p,[TG0 TG0], fillrange=[MG0-sc*SG0 MG0+sc*SG0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color2, label = class1)
plot!(p,1:k,[TG1 TG1], fillrange=[MG1-sc*SG1 MG1+sc*SG1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color1, label = class2)

savefig("FinalFigs/NormBias_$(k)_$(j)_sc_$sc.pdf")
savefig("FinalFigs/NormBias_$(k)_$(j)_sc_$sc.svg")

if k == 2
    ylab = "Affinity"
else
    ylab = ""
end
p = Plots.plot(yscale = :identity,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
plot!(p,[TA0 TA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color2, label = class1)
plot!(p,1:k,[TA1 TA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color1, label = class2)

savefig("FinalFigs/Affinity_$(k)_$(j)_sc_$sc.pdf")

end
end
