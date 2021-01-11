include("../src/hypergraph-affinity-functions.jl")
using MAT
using Plots
M = matread("Total_Ns.mat")
Ng = M["Ng"]
Nf = M["Nf"]
Nw = M["Nw"]

M = Nf

## Bootstrap
for k = 2:4

for j = 1:3
if  j == 1
    M = Nf
elseif j == 2
    M = Nw
else
    M = Ng
end
# Set of types of hyperedges
Nk = weights(M[k,1:k+1])

R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0 = Bootstrap_From_N(Nk,B,alpha)


# Plot
sc = 1
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
    yl = [.2,1.7]
elseif k ==3
    yl = [.2,2]
elseif k == 4
    yl = [.1,2.5]
end

if k == 2
    ylab = "Affinity / Baseline"
else
    ylab = ""
end
p = Plots.plot([1; k],[1; 1],yscale = :identity,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs,ylim = yl)
plot!(p,[MR0 MR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color2, label = class1)
plot!(p,1:k,[MR1 MR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color1, label = class2)

savefig("Figures/Ratios_$(k)_$(j)_sc_$sc.pdf")
end
end
