using MAT
D = matread("TripAdvisor_NAM_Europe_H.mat")
H = D["H"]
classes = D["classes"]


alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
class1 = "N. America"
class2 = "Europe"

## Plot
using Plots
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
for k = 3:3:9
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
title = "$k-hotel reviews"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,H2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,H1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("Figures/Affinities_$(k)_TA.pdf")
end

## Next
for k = 3:3:9
s1 = 400
s2 = 300
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline"
p = Plots.plot([1; k],[1; 1],yscale = :log10,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)

savefig("Figures/Ratios_$(k)_TA.pdf")
end

## Group Homophily Index
r = 15
S1,S2  = GHI(H,classes,r)
print("\n \$k\$ ")
start = 2
for k = start:13
    print("& $(k) ")
end
print("\\\\\n \\midrule \nN. America GHI ")
for k = start:13
    print("& $(round(Int64,S1[k])) ")
end
print("\\\\ \n Europe GHI")
for k = start:13
    print("& $(round(Int64,S2[k])) ")
end
print("\\\\ \n")


## use this table
# {\Large
# \begin{tabular}{lllllllllllllllll}
# \toprule
# group size $k$ & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 \\
#  \midrule
# N. America GHI & 1 & 1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 & 3 & 4 \\
#  Europe GHI& 1 & 1 & 1 & 1 & 1 & 1 & 2 & 3 & 3 & 3 & 3 \\
# \bottomrule
# \end{tabular}
# }
