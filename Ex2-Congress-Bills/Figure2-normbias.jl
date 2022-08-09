include("../src/hypergraph-affinity-functions.jl")
using MAT
using Plots

dataset = "congress-bills"
r = 25
H, classes = read_hypergraph_data(dataset,r)
# original labels:
# 1 = democrat
# 2 = republican

m,n = size(H)

# new labels
classes = classes .- 1
# 0 = democrat = 2 mod 0
# 1 = republican

alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
class1 = "Republican"
class2 = "Democrat"
hr, hd = clique_expansion_homophily(N,r)

G1, G2 = normalized_bias(H1, H2, B1, B2)



## Ratios

lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
for k = 5:5:15
color1 = :red
color2 = :blue
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline "
# p = Plots.plot(fontfamily = "Helvetica",yscale = :log10, linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/5:k/5:k, tickfontsize=tfs,guidefontsize = gfs)
if k == 5
    yt = [10^-.2, 10^0, 10^.2, 10^.4, 10^.6]
elseif k == 10
    yt = [10^-.5, 10^0, 10^.5, 10^1, 10^1.5, 10^2]
else
    yt =[10^-1, 10^0, 10^1, 10^2, 10^3, 10^4]
end
yaxis!(p,ylab, yticks = yt, tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Ratios_$(k)_Congress_x.svg")
savefig("FinalFigs/Ratios_$(k)_Congress_x.pdf") # for quick visualization, gets - signs wrong
end

## Plot
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
for k = 5:5:15
color1 = :red
color2 = :blue
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/5:k/5:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,H2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,H1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Affinities_$(k)_Congress_x.pdf")
end

## Normalized bias
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
for k = 5:5:15
color1 = :red
color2 = :blue
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Normalized Bias"
p = Plots.plot([1; k],[0; 0],fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,xticks = k/5:k/5:k,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,G2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,G1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/NormBias_$(k)_Congress.svg")
savefig("FinalFigs/NormBias_$(k)_Congress.pdf")
end



## MoHI and MaHI, Vertical table
S1,S2  = MaHI(H,classes,r)
T1,T2  = MoHI(H,classes,r)

start = 5
for k = start:20
    println("$(k) & $(round(Int64,S1[k])) & $(round(Int64,S2[k])) & $(round(Int64,T1[k])) & $(round(Int64,T2[k])) \\\\")
end

# Use this:
# \begin{tabular}{l |  c c  |c c} 
# 	\toprule
# 	& \multicolumn{2}{c}{MaHI} & \multicolumn{2}{c}{MoHI} \\
# %	\midrule
#  \cmidrule(lr){2-3}\cmidrule(lr){4-5}
# 	$k$ & Re & De & Re & De \\
# 	\midrule
# 5 & \textbf{2} & \textbf{1} & \textbf{2} & \textbf{3} \\
# 6 & 2 & 2 & 3 & 3 \\
# 7 & 2 & 2 & 3 & 4 \\
# 8 & 3 & 2 & 4 & 4 \\
# 9 & 3 & 3 & 5 & 4 \\
# 10 & \textbf{3} & \textbf{3} & \textbf{5} & \textbf{5} \\
# 11 & 4 & 3 & 5 & 6 \\
# 12 & 4 & 4 & 6 & 6 \\
# 13 & 4 & 4 & 6 & 7 \\
# 14 & 5 & 4 & 6 & 8 \\
# 15 & \textbf{5} & \textbf{5} & \textbf{8} & \textbf{7} \\
# 16 & 5 & 5 & 8 & 8 \\
# 17 & 6 & 5 & 9 & 8 \\
# 18 & 7 & 6 & 10 & 8 \\
# 19 & 7 & 6 & 10 & 9 \\
# 20 & 6 & 6 & 9 & 9 \\
# 	\bottomrule
# \end{tabular}

## MaHI horizontal table
S1,S2  = GHI(H,classes,r)
print("\n \$k\$ ")
start = 5
for k = start:20
    print("& $(k) & $(round(Int64,S1[k])) ")
end
print("\\\\\n \\midrule \nRep. MoHI ")
for k = start:20
    print("& $(round(Int64,S1[k])) ")
end
print("\\\\ \n Dem. MoHI")
for k = start:20
    print("& $(round(Int64,S2[k])) ")
end
print("\\\\ \n")
