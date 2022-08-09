using MAT
D = matread("TripAdvisor_NAM_Europe_H.mat")
H = D["H"]
classes = D["classes"]


alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
class1 = "N. America"
class2 = "Europe"

G1, G2 = normalized_bias(H1, H2, B1, B2)
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
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Affinities_$(k)_TA.pdf")
savefig("FinalFigs/Affinities_$(k)_TA.svg")

end

## Next
for k = 3:3:9
s1 = 350
s2 = 300
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline"
if k == 3
	yt = [10^-.6, 10^-.3, 10^0.0, 10^.3]
elseif k == 6
	yt = [10^-.5, 10^0, 10^0.5, 10^1.0]
else
	yt = [10^-1, 10^0, 10^1, 10^2.0]
end
p = Plots.plot([1; k],[1; 1],yscale = :log10,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab, yticks = yt, tickfontsize=tfs,guidefontsize = gfs)
plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)


savefig("FinalFigs/Ratios_$(k)_TA.svg")
savefig("FinalFigs/Ratios_$(k)_TA.pdf")
end


lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
for k =  3:3:9
    color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Normalized Bias"
p = Plots.plot([1; k],[0; 0],fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,xticks = k/3:k/3:k,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,G2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,G1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/NormBias_$(k)_TA.svg")
savefig("FinalFigs/NormBias_$(k)_TA.pdf")
end


## MoHI and MaHI scores
S1,S2  = MaHI(H,classes,r)
T1,T2  = MoHI(H,classes,r)

start = 2
for k = start:13
    println("$(k) & $(round(Int64,S1[k])) & $(round(Int64,S2[k])) & $(round(Int64,T1[k])) & $(round(Int64,T2[k])) \\\\")
end

# Use this:
# \begin{tabular}{l |  c c  |c c} 
# 	\toprule
# 	& \multicolumn{2}{c}{MaHI} & \multicolumn{2}{c}{MoHI} \\
# %	\midrule
#  \cmidrule(lr){2-3}\cmidrule(lr){4-5}
# 	$k$ & NA & Eu & NA & Eu \\
# 	\midrule
# 2 & 1 & 1 & 1 & 1 \\
# \textbf{3} & \textbf{1} & \textbf{1} & \textbf{2} & \textbf{1} \\
# 4 & 1 & 1 & 3 & 1 \\
# 5 & 1 & 1 & 3 & 2 \\
# \textbf{6} & \textbf{1} & \textbf{1} & \textbf{3} & \textbf{3} \\
# 7 & 2 & 1 & 4 & 3 \\
# 8 & 2 & 1 & 5 & 3 \\
# \textbf{9} & \textbf{2} & \textbf{2} & \textbf{5} & \textbf{4} \\
# 10 & 3 & 3 & 6 & 4 \\
# 11 & 3 & 3 & 6 & 5 \\
# 12 & 3 & 3 & 7 & 5 \\
# 13 & 4 & 3 & 6 & 4 \\
# 	\bottomrule
# \end{tabular}