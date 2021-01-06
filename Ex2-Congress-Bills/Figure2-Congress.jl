include("../src/hypergraph-affinity-functions.jl")

dataset = "congress-bills"
r = 25
H, classes = read_hypergraph_data(dataset,r)
# 1 = democrat
# 2 = republican
m,n = size(H)
classes = classes .- 1
# 0 = democrat = 2 mod 0
# 1 = republican

alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
class1 = "Republican"
class2 = "Democrat"
hr, hd = clique_expansion_homophily(N,r)

## Plot
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
k = 5
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
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("Figures/Affinities_$(k)_Congress.pdf")


## Next
k = 15
s1 = 400
s2 = 300
title = "$k-person bills"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline"
p = Plots.plot([1; k],[1; 1],yscale = :log10,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/5:k/5:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)

savefig("Figures/Ratios_$(k)_Congress.pdf")


## Group Homophily Index, print usable LaTeX code for building table
S1,S2  = GHI(H,classes,r)
print("\n \$k\$ ")
start = 5
for k = start:20
    print("& $(k) ")
end
print("\\\\\n \\midrule \nRep. GHI ")
for k = start:20
    print("& $(round(Int64,S1[k])) ")
end
print("\\\\ \n Dem. GHI")
for k = start:20
    print("& $(round(Int64,S2[k])) ")
end
print("\\\\ \n")
