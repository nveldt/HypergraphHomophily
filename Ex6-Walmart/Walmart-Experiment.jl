include("../src/hypergraph-affinity-functions.jl")

names = [
"Clothing, Shoes & Accessories"
"Electronics and Office"
"Home, Furniture & Appliances"
"Home Improvement & Patio"
"Baby"
"Toys, Games, and Video Games"
"Food, Household & Pets"
"Pharmacy, Health & Beauty"
"Sports, Fitness & Outdoors"
"Auto, Tires & Industrial"
"Other"
]
dataset = "walmart-trips"
r = 25
H, labels = read_hypergraph_data(dataset,r)

class1 = 10
namelabel = names[class1]

n = size(H,2)
classes = zeros(Int64,n)
C1 = findall(x->x==class1,labels)
classes[C1] .= 1

# Turn labels to two classes

alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
hr, hd = clique_expansion_homophily(N,r)
@show alpha
## Plot
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
k = 5
color1 = :blue
color2 = :gray
title = "$k-size hyperedge, $namelabel"
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

savefig("Figures/Label$class1/Affinities_$(k)_Label_$class1.pdf")

## Next
for k = 2:15
s1 = 400
s2 = 300
title = "$k-size hyperedge"
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

@show k
savefig("Figures/Label$class1/Ratios_$(k)_Label_$class1.pdf")
end

## Group Homophily Index, print usable LaTeX code for building table
include("../src/hypergraph-affinity-functions.jl")
S1,S2  = GHI(H,classes,k)
