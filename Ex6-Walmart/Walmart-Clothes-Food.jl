include("../src/hypergraph-affinity-functions.jl")

names = ["Clothing, Shoes & Accessories"
"Electronics and Office"
"Home, Furniture & Appliances"
"Home Improvement & Patio"
"Baby"
"Toys, Games, and Video Games"
"Food, Household & Pets"
"Pharmacy, Health & Beauty"
"Sports, Fitness & Outdoors"
"Auto, Tires & Industrial"
"Other"]

dataset = "walmart-trips"
r = 25
H, labels = read_hypergraph_data(dataset,r)

## Categories 1 and 7 make up over 50% of the dataset.
# At a high level, this suggests that a lot of shopping can be broadly
# categorized into grocery shopping or clothes shopping. We take the subset
# of the data corresponding to these two (almost balanced) categories, to obtain
# a two-class hypergraph

CF = findall(x->in(x,[1,7]),labels)
labels = labels[CF]
classes = round.(Int64,(labels .- 1) ./6)
n = length(CF)

H2 = H[:,CF]
order = vec(sum(H2,dims = 2))
nonempty = findall(x->x>0,order)
H = H2[nonempty,:]
# class 1 is food

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
k = 14
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

savefig("Figures/ClothesFood/Affinities_$(k)_Label_$class1.pdf")

## Next
for k = 4:4:12
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


savefig("Figures/ClothesFood/Ratios_$(k)_Label_$class1.pdf")
end

## Group Homophily Index, print usable LaTeX code for building table
include("../src/hypergraph-affinity-functions.jl")
S1,S2  = GHI(H,classes,14)
