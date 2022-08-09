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


## Bootstraps
alpha, B1, B2, H1, H2, tR1, tR2, N = Hypergraph_to_Scores(H,classes,25)
tG1, tG2 = normalized_bias(H1, H2, B1, B2)

for k = 4:4:12
B = 100
TR1 = tR1[k,1:k]  # True ratio scores
TR0 = tR2[k,1:k]

TA1 = H1[k,1:k]  # True affinity scores
TA0 = H2[k,1:k]

TG1 = tG1[k,1:k]  # True normalized bias scores
TG0 = tG2[k,1:k]

# Bootstraps
R0, R1, MR1, MR0, SR1, SR0, A0, A1, MA1, MA0, SA1, SA0, G0, G1, MG1, MG0, SG1, SG0  = Bootstrap_Affinities_NormBias(H,k,B,classes)

## Plot
sc = 2
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 400
s2 = 300
color1 = :purple
color2 = RGB(218/255,165/255,35/255)
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline"
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/4:k/4:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,[TR0 TR0], fillrange=[MR0-sc*SR0 MR0+sc*SR0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,[TR1 TR1], fillrange=[MR1-sc*SR1 MR1+sc*SR1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Ratios_$(k)_Walmart_Food_Clothing.pdf")

## Next
s1 = 400
s2 = 300
title = "$k-size hyperedge"
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = Plots.plot(yscale = :log10,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/4:k/4:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
plot!(p,1:k,[TA1 TA1], fillrange=[MA1-sc*SA1 MA1+sc*SA1],fillalpha=0.3,color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
plot!(p,1:k,[TA0 TA0], fillrange=[MA0-sc*SA0 MA0+sc*SA0], fillalpha=0.3,color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)


savefig("FinalFigs/Affinities_$(k)_Walmart_Food_Clothing_log.pdf")

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
savefig("FinalFigs/NormBias_$(k)_TA_Bootstrap.pdf")

end

