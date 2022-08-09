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
alpha, B1, B2, H1, H2, R1, R2, N = Hypergraph_to_Scores(H,classes,25)
G1, G2 = normalized_bias(H1, H2, B1, B2)

## Plot affinity scores
for k = 4:4:12


if k == 4
	yt = [10^-1.5, 10^-1, 10^-.5]
elseif k == 8
	yt = [10^-2.5, 10^-2, 10^-1.5, 10^-1, 10^-.5]
else
	yt = [10^-3, 10^-2.5,  10^-2., 10^-1.5, 10^-1, 10^-.5]
end
sc = 2
lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
color1 = :purple
color2 = RGB(218/255,165/255,35/255)
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity"
p = Plots.plot(fontfamily = "helvetica",xaxis = [1,k+.1],yaxis = [minimum(H1[k,1:k]), maximum(H1[k,1:k])+.1],linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, xticks = k/4:k/4:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab, yticks = yt,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,H2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2, label = class2)
plot!(p,1:k,H1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :diamond,markerstrokecolor = color1, label = class1)

savefig("FinalFigs/Affinity_$(k)_Walmart.pdf")
savefig("FinalFigs/Affinity_$(k)_Walmart.svg")
end


## Ratio scores

for k = 4:4:12

    if k == 4
        yt = [10^-1, 10^-.5, 10^0, 10^.5]
    elseif k == 8
        yt = [10^-1, 10^-.5, 10^0, 10^.5, 10^1, 10^1.5]
    else
        yt = [10^-1, 10^0, 10^1, 10^2, 10^3]
    end
    sc = 2
    lw = 2.5
    ms = 9
    gfs = 15
    tfs = 15
    titlesize = 14
    s1 = 350
    s2 = 300
    color1 = :purple
    color2 = RGB(218/255,165/255,35/255)
    title = ""
    leg = false
    xlab = "Affinity type t"
    ylab = "Affinity / Baseline"
    p = Plots.plot([1; k], [1;1], fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = false,grid = false, size = (s1,s2),color = :gray)
    title!(p,title,titlefontsize = titlesize)
    xaxis!(p,xlab, xticks = k/4:k/4:k,tickfontsize=tfs,guidefontsize = gfs)
    yaxis!(p,ylab, yticks = yt,tickfontsize=tfs,guidefontsize = gfs)
    
    plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
        markershape = :circle, markerstrokecolor = color2, label = class2)
    plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
        markershape = :diamond,markerstrokecolor = color1, label = class1)
    
    savefig("FinalFigs/Ratio_$(k)_Walmart.pdf")
    savefig("FinalFigs/Ratio_$(k)_Walmart.svg")
    end

## Next
for k = 4:4:12
    sc = 2
    lw = 2.5
    ms = 9
    gfs = 15
    tfs = 15
    titlesize = 14
    s1 = 350
    s2 = 300
    color1 = :purple
    color2 = RGB(218/255,165/255,35/255)
    title = ""
    leg = false
    xlab = "Affinity type t"
    ylab = "Normalized Bias"
    p = Plots.plot([1; k], [0;0], fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
    title!(p,title,titlefontsize = titlesize)
    xaxis!(p,xlab, xticks = k/4:k/4:k,tickfontsize=tfs,guidefontsize = gfs)
    yaxis!(p,ylab, tickfontsize=tfs,guidefontsize = gfs)
    
    plot!(p,1:k,G2[k,1:k],color = color2, linewidth = lw, markersize = ms,
        markershape = :circle, markerstrokecolor = color2, label = class2)
    plot!(p,1:k,G1[k,1:k],color = color1, linewidth = lw, markersize = ms,
        markershape = :diamond,markerstrokecolor = color1, label = class1)
    
    savefig("FinalFigs/NormBias_$(k)_Walmart.pdf")
    savefig("FinalFigs/NormBias_$(k)_Walmart.svg")
end

## Group Homophily Index, print usable LaTeX code for building table
# include("../src/hypergraph-affinity-functions.jl")
# S1,S2  = GHI(H,classes,14)


## MoHI and MaHI scores
S1,S2  = MaHI(H,classes,r)
T1,T2  = MoHI(H,classes,r)

for k = 2:14
    println("$(k) & $(round(Int64,S1[k])) & $(round(Int64,S2[k])) & $(round(Int64,T1[k])) & $(round(Int64,T2[k])) \\\\")
end


## Check simple homophily

for k = 2:12
    println("$(G1[k,k]) $(G2[k,k])")
end