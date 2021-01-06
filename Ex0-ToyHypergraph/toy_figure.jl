include("figure_helpers.jl")
include("../src/helpers.jl")

## Edit point embeddings
xy = [
1.3 1.7;
2.5 1.75;
1.2 .7;
2.1 0;
3.6 1.5;
2.6 .8;
3.5 .25
]
xy = [xy[:,2] xy[:,1]]
Edges = [
1 2 3;
2 5 7;
5 6 7;
2 3 4;
4 6 7;
3 4 6;
]
radius = .23
fsize = 24
lw = 2
al = .3
gverts = collect(1:4)
bverts = collect(5:7)
p = plot_hypergraph(xy,Edges,radius,fsize,lw,al,gverts,bverts)
savefig("Figures/ToyHyper.pdf")

## Plot the affinity scores
n = size(xy,1)
Hyp = Vector{Vector{Int64}}()
for i = 1:size(Edges,1)
    push!(Hyp,Edges[i,:])
end
H = elist2incidence(Hyp,n)
classes = zeros(Int64,n)
classes[gverts] .= 1
k = 3
N = get_hyperedge_counts(H,classes)
B1,B2 = n1_n2_baselines(k,sum(classes),n-sum(classes))
Hw, Hm = affinities_from_N(N)
H1 = Hw[k,:]
H2 = Hm[k,:]
R1 = H1./B1
R2 = H2./B2

lw = 2.5
ms = 9
gfs = 15
tfs = 15
titlesize = 14
s1 = 350
s2 = 300
k = 3
color1 = :green
color2 = :blue
title = ""
leg = false
xlab = "Affinity type t"
ylab = "Affinity / Baseline "
p = Plots.plot([1; k],[1; 1],fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab, 1:k, 1:k,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:k,R2[1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :circle, markerstrokecolor = color2)
plot!(p,1:k,R1[1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1)

savefig("Figures/ToyRatioPlot.pdf")


## Table of types degrees, for LaTeX code
d = zeros(Int64,n,3)
for i = 1:n
    ei = findall(x->x==1,H[:,i])
    di = length(ei)
    for j = 1:di
        J = ei[j]
        edge = findall(x->x==1,H[J,:])
        t = sum(classes[edge])
        if classes[i] == 0
            t = k-t
        end
        # @show i,t, edge
        d[i,t]+=1
    end
end

## Sum differently

print("
{\\footnotesize
\\centering
\\begin{tabular}{lc c c c c c c c}
    \\toprule
    ")
print("& ")
for i = 1:n-1
    if i < 5
        print("{\\color{dgreen} \\bf $i} & ")
    else
        print("{\\color{blue} \\bf $i} & ")
    end
end
print("{\\color{blue}\\bf $n}\\\\ \n")
print("\\midrule\n")
for j = 1:4
    if j == 4
        print("\$ d_{v}\$ & ")
        for i = 1:n-1
            print("$(sum(d[i,:])) & ")
        end
        print("$(sum(d[n,:])) \\\\\n")
    else
        print(" type-$j & ")
        for i = 1:n-1
            print("$(d[i,j]) & ")
        end
        print("$(d[n,j]) \\\\\n")
    end
end
sumd = sum(d,dims =2)

print("\\bottomrule \\\\ \n \\end{tabular} \\\\ }")


## Now print out the array of affinity and baseline scores

G = "{\\color{dgreen} \\text{G}}"
B = "{\\color{blue} \\text{B}}"
hg = ["\\textbf{h}_1($G)","\\textbf{h}_2($G)","\\textbf{h}_3($G)"]
bg = ["\\textbf{b}_1($G)","\\textbf{b}_2($G)","\\textbf{b}_3($G)"]
hb = ["\\textbf{h}_1($B)","\\textbf{h}_2($B)","\\textbf{h}_3($B)"]
bb = ["\\textbf{b}_1($B)","\\textbf{b}_2($B)","\\textbf{b}_3($B)"]

Hw = round.(H1,digits = 2)
Bw = round.(B1,digits = 2)
Hm = round.(H2,digits = 2)
Bm = round.(B2,digits = 2)
print("\\centering \n {\\footnotesize \n \$	\\begin{array}{lll} \n")
print("$(hg[1]) = $(Hw[1]) & $(hg[2]) = $(Hw[2]) & $(hg[3]) = $(Hw[3]) \\\\ \n")
print("$(bg[1]) = $(Bw[1]) & $(bg[2]) = $(Bw[2]) & $(bg[3]) = $(Bw[3]) \\\\ \n")
print("$(hb[1]) = $(Hm[1]) & $(hb[2]) = $(Hm[2]) & $(hb[3]) = $(Hm[3])\\\\ \n")
print("$(bb[1]) = $(Bm[1]) & $(bb[2]) = $(Bm[2]) & $(bb[3]) = $(Bm[3])\n")
print("\\end{array} \$ }")
