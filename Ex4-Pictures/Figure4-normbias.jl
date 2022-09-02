include("../src/hypergraph-affinity-functions.jl")
using MAT

using Plots
M = matread("Total_Ns.mat")
Ng = M["Ng"]   
Nf = M["Nf"]
Nw = M["Nw"]
r = 10

## Run all
for j = 1:3
for k = 2:4
if  j == 1
    N = Nf
elseif j == 2
    N = Nw
else
    N = Ng
end

alpha = 1/2
Hw, Hm = affinities_from_N(N)
alpha = 1/2
Bw = baselines(r,alpha)
Bm = baselines(r,1-alpha)
R1 = Hw./Bw
R2 = Hm./Bm
G1, G2 = normalized_bias(Hw, Hm, Bw, Bm)


## Counting balance between men and women

N = Nw
T = 2*sum(N[2,:]) + 3*sum(N[3,:]) + 4*sum(N[4,:])
women = N[2,2] + 2*N[2,3] + N[3,2] + 2*N[3,3] + 3*N[3,4] + N[4,2] + 2*N[4,3] + 3*N[4,4] + 4*N[4]
@show women/T

## Plot
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
plot!(p,1:k,R1[k,1:k],color = color1, linewidth = lw, markersize = ms,
    markershape = :circle,markerstrokecolor = color1, label = class1)
plot!(p,1:k,R2[k,1:k],color = color2, linewidth = lw, markersize = ms,
    markershape = :diamond, markerstrokecolor = color2, label = class2)

savefig("FinalFigs/Ratios_$(k)_$(j).pdf")

end

end

for j = 1:3
    for k = 2:4
    if  j == 1
        N = Nf
    elseif j == 2
        N = Nw
    else
        N = Ng
    end
    
    alpha = 1/2
    Hw, Hm = affinities_from_N(N)
    alpha = 1/2
    Bw = baselines(r,alpha)
    Bm = baselines(r,1-alpha)
    R1 = Hw./Bw
    R2 = Hm./Bm
    G1, G2 = normalized_bias(Hw, Hm, Bw, Bm)

    ## Counting balance between men and women
    
    N = Nw
    T = 2*sum(N[2,:]) + 3*sum(N[3,:]) + 4*sum(N[4,:])
    women = N[2,2] + 2*N[2,3] + N[3,2] + 2*N[3,3] + 3*N[3,4] + N[4,2] + 2*N[4,3] + 3*N[4,4] + 4*N[4]
    @show women/T
    
    ## Plot
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
        ylab = "Normalized Bias"
    else
        ylab = ""
    end
    p = Plots.plot([1; k],[0;0],yscale = :identity,fontfamily = "helvetica",linewidth = 2,legend = false,grid = false, size = (s1,s2),color = :gray)
    title!(p,title,titlefontsize = titlesize)
    xaxis!(p,xlab, xticks = 1:k,tickfontsize=tfs,guidefontsize = gfs)
    yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
    plot!(p,1:k,G1[k,1:k],color = color1, linewidth = lw, markersize = ms,
        markershape = :circle,markerstrokecolor = color1, label = class1)
    plot!(p,1:k,G2[k,1:k],color = color2, linewidth = lw, markersize = ms,
        markershape = :diamond, markerstrokecolor = color2, label = class2)
    
    savefig("FinalFigs/Normbias$(k)_$(j).pdf")
    
end

end

## Compute pairs in the clique projection, where we project all at once
N = Nw
In1 = 0
In2 = 0
Cross = 0
r = 10
numw = 0
numm = 0
for k = 2:r
    global In1, In2, Cross, numw, numm
    for j = 1:k+1
        i = j-1         # number of women in this type hyperedge
        Hik = N[k,j]    # number of hyperedges of this type
        In1 += binomial(i,2)*Hik    # number of extra edges of this type
        In2 += binomial(k-i,2)*Hik
        Cross += Hik*i*(k-i)
        # println("$k, $i")

        # get the number of women and men in pictures to see if this is roughly even
        numw += Hik*i
        numm += Hik*(k-i)
    end
end

hw = round(2*In1/(2*In1 + Cross),digits = 3)
hm = round(2*In2/(2*In2 + Cross),digits = 3)
println("Clique expansion homophily scores: women = $hw, men = $hm")
NN = numm+numw
rw = round(numw/NN,digits = 2)
rm = round(numm/NN,digits = 2)
println("Baseline scores: $rw $rm")
println("\\midrule \n 2-$r & $hw & $hm & $rw\\% & $rm\\% \\\\")


## Compute pairs in the clique projection, where we project individual hyperedge sizes
N = Nw
r = 10
println("Clique expansion homophily scores")
for k = 2:4
    In1 = 0
    In2 = 0
    numw = 0
    numm = 0
    Cross = 0
    for j = 1:k+1
        i = j-1         # number of women in this type hyperedge
        Hik = N[k,j]    # number of hyperedges of this type
        In1 += binomial(i,2)*Hik    # number of extra edges of this type
        In2 += binomial(k-i,2)*Hik
        Cross += Hik*i*(k-i)
        
        # get the number of women and men in pictures to see if this is roughly even
        numw += Hik*i
        numm += Hik*(k-i)
    end
    hw = round(2*In1/(2*In1 + Cross),digits = 3)
    hm = round(2*In2/(2*In2 + Cross),digits = 3)
    # println("Clique expansion homophily scores: women = $hw, men = $hm")
    # println("$k \t women = $hw \t men = $hm, nums = ($numw, $numm)")
    # rw = numw/(numw+numm)*100
    NN = numw+numm
    rw = round(numw/NN,digits = 2)
    rm = round(numm/NN,digits = 2)
    # println("\\midrule \n $k & $hw & $hm & $rw \\% \\\\")
    println("\\midrule \n $k & $hw & $hm & $rw & $rm \\\\")

end


## Final table, more concise

for k = 2:4
    print("\\midrule \n $k ")
    for j = 1:3
        if  j == 1
            N = Nf
        elseif j == 2
            N = Nw
        else
            N = Ng
        end
        hw,hm,rw,rm = clique_expansion_homophily_proportions(N,k)
        print("& $hw & $hm & $rw / $rm ")
    end
    println("\\\\")
end
for t = [4 10]
print("\\midrule \n 2-$t ")
for j = 1:3
    if  j == 1
        N = Nf
    elseif j == 2
        N = Nw
    else
        N = Ng
    end
    hw,hm,rw,rm = clique_expansion_homophily_proportions(N,2:t)
    print("& $hw & $hm & $rw / $rm ")
end
println("\\\\")
end

