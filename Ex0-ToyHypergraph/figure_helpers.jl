include("../src/helper.jl")
using LaTeXStrings
using Plots

function circle(h,k,r)
    th = LinRange(0,2*pi,500)
    h .+ r*sin.(th), k.+ r*cos.(th)
end

function plot_hypergraph(xy,Edges,radius,fsize,lw,al,gverts,bverts)

    k = size(Edges,2)
    n = size(xy,1)
    xy_tuples = Vector{Tuple{Float64,Float64}}()
    for i = 1:size(xy,1)
        push!(xy_tuples,(xy[i,1],xy[i,2]))
    end

    # Hyperedges as triangles
    p = plot(grid = false,legend = false,axis = false)
    for i = 1:size(Edges,1)
        sp = Shape(xy_tuples[Edges[i,:]])
        plot!(sp,color = :black, alpha = al,linewidth = lw)
    end

    # Now vertices
    for i = gverts
    plot!(p,circle(xy[i,1],xy[i,2],radius),linewidth = 0,aspect_ratio = 1,fillalpha = 1, c = :green,seriestype = [:shape])
    annotate!(p,xy[i,1],xy[i,2],text("$i",fsize,:white))
    end

    for i = bverts
    plot!(p,circle(xy[i,1],xy[i,2],radius),linewidth = 0,aspect_ratio = 1,fillalpha = 1, c = :blue,seriestype = [:shape])
    annotate!(p,xy[i,1],xy[i,2],text("$i",fsize,:white))
    end

    return p

end


function plot_bipartite(a,dg,db,dif,Edges,radius,grad,fsize,lw,gverts,bverts)
p1 = plot(grid = false,legend = false,axis = false)


xy =[-a 3*dg;
-a dg;
-a -dg;
-a -3*dg;
a db;
a 0;
a -db
]

m = size(Edges,1)
spn = m*dif/2
xyg =[ zeros(m,1) LinRange(spn,-spn,m)]

# Plot lines for edges
for j = 1:m
    edge = Edges[j,:]
    gp_xy = xyg[j,:]
    for i = 1:3
        nd = edge[i]
        if nd > 4
            cl = :blue
        else
            cl = :green
        end
        nd_xy = xy[edge[i],:]
        plot!(p1,[gp_xy[1];nd_xy[1]],[gp_xy[2];nd_xy[2]], linewidth = lw, linecolor = cl)
    end
end

# Nodes on top
n = size(xy,1)
alphabet = ["a";"b";"c";"d";"e";"f";"g";"h";"i";"j";"k";"l"]
# Plot
xyg_tuples = Vector{Tuple{Float64,Float64}}()
for i = 1:size(xyg,1)
    push!(xyg_tuples,(xyg[i,1],xyg[i,2]))
end
xy_tuples = Vector{Tuple{Float64,Float64}}()
for i = 1:size(xy,1)
    push!(xy_tuples,(xy[i,1],xy[i,2]))
end


# Now vertices
for i = 1:n
    if in(i,gverts)
        plot!(p1,circle(xy[i,1],xy[i,2],radius),linewidth = 0,aspect_ratio = 1,fillalpha = 1, c = :green,seriestype = [:shape])
        annotate!(p1,xy[i,1],xy[i,2],text("$i",fsize,:white))
    else
        plot!(p1,circle(xy[i,1],xy[i,2],radius),linewidth = 0,aspect_ratio = 1,fillalpha = 1, c = :blue,seriestype = [:shape])
        annotate!(p1,xy[i,1],xy[i,2],text("$i",fsize,:white))
    end
end

# Group Vertices
for i = 1:m
plot!(p1,circle(xyg[i,1],xyg[i,2],grad),linewidth = 0,aspect_ratio = 1,fillalpha = 1, c = :gray,seriestype = [:shape])
annotate!(p1,xyg[i,1],xyg[i,2],text("$(alphabet[i])",fsize,:white))
end

return p1
end


function vert_table()
    ## Vertical table
    for i = 1:n
        if i < 5
            print("{\\color{dgreen}$i} & ")
        else
            print("{\\color{blue}$i} & ")
        end
        for j = 1:2
            print("$(d[i,j]) & ")
        end
        print("$(d[n,3]) \\\\\n")
    end
end
