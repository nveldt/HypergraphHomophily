using MAT
using Random
using StatsBase
include("../src/hypergraph-affinity-functions.jl")
D = matread("DBLP-Author-Gender_H.mat")
H = D["H"]
classes = D["class"]


Edges = incidence2elist(H)

open("hyperedges-dblp.txt","w") do f
    for edge in Edges
        for i = 1:length(edge)-1 
            node = edge[i]
            write(f,"$node,")
        end
        node = edge[end]
        write(f,"$node \n")
    end
end

##

open("node-labels-dblp.txt","w") do f
    for i = 1:length(classes)
        node = round(Int64,1-(classes[i]-1))
        write(f,"$node \n")
    end
end

open("names-node-labels-dblp.txt","w") do f
    write(f,"1 female \n")
    write(f,"2 male")
end