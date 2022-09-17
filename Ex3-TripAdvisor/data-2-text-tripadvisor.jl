using MAT
using Random
using StatsBase
include("../hypergraph-affinity-functions.jl")
D = matread("TripAdvisor_NAM_Europe_H.mat")
H = D["H"]
classes = D["classes"]


##
Edges = incidence2elist(H)

open("hyperedges-tripadvisor.txt","w") do f
    for edge in Edges
        if length(edge) > 1
            for i = 1:length(edge)-1 
                node = edge[i]
                write(f,"$node,")
            end
            node = edge[end]
            write(f,"$node \n")
        end
    end
end

##

open("node-labels-tripadvisor.txt","w") do f
    for i = 1:length(classes)
        node = round(Int64,classes[i])
        write(f,"$node \n")
    end
end

open("names-node-labels-tripadvisor.txt","w") do f
    write(f,"1 North America \n")
    write(f,"0 Europe")
end