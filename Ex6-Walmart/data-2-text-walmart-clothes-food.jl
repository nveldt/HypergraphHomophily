include("../hypergraph-affinity-functions.jl")
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

##
Edges = incidence2elist(H)

open("hyperedges-groceries-clothes.txt","w") do f
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

open("node-labels-groceries-clothes.txt","w") do f
    for i = 1:length(classes)
        node = round(Int64,classes[i])
        write(f,"$node \n")
    end
end

open("names-node-labels-walmart-groceries-clothes.txt","w") do f
    write(f,"1 Groceries \n")
    write(f,"0 Clothing")
end