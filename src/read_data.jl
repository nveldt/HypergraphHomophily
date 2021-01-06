using Parameters

Parameters.@with_kw mutable struct hyp
    """
    A very simple hypergraph composite type, designed to hold a node list N, an edge list E, a degree sequence D,
    """
    N::Vector{Int64}
    E::Dict{Int64, Dict}
    D::Array{Int64, 1} = Array{Int64, 1}()
end

function read_hypergraph_data(dataname::String, maxsize::Int64=25)
    labels = Int64[]
    open("../data/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    n = length(labels)

    E = Dict{Integer, Dict}()
    open("../data/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if length(edge) > maxsize; continue; end
            sz = length(edge)
            if !haskey(E, sz)
                E[sz] = Dict{}()
            end
            E[sz][edge] = 1
        end
    end

    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end

    N = 1:n
    H = hyp(N, E, D)

    Hin, weights = hypergraph2incidence(H)
    return Hin, weights, labels
end
