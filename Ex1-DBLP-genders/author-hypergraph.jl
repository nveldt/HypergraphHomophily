## Put the file in a nicer form

L = readlines("all_author_data.txt")

Author2ID = Dict()
Author2papers = Vector{Vector{Int64}}()
Author2Gender = Vector{Int64}()         # 1 is female, 0 is male

## Process line by line
nextID = 1  # current node id
removelist = Set{Int64}()   # set of papers to remove because we don't have gender data for an author
known = Set{Int64}()        # set of authors for which we have gender data
numpapers = 148251          # number of papers, known from processing data
GenderCounts = zeros(numpapers,3)
total_n = 0
total_w = 0
m = 148251
for j = 2:length(L)
    global nextID
    Lpiece = split(L[j],",")
    author_name = String(strip(Lpiece[4]))
    paperID = parse(Int64,Lpiece[1])
    gender = Lpiece[6]
    if occursin('M',gender)
        gid = 0
    elseif occursin('F',gender)
        gid = 1
    elseif occursin('-',gender)
        push!(removelist,paperID)
        gid = -1
    else
        @show gender
        println("oops")
    end

    id = get(Author2ID,author_name,0)
    if id == 0
        Author2ID[author_name] = nextID
        id = nextID
        nextID += 1

        # if new author, start new paper list
        paperlist = [paperID]
        push!(Author2papers,paperlist)

        # and add their gender
        push!(Author2Gender,gid)

        if gid == 0 || gid == 1
            push!(known,id)
        end
    else
        # add to existing paper list for this person
        push!(Author2papers[id], paperID)
    end
end


## Turn into hypergraph
using SparseArrays
n = nextID-1
I = Vector{Int64}()
J = Vector{Int64}()

for i = 1:n
    papers = Author2papers[i]

    for p in papers
        # put an edge in the bipartite representation for this author-paper
        push!(I,i)
        push!(J,p)
    end

end

B = sparse(J,I,ones(length(I)),m,n)

# Some authors show up in the author list twice, leading to a 2 in an
# entry of B. On two occassions, this is because two authors of the
# same name were on the same paper (often this is avoided by using numbers next
# to author, but evidently there's still some noise). In all 5 other cases,
# I checked the paper, and there's just a mistake in the data. So turning the
# 2 into a 1 is the right thing to do.

## Get rid of papers that have an unknown author on them, and get rid of
# authors of unknown gender. Save author names.
names = Vector{String}()
for j = 1:n
    push!(names,"")
end

for k in keys(Author2ID)
    id = Author2ID[k]
    names[id] = k
    # @show id,k
end

## save the data as a hypergraph

keepedges = setdiff(1:m,removelist)
H = B[keepedges,collect(known)]
m,n = size(H)
I,J,V = findnz(H)
H = sparse(I,J,ones(length(I)),m,n)
name = names[collect(known)]
class = Author2Gender[collect(known)]

matwrite("DBLP-Author-Gender_H.mat", Dict("H"=>H, "name"=>name,"class"=>class))
