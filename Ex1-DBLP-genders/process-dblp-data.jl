## Read in original dblp data and store in nicer form

using JSON

L = readlines("../original-data/author-dblp/dblp_authors.txt")

open("all_author_data.txt", "w") do f
    write(f,"Paper ID, URL, Conference tag, author, author place, gender, gender confidence\n")
end

Paper2id = Dict()               # paper string ID to new interger ID
PaperURLs = Vector{String}()    # from ID to URL
PaperAnum = Vector{Int64}()     # from paper ID to number of authors
AuthorList = Set{String}()

## Process line by line
tagset = Set{String}()
nextID = 1
numpapers = 148251
GenderCounts = zeros(numpapers,3)
total_n = 0
total_w = 0
for j = 1:length(L)
    global nextID, total_n, total_w
    Lpiece = split(L[j],"),(")

    for i = 1:length(Lpiece)

        ln = split(Lpiece[i],",")
        conf = ln[2]
        anum = parse(Int64,ln[3]) + 1
        name = ln[4]
        gender = string(ln[5])

        gp1 = ln[6]
        if occursin(");",gp1)
            gp1 = gp1[1:end-2]
        end
        genp = parse(Float64,gp1)

        id = get(Paper2id,conf,0)
        if id == 0
            # new paper, set author count to at least one, put it in the system
            Paper2id[conf] = nextID
            id = nextID
            nextID += 1
            push!(PaperURLs,conf)
            push!(PaperAnum,1)
        else
            # update the author count
            PaperAnum[id] += 1
        end

        confpieces = split(conf,"/")
        tag = confpieces[2]
        push!(tagset,tag)

        # Update gender count for the paper
        if occursin('M',gender)
            gnum = 1
        elseif occursin('F',gender)
            gnum = 2
        elseif occursin('-',gender)
            gnum = 3
        else
            println("error: $gender")
            @show typeof(gender),typeof('M')
        end
        GenderCounts[id,gnum] += 1

        # Count number of unique authors and unique women
        # excluding the set of authors with unknown gender
        if ~in(name,AuthorList)
            push!(AuthorList,name)
            if gnum != 3
                total_n += 1
            end
            if gnum == 2
                total_w += 1
            end
        end

        open("all_author_data.txt", "a") do f
            write(f,"$id, $conf, $tag, $name, $anum, $gender, $genp\n")
        end
    end
end

## Sort through and count each type of hyperedge
# extract the ones where we know the gender of the authors with high confidence

know = findall(x->x<1,GenderCounts[:,3])
G = GenderCounts[know,1:2]
Kvals = vec(round.(Int64,sum(G,dims=2)))

kmax = maximum(Kvals)

N = zeros(21,22)
for i = 1:size(G,1)
    k = Kvals[i]
    t = round(Int64,G[i,2])  # number of female authors on paper
    N[k,t+1] += 1
end


## Save the hyperedge gender count
# using MAT
# num_f = total_w
# num_m = total_n-total_w
# matwrite("DBLP-Genders-M-F-Unknown.mat", Dict("GenderCounts"=>GenderCounts,
# "N"=>N, "num_f"=>num_f, "num_m"=>num_m))
