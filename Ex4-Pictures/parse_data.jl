sets = ["Fam2a","Fam4a","Fam8a","Group2a","Group4a","Group5a","Wed2a","Wed3a","Wed5a"]

#http://chenlab.ece.cornell.edu/people/Andy/ImagesOfGroups.html

# 1		Female
# 2		Male
for i = 1:length(sets)

f = readlines("../original-data/Group_Pic_Data/$(sets[i])/PersonData.txt")
Men = Vector{Int64}()
People = Vector{Int64}()
r = 10
N = zeros(r,r+1)
for j = 1:length(f)
    global m, p
    if j == 1
        p = 0
        m = 0
        continue
    end
    line = split(f[j],"\t")
    if length(line) == 1
        @assert(occursin("@",string(line)))
        @assert(occursin(".jpg",string(line)))
        # Whenever you get to a new picture
        push!(Men,m)
        push!(People,p)
        if p < 11
            N[p,p-m+1] += 1
        end
        m = 0
        p = 0
    else
        p += 1
        m += parse(Int64,line[end])-1
    end
end
matwrite("N_$(sets[i]).mat",Dict("N"=>N))

end


## Combine them together

wed = ["Wed2a","Wed3a","Wed5a"]
Nw = zeros(10,11)
for i = 1:3
    global Nw
    M = matread("N_$(wed[i]).mat")
    Nw += M["N"]
end


group = ["Group2a","Group4a","Group5a"]
Ng = zeros(10,11)
for i = 1:3
    global Ng
    M = matread("N_$(group[i]).mat")
    Ng += M["N"]
end

fam = ["Fam2a","Fam4a"]
Nf = zeros(10,11)
for i = 1:2
    global Nf
    M = matread("N_$(fam[i]).mat")
    Nf += M["N"]
end

matwrite("Total_Ns.mat",Dict("Nf"=>Nf,"Nw"=>Nw,"Ng"=>Ng))
