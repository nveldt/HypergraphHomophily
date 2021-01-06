REAME

In order to run the entire experiment from scratch, original data needs to be obtained from https://data.mendeley.com/datasets/3p9w84t5mr/1 and stored in original-data folder. The data is processed in two steps before affinity scores are computed and the results are plotted.

The processed hypergraph is already stored in DBLP-Author-Gender_H.mat, so it is possible to directly run experiments in the last .jl file.

process-dblp-data.jl 
--------------------------
Reads in original data and stores it in a more convenient format


author-hypergraph.jl
--------------------------
Converts data into a hypergraph form, cleans data somewhat and restricts to authors whose gender is known with high confidence

Figure1-DBLP.jl
-------------------------
Computes affinity scores for the hypergraph and plots results
