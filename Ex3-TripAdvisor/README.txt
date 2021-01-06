README

Original TripAdvisor data can be obtained from http://times.cs.uiuc.edu/~wang296/Data/. This was converted into a hypergraph where nodes are hotels and a hyperedge indicates a set of hotels reviewed by the same user account. Reviews associated with a non-unique account were filtered out (e.g. reviews by "A TripAdvisor Member". This was finally stored in a hypergraph with two class labels, indicating whether the hotel is in Europe or North America (all other hotels discarded). This hypergraph is stored in TripAdvisor_NAM_Europe_H.mat.

