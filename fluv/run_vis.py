""" Execute this file to generate all the embedding plots from 2012 North America sequences. """

from visualization import generate_dist_mat, generate_MDS, generate_tSNE

methods = ["PAM250_hybrid", "PAM250", "FLU_hybrid", "FLU", "Hamming"]
for method in methods:
    dist_mat = generate_dist_mat(method)
    