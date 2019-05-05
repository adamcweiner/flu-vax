""" This file performs analysis on 2012 sequences in order to generate embedding plots. """
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from distance_calculation import Distance
from sequence_processing import trim
from sklearn import manifold

def generate_dist_mat(method):
    """ This function generates a distance matrix between all the virus strains.
    method: method for distance calculation. can be PAM250 (hybrid), FLU (hybrid), or Hamming. """
    # this is sloppy code but it'll do for now
    if method == "PAM250_hybrid":
        mat = "PAM250"
        hybrid = True
        ham = False
    elif method == "PAM250":
        mat = "PAM250"
        hybrid = True
        ham = False
    elif method == "FLU_hybrid":
        mat = "FLU"
        hybrid = True
        ham = False
    elif method == "FLU":
        mat = "FLU"
        hybrid = True
        ham = False
    elif method == "Hamming":
        mat = "FLU"  # this doesn't get used
        hybrid = False  # this doesn't get used
        ham = True
    
    D_2012 = Distance(mat)

    # this is the long calculation that generates a full pairwise distance matrix
    distMat_2012 = D_2012.dist_mat("NA_2012_aligned_seq.fa", hybrid=hybrid, hamming=ham)
    return distMat_2012

def generate_MDS(dist_mat, method, out_file):
    """ Runs MDS and saves figures given a specific dist_mat. 
        method represents the method in which the dist_mat was generated.
        out_file is the name of the figure that gets saved."""
    labels_2012, sequences_2012 = trim("NA_2012_aligned_seq.fa")
    numSeq_2012 = sequences_2012.shape[0]
    
    mds = manifold.MDS(n_components=2, max_iter=8000, dissimilarity="precomputed", n_jobs=1)

    results = mds.fit(distMat_2012)
    pos = results.embedding_
    stress = results.stress_
    print('stress: ' +str(stress))

    cmap = mpl.cm.autumn # more options found here: https://matplotlib.org/tutorials/colors/colormaps.html
    autumn_map = plt.get_cmap('autumn')
    color_array = np.zeros((numSeq_2012,4))
    for ii in range(0,numSeq_2012):
        color_array[ii,:] = cmap(ii / float(numSeq_2012))

    pts = plt.scatter(pos[:, 0], pos[:, 1], color=color_array, cmap=autumn_map, s=10, alpha=0.5)
    plt.scatter(pos[1, 0], pos[1, 1], color='b', marker='*', s=50, alpha=0.9, label=('vaccine target: ' + str(labels_2012[1])))
    plt.legend()
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.grid(alpha=0.3)
    plt.title('HA1 sequences of H3N2 virus: North America - 2012 - ' + str(method))
    plt.savefig(out_file)
