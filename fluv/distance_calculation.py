from sequence_processing import trim
from pymsa_material import SubstitutionMatrix, PAM250, FLU_sub
import numpy as np
from sklearn import manifold
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

class Distance:
    def __init__(self, seqFile, subMat):
        self.labels, self.sequences = trim(seqFile)
        self.subMat = subMat
        if subMat is "PAM250":
            self.M = PAM250()
        elif subMat is "FLU":
            self.M = FLU_sub()
        self.numSeq = self.sequences.shape[0]

        
    def seq_dist(self, seq1, seq2):
        if (len(seq1) != len(seq2)):
            print('the sequences are of different length')
            return -1
        else:
            dist = np.zeros((len(seq1)))
            for ii in range(len(seq1)):
                dist[ii] = self.M.get_distance(seq1[ii], seq2[ii])
            avg_dist = np.sum(dist) / 317.0
            
            # PAM250 distances are log-scaled... convert to true value
            if self.subMat is "PAM250":
                exp_dist = np.exp(avg_dist)
            
            # find inverse 'K'->'K' gives a highly positive score
            inv_dist = 1 / exp_dist
            
            return inv_dist
   
    def test_mat(self):
        """ function is the same as "dist_mat()" except that it only looks at first 10 sequences
        in order to get a proof of concept for all my functions before scaling up to the full dataset 
        """
        testMat = np.zeros((1000,1000))
        print('calculating the test distance matrix based on PAM250')
        for i in range(0,1000):
            for j in range(i,1000):
                testMat[i,j] = self.seq_dist(self.sequences[i], self.sequences[j])
                # plug in values for mirror images
                testMat[j,i] = testMat[i,j]

        return testMat
        
    def dist_mat(self):
        distMat = np.zeros((self.numSeq, self.numSeq))
        print('calculating the full distance matrix based on PAM250')
        for i in range(0,self.numSeq):
            for j in range(i,self.numSeq):
                distMat[i,j] = self.seq_dist(self.sequences[i], self.sequences[j])
                distMat[j,i] = distMat[i,j] # plug in mirror image values

        return distMat
    

    
