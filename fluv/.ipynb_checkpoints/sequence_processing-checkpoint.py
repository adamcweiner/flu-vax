from Bio import SeqIO
import os
import numpy as np
import random

def trim(seqFile):
    cwd = os.getcwd()
    homeDir = cwd[1:5]
    if (homeDir == 'home'):
        print('using path from server to load sequences')
        pathToFile = os.path.join("/home","adamw","flu-vax","fluv", str(seqFile)) #aretha server
    else:
        print('using path from windows machine to load sequences')
        pathToFile = os.path.join("C:\\","Users","Adam","Documents","flu-vax", str(seqFile))  #windows machine
    
    allSeqs = []
    allLabels = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
            allSeqs.append(seq_record.seq)
            allLabels.append(seq_record.id)
    
    print((allSeqs))
    
    seqMat = np.array(allSeqs)
    label = np.array(allLabels)
    
    print(seqMat)
    numSeq = seqMat.shape[0]
    print('numSeq = ' +str(numSeq))
    # find label of sequences that are smaller than 317
    for ii in range(numSeq):
        #print(seqMat[ii])
        #print(seqMat[ii, 0])
        if len(seqMat[ii]) < 317:
            print(label[ii])
    
    sequence = seqMat[:, 0:317]
    
    # filtering out residues not included in PAM250 pymsa distance matrix (http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html)
    for i in range(0, sequence.shape[0]):
        for j in range(0,sequence.shape[1]):
            if (sequence[i,j] == 'J'):
                sequence[i,j] = random.choice(['I', 'L'])
    
    return (label, sequence)
