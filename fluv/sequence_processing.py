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

    numSeq = len(allSeqs)
    sequence = np.empty((numSeq, 317), dtype=str) # create an empty array of strings
    for ii in range(numSeq):
        for jj in range(317): # truncate at 317 residues to focus on HA1 region
            sequence[ii, jj] = allSeqs[ii][jj]
            if sequence[ii, jj] is "J":
                print("caught a J in new loop")
                sequence[ii, jj] = random.choice(['I', 'L'])
    
    label = np.array(allLabels)
    
    print(sequence)
    print("sequence.shape[0]: " + str(sequence.shape[0]))
    print("sequence.shape[1]: " + str(sequence.shape[1]))
    
    # filtering out residues not included in PAM250 pymsa distance matrix (http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html)
    #for i in range(sequence.shape[0]):
    #    for j in range(sequence.shape[1]):
    #        if sequence[i,j] is "J":
    #            print("caught at J in original loop")
    #            sequence[i,j] = random.choice(['I', 'L'])
    
    return (label, sequence)
