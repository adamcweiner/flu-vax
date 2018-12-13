from Bio import SeqIO
import os
import numpy as np
import random

def trim():
    cwd = os.getcwd()
    homeDir = cwd[1:5]
    if (homeDir == 'home'):
        print('hi')
        pathToFile = os.path.join("/home","adamw","final-project-adamcweiner","time_aligned_seqs.fa") #aretha server
    else:
        print('hello')
        pathToFile = os.path.join("C:\\","Users","Adam","Documents","final-project-adamcweiner", "time_aligned_seqs.fa")  #windows machine
    
    #try:
    #    print('hello')
    #    pathToFile = os.path.join("C:\\","Users","Adam","Documents","final-project-adamcweiner", "time_aligned_seqs.fa")  #windows machine
    #except:
    #    print('hi')
    #    pathToFile = os.path.join("home","adamw","final-project-adamcweiner","time_aligned_seqs.fa") #aretha server
        
    
    allSeqs = []
    allLabels = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
            allSeqs.append(seq_record.seq)
            allLabels.append(seq_record.id)
    
    seqMat = np.array(allSeqs)
    label = np.array(allLabels)
    
    sequence = seqMat[:, 0:317]
    
    # filtering out residues not included in PAM250 pymsa distance matrix (http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html)
    for i in range(0, sequence.shape[0]):
        for j in range(0,sequence.shape[1]):
            if (sequence[i,j] == 'J'):
                sequence[i,j] = random.choice(['I', 'L'])
            elif (sequence[i,j] == 'B'):
                sequence[i,j] = random.choice(['D', 'N'])
            elif (sequence[i,j] == 'Z'):
                sequence[i,j] = random.choice(['E', 'Q'])
    
    return (label, sequence)
