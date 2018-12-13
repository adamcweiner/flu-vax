def FASTA():
    try:
        f = open('human_aligned_seqs.fa', 'r')
    except IOError:                     
        print("The file does not exist")
        return

    sequences = {}
    
    for line in f:
        if line.startswith('>'):
            #print('\n')
            name = line[1:].rstrip('\n')
            name = name.replace('_', ' ')

            sequences[name] = ''
        else:
            sequences[name] += line.rstrip('\n').rstrip('*')

    for key, value in sequences.items():
        sequences[key] = value[:317]

    #print(str(len(sequences)) + " sequences found")
    print(sequences)

    return sequences

FASTA()
