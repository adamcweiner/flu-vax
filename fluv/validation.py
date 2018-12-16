from sequence_processing import trim

vax_labels, vax_seqs = trim("vaccine_strains.fa")
print(vax_labels)
print(vax_seqs)
circ_labels, circ_seqs = trim("circulating_strains.fa")
