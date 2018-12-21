from sequence_processing import trim
from distance_calculation import Distance


class historical_validation:
    def __init__(self):
        """ Load sequences and vaccine efficacy scores according to literature. """
        self.vax_labels, self.vax_seqs = trim("vaccine_strains.fa")
        self.circ_labels, self.circ_seqs = trim("circulating_strains.fa")
        # TODO: add array of measured vaccine efficacies... make sure indexing matches sequence ordering
        assert(len(vax_seqs) == len(circ_seqs))
        self.size = len(vax_seqs)
        
    def calc_PAM(self):
        """ Caclulate distances between circulating strain and vaccine in each given year. """
        D_PAM = Distance("PAM250")
        
        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = 
    
    
