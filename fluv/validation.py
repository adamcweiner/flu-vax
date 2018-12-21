from sequence_processing import trim


class historical_validation:
    def __init__(self):
        """ Load sequences and vaccine efficacy scores according to literature. """
        self.vax_labels, self.vax_seqs = trim("vaccine_strains.fa")
        self.circ_labels, self.circ_seqs = trim("circulating_strains.fa")
        # TODO: add array of measured vaccine efficacies... make sure indexing matches sequence ordering
        
    def calc_PAM(self):
        """"""
    
    
