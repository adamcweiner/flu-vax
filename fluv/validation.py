import numpy as np
from sequence_processing import trim
from distance_calculation import Distance


class historical_validation:
    def __init__(self):
        """ Load sequences and vaccine efficacy scores according to literature. """
        print("class being initialized")
        _, _ = trim("NA_2012_aligned_seq.fa") # test that sequence_processing still works
        print("done trimming NA_2012 sequences")
        self.circ_labels, self.circ_seqs = trim("circulating_strains.fa")
        self.vax_labels, self.vax_seqs = trim("vaccine_strains.fa")
        self.year = np.array([1971, 1972, 1973, 1975, 1984, 1985, 1987, 1989, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2001, 2003])
        self.eff = np.array([7, 15, 11, -3, -6, -2, 17, -5, 59, 38, 25, 45, 28, -17, 34, 43, 55, 12]) # https://doi.org/10.1016/j.vaccine.2006.01.010
        assert(len(self.vax_seqs) == len(self.circ_seqs) == len(self.year) == len(self.eff))
        self.size = len(self.vax_seqs)

    def calc_PAM(self):
        """ Caclulate distances between circulating strain and vaccine in each given year. Using PAM250 matrix here. """
        D_PAM = Distance("PAM250")

        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = D_PAM.seq_dist(self.vax_seqs[ii], self.circ_seqs[ii])
        
        print(distances)
        return distances # this is temporary
        # add plot where x-axis is `self.eff` and y-axis is `distances`... maybe use `circ_labels` as legend

    def calc_FLU(self):
        """ Caclulate distances between circulating strain and vaccine in each given year. Using FLU sub. matrix here."""
        D_FLU = Distance("FLU")

        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = D_FLU.seq_dist(self.vax_seqs[ii], self.circ_seqs[ii])
        
        print(distances)
        print(self.eff)
        return distances # this is temporary
        # add plot where x-axis is `self.eff` and y-axis is `distances`... maybe use `circ_labels` as legend