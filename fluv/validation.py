import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sequence_processing import trim
from distance_calculation import Distance, Hamming_dist


class historical_validation:
    def __init__(self):
        """ Load sequences and vaccine efficacy scores according to literature. """
        print("class being initialized")
        _, _ = trim("data/NA_2012_aligned_seq.fa") # test that sequence_processing still works
        print("done trimming NA_2012 sequences")
        self.circ_labels, self.circ_seqs = trim("data/aligned_circulating_strains.fa")
        self.vax_labels, self.vax_seqs = trim("data/aligned_vaccine_strains.fa")
        self.year = np.array([1971, 1972, 1973, 1975, 1984, 1985, 1987, 1989, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2001, 2003])
        self.eff = np.array([7, 15, 11, -3, -6, -2, 17, -5, 59, 38, 25, 45, 28, -17, 34, 43, 55, 12]) # https://doi.org/10.1016/j.vaccine.2006.01.010
        assert(len(self.vax_seqs) == len(self.circ_seqs) == len(self.year) == len(self.eff))
        self.size = len(self.vax_seqs)

    def calc_PAM(self, hybrid=False):
        """ Caclulate distances between circulating strain and vaccine in each given year. Using PAM250 matrix here. """
        D_PAM = Distance("PAM250")

        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = D_PAM.seq_dist(self.vax_seqs[ii], self.circ_seqs[ii], hybrid)

        return distances

    def calc_FLU(self, hybrid=False):
        """ Caclulate distances between circulating strain and vaccine in each given year. Using FLU sub. matrix here."""
        D_FLU = Distance("FLU")

        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = D_FLU.seq_dist(self.vax_seqs[ii], self.circ_seqs[ii], hybrid)

        return distances
    
    def calc_Ham(self):
        """ Calculates distances between circulating and vaccine strain via Hamming distance. Each mistmatch is +1 and each match is 0. """
        distances = np.zeros((self.size))
        for ii in range(self.size):
            distances[ii] = Hamming_dist(self.vax_seqs[ii], self.circ_seqs[ii])
        
        return distances

    def plot(self):
        """ Generates plots using different calc functions. """
        flu_dist = self.calc_FLU()
        pam_dist = self.calc_PAM()
        ham_dist = self.calc_Ham()
        print("Hamming distance:", ham_dist)

        # distance calculations that only reference sub matrix if two characters don't match
        hyb_flu_dist = self.calc_FLU(hybrid=True)
        hyb_pam_dist = self.calc_PAM(hybrid=True)
        print("FLU hybrid", hyb_flu_dist)
        print("PAM hybrid", hyb_pam_dist)

        plt.figure(figsize=(8,5)) # set up figure for plotting, width and height in inches
        cmap_1 = cm.autumn(np.linspace(0, 1, self.size))
        cmap_2 = cm.winter(np.linspace(0, 1, self.size))
        cmap_3 = cm.summer(np.linspace(0, 1, self.size))
        
        # plot standard form of FLU matrix
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.eff , np.log(flu_dist))
        plt.scatter(self.eff, np.log(flu_dist), c=cmap_1, label=r"$\mathrm{standard: R^{2} = }$"+str(round(r_value**2, 3)))
        abline(slope, intercept, plt, cmap_1[0])
        
        # plot hybrid form of FLU matrix
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.eff , np.log(hyb_flu_dist))
        plt.scatter(self.eff, np.log(hyb_flu_dist), c=cmap_2, label=r"$\mathrm{hybrid: R^{2} = }$"+str(round(r_value**2, 3)))
        abline(slope, intercept, plt, cmap_2[0])

        plt.title("FLU Distance Prediction Performance")
        plt.xlabel("Observed Efficacy")
        plt.ylabel("Relative Distance (log scaled)")
        plt.legend()
        plt.savefig('Figures/FLU_validation.pdf')
        
        plt.figure(figsize=(8,5)) # set up figure for plotting, width and height in inches
        plt.scatter(self.eff, pam_dist, c=cmap_1, label="standard")
        plt.scatter(self.eff, hyb_pam_dist, c=cmap_2, label="hybrid")
        plt.title("PAM250 Distance Prediction Performance")
        plt.xlabel("Observed Efficacy")
        plt.ylabel("Relative Distance")
        plt.legend()
        plt.savefig('Figures/PAM_validation.pdf')
        
        plt.figure(figsize=(8,5)) # set up figure for plotting, width and height in inches
        plt.scatter(self.eff, ham_dist, c=cmap_1, label="Ham")
        plt.title("Hamming Distance Prediction Performance")
        plt.xlabel("Observed Efficacy")
        plt.ylabel("Relative Distance")
        plt.legend()
        plt.savefig('Figures/Hamming_validation.pdf')

def abline(slope, intercept, ax, color):
    """Plot a line from slope and intercept"""
    axes = ax.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    ax.plot(x_vals, y_vals, '--', c=color)
