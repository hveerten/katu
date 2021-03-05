import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import corner

def print_medians(i, samples, ndim):
    print("{:03d}".format(i+1), end='\t')
    for j in range(ndim):
        medians = np.percentile(samples[:,j] , [16,50,84])
        medians[0] = medians[1] - medians[0]
        medians[2] = medians[2] - medians[1]
        print("{:.3f} ± {:.3f}".format(medians[1], (medians[0] + medians[2]) / 2), end='\t')
    print("", end='\n')

def print_best(i, best, best_lnprob):
    print("BEST:", end='\t')
    for b in best:
        print("{:.3f}".format(b), end='\t')
    print("{:.3f}".format(best_lnprob))

# sim_data is a 2,n ndarray
# obs_data is a 3,n ndarray
#
# No unit conversion done internally, so they must be compatible.
# The third component of obs_data is the error
def create_fig(filename, sim_data, obs_data):
    plt.figure(figsize=(24,12), dpi=128)
    plt.ylim([1e-17, 1e-9])
    plt.xlim([1e-6, 1e12])
    plt.grid()

    plt.loglog(sim_data[0], sim_data[1])
    plt.errorbar(obs_data[0], obs_data[1], fmt='o', yerr=obs_data[2])

    plt.savefig(filename, dpi='figure')
    plt.clf()
    plt.close()

def create_corner(filename, samples, labels, truths):
    print("Creating Corner Plot")
    fig = corner.corner(samples,
                        labels=labels,
                        truths=truths,
                        quantiles=[0.16,0.5,0.84],
                        show_titles=True,
                        quiet=True)
    fig.savefig(filename)
    plt.close(fig)

def save_data(filename, data, lnprobs):
    print("Saving Data")
    np.savez(filename, data=data, lnprobs=lnprobs)
