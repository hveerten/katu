import numpy as np
import matplotlib.pyplot as plt

import corner

def create_report(probability, parameters, labels, best):
    ndim = parameters.shape[1]

    for dim in range(ndim):
        a = np.c_[probability, parameters[:,dim]]

        a = a[a[:,1].argsort()]

        a = a.T

        a[0] = a[0].cumsum()
        # a[0] /= a[0,-1]

        median      = np.interp(0.50, a[0], a[1])
        quantile_16 = np.interp(0.16, a[0], a[1])
        quantile_84 = np.interp(0.84, a[0], a[1])

        print("{}\t${:.2f}$\t${:.2f}_{{{:+.2f}}}^{{{:+.2f}}}".format(
            labels[dim],
            best[dim],
            median, median-quantile_16, quantile_84-median))

def create_corner(samples, weights, labels, best, levels, filename, extra = lambda x: None):
    if samples is not None:
        print("Creating {}".format(filename))
        plt.rcdefaults()
        fig = corner.corner(
                samples,
                weights=weights,
                labels=labels,
                quantiles=[0.16, 0.5, 0.84],
                show_titles=True,
                truths=best,
                fill_contours=True,
                levels=levels,
                quiet=True,
                contourf_kwargs={'cmap':'viridis','colors':None})

        extra(fig)
        fig.savefig(filename)
        plt.close(fig)
