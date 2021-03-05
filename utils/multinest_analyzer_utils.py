import os
import numpy as np
import matplotlib.pyplot as plt

import scipy.io

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

def extract_points(prefix, ndim, prior):
    points_filename = "{}IS.points".format(prefix)
    points = scipy.io.FortranFile(points_filename, 'r')

    # (ndim + 1) doubles + 1 int
    # parameters + loglikelihood + cluster
    guesstimated_points = os.stat(points_filename).st_size / (8 * (ndim + 1) + 4)

    parameters_record = "({},)f8".format(n_dim + 1)

    all_points = np.zeros((int(guesstimated_points), ndim + 1))

    i = 0
    while True:
        try:
            pp, _ = points.read_record(parameters_record, 'i4')
        except Exception as e:
            break

        if prior is not None:
            all_points[i] = prior(pp.copy())
        else:
            all_points[i] = pp.copy()

        i += 1

        if i % 1000 == 0: print(i)

    all_points.resize((i, n_dim + 1))
    np.savetxt("{}all_data_physical.tsv".format(prefix), all_points)
    np.savez("{}all_data_physical.npz".format(prefix), all_points)
