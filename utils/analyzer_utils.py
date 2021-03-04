import numpy as np
import matplotlib.pyplot as plt

import corner

def plot_bin_evolution(raw_data, labels, bins):
    walkers, total_iterations, ndim = raw_data.shape
    step = total_iterations // bins

    for dim in range(ndim):
        medians  = []
        averages = []
        for i in range(bins):
            data = raw_data[:, i * step : (i+1) * step, dim]

            medians.append(np.percentile(data, [16, 50, 84]))
            averages.append([np.average(data), np.std(data)])

        medians = np.array(medians).T
        medians[0] = medians[1] - medians[0]
        medians[2] = medians[2] - medians[1]
        averages = np.array(averages).T


        plt.errorbar((np.arange(bins) + 1) / bins, medians[1],  yerr=[medians[0], medians[2]], capsize=4, label=labels[dim])
        plt.errorbar((np.arange(bins) + 1) / bins, averages[0], yerr=averages[1],              capsize=4, label=labels[dim])
        plt.legend()
        plt.show()

def plot_lnprob_evolution(lnprobs, bins):
    if lnprobs is None: return

    walkers, total_iterations = lnprobs.shape
    step = total_iterations // bins

    medians  = []
    averages = []

    # Do not check the first bin because it is rubbish
    for i in range(1, bins):
        data = lnprobs[:, i * step : (i+1) * step]

        medians.append(np.percentile(data, [16, 50, 84]))
        averages.append([np.average(data), np.std(data)])

    medians = np.array(medians).T
    medians[0] = medians[1] - medians[0]
    medians[2] = medians[2] - medians[1]
    averages = np.array(averages).T

    plt.errorbar((np.arange(1, bins) + 1) / bins, medians[1],  yerr=[medians[0], medians[2]], capsize=4, label='lnprob')
    plt.errorbar((np.arange(1, bins) + 1) / bins, averages[0], yerr=averages[1],              capsize=4, label='lnprob')
    plt.legend()
    plt.show()

def plot_random_walkers(raw_data, num):
    walkers, total_iterations, ndim = raw_data.shape

    for dim in range(ndim):
        for i  in np.arange(0, walkers, walkers//num):
            plt.plot(raw_data[i, :, dim])
        plt.show()

def plot_histogram(data, truth, label):
    plt.hist(data, color="#aaaaaa", histtype="step", bins=20)
    plt.axvline(truth, color="#4682b4")

    q_16, q_50, q_84 = np.quantile(data, [0.16, 0.5, 0.84])
    q_m, q_p = q_50-q_16, q_84-q_50

    fmt = "{{0:{0}}}".format(".2f").format
    title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
    title = title.format(fmt(q_50), fmt(q_m), fmt(q_p))
    title = "{0} = {1}".format(label, title)

    plt.title(title)
    plt.show()

def plot_correlation(data1, data2, truth1, truth2, label1, label2, levels):
    if data1 is not None and \
       data2 is not None:

        H, XX, YY = np.histogram2d(data1, data2, 20)

        # X, Y = np.meshgrid(XX, YY)

        X = (XX[:-1] + XX[1:]) / 2
        Y = (YY[:-1] + YY[1:]) / 2

        plt.contourf(X, Y, H.T)
        # plt.hist2d(data1, data2, bins=32)
        # plt.hist2d(data1, data2, 
                # levels=levels,
                # contourf_kwargs={'cmap':'viridis','colors':None})

        if truth1 is not None:
            plt.axvline(truth1, color="#4682b4")
        if truth2 is not None:
            plt.axvline(truth1, color="#4682b4")
        if label1 is not None:
            plt.xlabel(label1)
        if label2 is not None:
            plt.ylabel(label2)

        plt.show()

def create_histogram(data, truth, label, filename):
    if data is not None:
        print("Creating {}".format(filename))
        plt.rcdefaults()
        plt.rcParams.update({'font.size': 20})
        plt.hist(data, density=True, color="k", histtype="step", bins=20)

        if truth is not None:
            plt.axvline(truth, color="#4682b4")

        q_16, q_50, q_84 = np.quantile(data, [0.16, 0.5, 0.84])
        q_m, q_p = q_50-q_16, q_84-q_50

        plt.axvline(q_16, ls="dashed", color="#aaaaaa")
        plt.axvline(q_50, ls="dashed", color="#aaaaaa")
        plt.axvline(q_84, ls="dashed", color="#aaaaaa")

        fmt = "{{0:{0}}}".format(".2f").format
        title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
        title = title.format(fmt(q_50), fmt(q_m), fmt(q_p))
        title = "{0} = {1}".format(label, title)

        plt.title(title)
        plt.savefig(filename)
        plt.close()

def create_correlation(data1, data2, truth1, truth2, label1, label2, levels, filename):
    if data1 is not None and \
       data2 is not None:
        print("Creating {}".format(filename))
        plt.rcdefaults()
        # plt.rcParams.update({'font.size': 20})
        # plt.rcParams.update({'figure.autolayout': True})

        H, XX, YY = np.histogram2d(data1, data2, 20)

        X = (XX[:-1] + XX[1:]) / 2
        Y = (YY[:-1] + YY[1:]) / 2

        plt.contourf(X, Y, H.T)

        if truth1 is not None:
            plt.axvline(truth1, color="#4682b4")
        if truth2 is not None:
            plt.axvline(truth1, color="#4682b4")
        if label1 is not None:
            plt.xlabel(label1)
        if label2 is not None:
            plt.ylabel(label2)

        plt.savefig(filename)
        plt.close()

def create_corner(samples, labels, best, levels, filename, extra = lambda x: None):
    if samples is not None:
        print("Creating {}".format(filename))
        plt.rcdefaults()
        fig = corner.corner(
                samples,
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

def create_report(samples, labels, best):
    data_len, ndim = samples.shape

    for dim in range(ndim):
        data = samples[:,dim]
        medians = np.percentile(data, [16, 50, 84])

        print("{}\t${:.2f}$\t${:.2f}_{{{:+.2f}}}^{{{:+.2f}}}".format(
            labels[dim],
            best[dim],
            medians[1], medians[0]-medians[1], medians[2]-medians[1]))
