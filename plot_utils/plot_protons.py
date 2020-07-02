import matplotlib.pyplot as plt
import numpy as np

def plot_protons(i):
    data_filename = "../data/hadron_data_{:04d}.tsv".format(i)

    data = np.loadtxt(data_filename, unpack=True)

    norm = data[1] / np.sum(data[3:5],0)

    plt.rcdefaults()

    # plt.xlim([1e4, 1e9])
    plt.ylim([1e-7, 1e5])
    plt.grid()
    plt.xlabel("Normalized Energy")
    plt.ylabel("$p^+$ Population $(1/cm^3)$")

    plt.loglog(data[0],        data[1], label='Population')
    plt.loglog(data[0], norm * data[3], label='Injection')
    plt.loglog(data[0], norm * data[4], label='Multi Resonances')
    plt.loglog(data[0], norm * data[5], label='Direct Resonance')

    plt.gcf().set_tight_layout(True)

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.legend(loc='best')
    plt.show()
