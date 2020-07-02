import matplotlib.pyplot as plt
import numpy as np

def plot_electrons(i):
    data_filename = "../data/lepton_data_{:04d}.tsv".format(i)

    data = np.loadtxt(data_filename, unpack=True)

    norm = data[1] / (data[2] + data[3] + data[4] + data[5] + data[6])

    plt.rcdefaults()
    # plt.rcParams.update({'font.size': 20})

    # plt.xlim([1e4, 1e9])
    plt.ylim([1e-30, 1e2])
    plt.grid()
    plt.xlabel("Normalized Energy")
    plt.ylabel("$e^-$ Population $(1/cm^3)$")

    plt.loglog(data[0],      data[1], label='Population')
    plt.loglog(data[0], norm*data[2], label='Injection')
    plt.loglog(data[0], norm*data[4], label='Pair Production')
    plt.loglog(data[0], norm*data[5], label='Bethe-Heitler')
    plt.loglog(data[0], norm*data[6], label='Muon Decay')

    plt.gcf().set_tight_layout(True)

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.legend(loc='best')
    plt.show()
