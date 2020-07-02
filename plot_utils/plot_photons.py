import matplotlib.pyplot as plt
import numpy as np

def plot_photons(i):
    data_filename = "../data/photon_data_{:04d}.tsv".format(i)

    data = np.loadtxt(data_filename, unpack=True)

    norm = data[0]**2 * data[2] / np.sum(data[3:16],0)

    plt.rcdefaults()

    # plt.xlim([1e4, 1e9])
    plt.ylim([1e-6, 1e3])
    plt.grid()
    plt.xlabel("Normalized Energy")
    plt.ylabel("$\\nu F_{\\nu}$ $(Erg/cm^2 \\times s)$")

    plt.loglog(data[0], data[0]**2 * data[ 2], label='Population')
    plt.loglog(data[0],       norm * data[ 3], label='Injection')
    plt.loglog(data[0],       norm * data[ 4], label='$e^-$ Synchrotron')
    plt.loglog(data[0],       norm * data[ 5], label='$e^+$ Synchrotron')
    plt.loglog(data[0],       norm * data[ 6], label='$p^+$ Synchrotron')
    plt.loglog(data[0],       norm * data[13], label='Inverse Compton (up)')
    plt.loglog(data[0],       norm * data[14], label='Inverse Compton (down)')
    plt.loglog(data[0],       norm * data[15], label='$\\pi^0$ Decay')
    plt.loglog(data[0],       norm * data[16], label='Pair Annihilation')

    plt.gcf().set_tight_layout(True)

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.legend(loc='best')
    plt.show()
