import matplotlib.pyplot as plt
import numpy as np


plt.close()
plt.rcdefaults()

# cm = plt.cm.get_cmap('OrRd')
cm = plt.cm.get_cmap('tab20')

i = 16

filename = "./data/photon_data_{:04d}.tsv".format(i)

try:
    data = np.loadtxt(filename,unpack=True)
    data_0_eV = data[0] * 511000

    photon_gains = [1e-307] * len(data[0])
    photon_gains += np.sum(data[4:14], 0)

    plt.gca().set_prop_cycle(color=[cm(i /10.) for i in range(10)])

    aux = 1 / photon_gains
    plt.fill_between(data_0_eV, aux * np.sum(data[4:4], 0),  aux * np.sum(data[4:5], 0),  label="$e^- Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:5], 0),  aux * np.sum(data[4:6], 0),  label="$p^+ Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:6], 0),  aux * np.sum(data[4:7], 0),  label="$\pi^+ Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:7], 0),  aux * np.sum(data[4:8], 0),  label="$\pi^- Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:8], 0),  aux * np.sum(data[4:9], 0),  label="$\mu^+_L Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:9], 0),  aux * np.sum(data[4:10], 0), label="$\mu^+_R Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:10], 0), aux * np.sum(data[4:11], 0), label="$\mu^-_L Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:11], 0), aux * np.sum(data[4:12], 0), label="$\mu^-_R Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:12], 0), aux * np.sum(data[4:13], 0), label="$IC\ up$")
    plt.fill_between(data_0_eV, aux * np.sum(data[4:13], 0), aux * np.sum(data[4:14], 0), label="$IC\ down$")

    plt.xscale('log')
    plt.grid()

    plt.tight_layout(0)
    plt.legend()
    plt.show()

    photon_losses = [1e-307] * len(data[0])
    photon_losses += np.sum(data[14:25], 0)

    plt.gca().set_prop_cycle(color=[cm(i /11.) for i in range(11)])

    aux = 1 / photon_losses
    plt.fill_between(data_0_eV, aux * np.sum(data[14:14], 0), aux * np.sum(data[14:15], 0),  label="$e^- Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:15], 0), aux * np.sum(data[14:16], 0),  label="$p^+ Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:16], 0), aux * np.sum(data[14:17], 0),  label="$\pi^+ Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:17], 0), aux * np.sum(data[14:18], 0),  label="$\pi^- Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:18], 0), aux * np.sum(data[14:19], 0),  label="$\mu^+_L Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:19], 0), aux * np.sum(data[14:20], 0), label="$\mu^+_R Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:20], 0), aux * np.sum(data[14:21], 0), label="$\mu^-_L Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:21], 0), aux * np.sum(data[14:22], 0), label="$\mu^-_R Synchrotron$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:22], 0), aux * np.sum(data[14:23], 0), label="$IC$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:23], 0), aux * np.sum(data[14:24], 0), label="$Pair\ Production$")
    plt.fill_between(data_0_eV, aux * np.sum(data[14:24], 0), aux * np.sum(data[14:25], 0), label="$Escape$")

    plt.xscale('log')
    plt.grid()

    plt.tight_layout(0)
    plt.legend()
    plt.show()

except IOError as e:
    raise e

filename = "./data/lepton_data_{:04d}.tsv".format(i)
print(filename)

try:
    data = np.loadtxt(filename,unpack=True)
    data_0 = data[0]

    data[np.isnan(data)] = 1e-307

    lepton_gains = [1e-307] * len(data[0])
    lepton_gains += np.sum(data[2:5], 0)

    aux = 1 / lepton_gains
    plt.fill_between(data_0, aux * np.sum(data[2:2], 0),  aux * np.sum(data[2:3], 0),  label="$External\ Injection$")
    plt.fill_between(data_0, aux * np.sum(data[2:3], 0),  aux * np.sum(data[2:4], 0),  label="$Acceleration$")
    plt.fill_between(data_0, aux * np.sum(data[2:4], 0),  aux * np.sum(data[2:5], 0),  label="$Pair\ Production$")

    plt.xscale('log')
    plt.grid()
    plt.tight_layout(0)
    plt.legend()
    plt.show()

    lepton_losses = [1e-307] * len(data[0])
    lepton_losses += np.sum(data[5:8], 0)

    aux = 1 / lepton_losses
    plt.fill_between(data_0, aux * np.sum(data[5:5], 0),  aux * np.sum(data[5:6], 0),  label="$Synchrotron\ Cooling$")
    plt.fill_between(data_0, aux * np.sum(data[5:6], 0),  aux * np.sum(data[5:7], 0),  label="$IC$")
    plt.fill_between(data_0, aux * np.sum(data[5:7], 0),  aux * np.sum(data[5:8], 0),  label="$Escape$")

    plt.xscale('log')
    plt.grid()
    plt.tight_layout(0)
    plt.legend()
    plt.show()

except IOError as e:
    raise e
