import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integ

from data.obs_data import     \
    radio_obs_data,      \
    optical_obs_data,    \
    x_rays_obs_data,     \
    gamma_rays_obs_data, \
    optical_obs_errors,  \
    x_rays_obs_errors,   \
    gamma_rays_obs_errors, \
    optical_archival_data, \
    x_rays_archival_data,  \
    gamma_rays_archival_data

import toml

i    = 17
step = 1

volume = "sphere"
model  = "injection"

assert(volume == "shell"        or volume == "sphere")
assert(model  == "steady_state" or model  == "injection")

photon_data_file   = "./data/saved/photon_{}_{}.tsv".format(volume, model)
lepton_data_file   = "./data/saved/lepton_{}_{}.tsv".format(volume, model)
neutrino_data_file = "./data/saved/neutrino_{}_{}.tsv".format(volume, model)
config_file        = "./data/saved/config_{}_{}.toml".format(volume, model)

# config_file = "configs/{}_{}.toml".format(volume, model)
# config_file = "configs/config_gao.toml".format(volume, model)
config = toml.load(config_file)

gamma = config['general'].get('gamma', 1)
B     = config['general']['magnetic_field']
R     = config['volume']['R']
g_min = config['volume'].get('h', 1)
g_min = config['electrons'].get('gamma_min', 1)
g_c   = config['electrons'].get('break_point', 1)
g_max = config['electrons'].get('gamma_max', 1)

ELECTRON_MASS   = 9.109e-28
ELECTRON_CHARGE = 4.8032045e-10
LIGHT_SPEED     = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED
ELECTRON_RADIUS = ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_ENERGY
THOMSON_CROSS_SECTION = (8 * np.pi / 3) * ELECTRON_RADIUS * ELECTRON_RADIUS

q = 4.8032045e-10
m = ELECTRON_MASS
c = LIGHT_SPEED
h = 6.626e-27

U_b = B * B / (8 * np.pi)
nu_0 = 3 / (4 * np.pi) * B * q / (m*c)

synchrotron_timescale = 4/3 * LIGHT_SPEED * THOMSON_CROSS_SECTION * U_b / ELECTRON_ENERGY

escape_timescale = 1 / LIGHT_SPEED
if   volume == "shell":  escape_timescale *= h * np.pi / 4 # Shell
elif volume == "sphere": escape_timescale *= R * 3     / 4 # Sphere
else: assert(0)

plt.close()
plt.rcdefaults()
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax3 = ax1.twiny()

def tic_function(x):
    x2 = x / 4.13e-15
    return x2

ax1.set_ylim([1e2, 1e30])
ax1.grid()
ax1.set_xlabel("Energy (eV)")
ax1.set_ylabel("$Timescale (s)$")

first_plot = True
while True:
    # photon_data_filename = "./data/photon_data_{:04d}.tsv".format(i)

    photon_data_filename   = photon_data_file

    try:
        data = np.loadtxt(photon_data_filename, unpack=True)
        data_0_eV = data[0] * 511000

        with np.errstate(divide='ignore', over='ignore'):
            ax1.loglog(data_0_eV, - data[2] / data[15], label="$e^- Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[16], label="$p^+ Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[17], label="$\pi^+ Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[18], label="$\pi^- Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[19], label="$\mu^+_L Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[20], label="$\mu^+_R Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[21], label="$\mu^-_L Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[22], label="$\mu^-_R Synchrotron$")
            ax1.loglog(data_0_eV, - data[2] / data[23], label="$IC$")
            ax1.loglog(data_0_eV, - data[2] / data[24], label="$PP$")
            ax1.loglog(data_0_eV, - data[2] / data[25], label="$Escape$")

            ax1.hlines(escape_timescale, *ax1.get_xbound(), label="$Escape$")

            if first_plot:
                ax2.set_xlabel("Frequency (Hz)")
                ax2.xaxis.set_label_coords(0.5,1.05)
                ax2.set_xlim(tic_function(np.array(ax1.get_xbound())))
                ax2.loglog(tic_function(data_0_eV), - data[2] / data[15], 'b')

                ax3.set_xlabel("Normalized Energy")
                ax3.set_xlim(np.array(ax1.get_xbound()) / 511000)
                ax3.loglog(data[0], - data[2] / data[15], 'b')
                ax3.spines['top'].set_position(('outward', 60))

                first_plot = False

    except IOError as e:
        break

    break
    i += step

ax1.legend()
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.tight_layout(0)
ax2.set_xlim(tic_function(np.array(ax1.get_xbound())))
ax3.set_xlim(np.array(ax1.get_xbound()) / 511000)
plt.show()

plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.set_ylim([1e1, 1e13])
ax1.grid()
ax1.set_xlabel("Energy (eV)")
ax1.set_ylabel("$Timescale (s)$")

first_plot = True
while True:
    # lepton_data_filename = "./data/lepton_data_{:04d}.tsv".format(i)

    lepton_data_filename = lepton_data_file

    try:
        data = np.loadtxt(lepton_data_filename, unpack=True)
        data_0_eV = data[0] * 511000

        with np.errstate(divide='ignore', over='ignore'):
            ax1.loglog(data_0_eV, 1 / (synchrotron_timescale * data[0]), label="$Synchrotron$")
            ax1.loglog(data_0_eV, 0.1 / (synchrotron_timescale * data[0]), label="$Synchrotron$")
            ax1.loglog(data_0_eV, - data[1] / data[5], label="$Synchrotron$")
            ax1.loglog(data_0_eV, - data[1] / data[6], label="$IC$")
            ax1.loglog(data_0_eV, - data[1] / data[7], label="$Escape$")

            ax1.hlines(escape_timescale, *ax1.get_xbound(), label="$Escape$")

            if first_plot:
                ax2.set_xlabel("Normalized Energy")
                ax2.xaxis.set_label_coords(0.5,1.05)
                ax2.set_xlim(np.array(ax1.get_xbound()) / 511000)
                ax2.loglog(data[0], - data[1] / data[5], 'b')

                first_plot = False

    except IOError as e:
        break

    break
    i += step

ax1.legend()
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.tight_layout(0)
ax2.set_xlim(np.array(ax1.get_xbound()) / 511000)
plt.show()
