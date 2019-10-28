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

i    = 121
step = 1

volume = "sphere"
model  = "injection"

assert(volume == "shell"        or volume == "sphere")
assert(model  == "steady_state" or model  == "injection")

photon_data_file   = "./data/saved/photon_{}_{}.tsv".format(volume, model)
neutrino_data_file = "./data/saved/neutrino_{}_{}.tsv".format(volume, model)
config_file        = "./data/saved/config_{}_{}.toml".format(volume, model)

# config_file = "configs/{}_{}.toml".format(volume, model)
config_file = "configs/config_gao.toml"
# config_file = "configs/config_cerruti_steady_state.toml"
# config_file = "configs/config_cerruti_injection.toml"
config = toml.load(config_file)

gamma = config['general'].get('gamma', 1)
B     = config['general']['magnetic_field']
R     = config['volume']['R']
g_min = config['electrons'].get('gamma_min', 1)
g_c   = config['electrons'].get('break_point', 1)
g_max = config['electrons'].get('gamma_max', 1)

# print(config)

ELECTRON_MASS = 9.109e-28
LIGHT_SPEED   = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED

q = 4.8032045e-10
m = ELECTRON_MASS
c = LIGHT_SPEED
h = 6.626e-27

z = 0.3365
H_0 = 67.66
Omega_m = 0.3111
Omega_l = 0.6889

D_l = (1 + z) * LIGHT_SPEED / (H_0 * 1e5) * integ.quad(lambda x: 1 / np.sqrt((1 + x)**3 * Omega_m + Omega_l), 0, z)[0] * 3.086e24

S_f  = ELECTRON_ENERGY * LIGHT_SPEED / D_l**2 / (1 + z)
S_f *= R**2 * gamma**4
if   volume == "shell":  S_f *= 1 / np.pi    # Shell
elif volume == "sphere": S_f *= 4./ 9        # Sphere
else: assert(0)


nu_0 = 3 / (4 * np.pi) * B * q / (m*c)

plt.close()
plt.rcdefaults()
# fig = plt.figure()
plt.rcParams.update({'font.size': 20})
fig = plt.figure(figsize=(18,9))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax3 = ax1.twiny()

def tic_function(x):
    x2 = x / 4.13e-15
    return x2

ax1.set_ylim([1e-17, 1e-9])
ax1.set_ylim([1e-15, 1e-9])
ax1.grid()
ax1.set_xlabel("Energy (eV)")
ax1.set_ylabel("$\\nu F_{\\nu} (erg \\cdot s^{-1} \\cdot cm^{-2})$")

def lnlike(data):

    sim_data_x = data[0] * 511000
    sim_data_y = data[0] * data[0] * data[2]

    # Transform the simulated data from the fluid frame to the obs frame
    sim_data_x *= gamma / (1 + z)
    sim_data_y *= S_f

    radio_sim_data      = np.interp(radio_obs_data[0],      sim_data_x, sim_data_y)
    optical_sim_data    = np.interp(optical_obs_data[0],    sim_data_x, sim_data_y)
    x_rays_sim_data     = np.interp(x_rays_obs_data[0],     sim_data_x, sim_data_y)
    gamma_rays_sim_data = np.interp(gamma_rays_obs_data[0], sim_data_x, sim_data_y)

    optical_lnlike = -np.sum(0.5*np.power((optical_sim_data - optical_obs_data[1]) / (optical_obs_errors), 2)
                             + np.log(optical_obs_errors))
    x_rays_lnlike = -np.sum(0.5*np.power((x_rays_sim_data - x_rays_obs_data[1]) / (x_rays_obs_errors), 2)
                             + np.log(x_rays_obs_errors))

    # Gamma rays are trickier because the errors are asymmetric
    gamma_rays_obs_errors_temp = gamma_rays_obs_errors[0].copy()
    gamma_rays_obs_errors_temp[gamma_rays_sim_data > gamma_rays_obs_data[1]] = \
        gamma_rays_obs_errors[1][gamma_rays_sim_data > gamma_rays_obs_data[1]]

    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2)
                                + np.log(gamma_rays_obs_errors_temp))

    optical_lnlike    = -np.sum(0.5*np.power((optical_sim_data    - optical_obs_data[1])    / (optical_obs_errors), 2))
    x_rays_lnlike     = -np.sum(0.5*np.power((x_rays_sim_data     - x_rays_obs_data[1])     / (x_rays_obs_errors), 2))
    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2))

    return optical_lnlike + x_rays_lnlike + gamma_rays_lnlike,  \
            (optical_lnlike, x_rays_lnlike, gamma_rays_lnlike)


first_plot = True
while True:
    photon_data_filename = "./data/photon_data_{:04d}.tsv".format(i)
    neutrino_data_filename = "./data/neutrino_data_{:04d}.tsv".format(i)

    # photon_data_filename   = photon_data_file
    # neutrino_data_filename = neutrino_data_file

    try:
        data = np.loadtxt(photon_data_filename, unpack=True)
        data_0_eV = data[0] * 511000
        # ax1.loglog(data_0_eV * gamma / (1 + z), data[0]*data[0]*data[2] * S_f, color='gray', linestyle=':')
        ax1.loglog(data_0_eV * gamma / (1 + z), data[0]*data[0]*data[2] * S_f, color='b')

        neutrino_data = np.loadtxt(neutrino_data_filename, unpack=True)
        neutrino_data_0_eV = neutrino_data[0] * 511000
        ax1.loglog(neutrino_data_0_eV * gamma / (1 + z), neutrino_data_0_eV**2 *(neutrino_data[3]+neutrino_data[4]) * S_f, 'r')

        photon_gains = [1e-307] * len(data[0])
        photon_gains += np.sum(data[4:15], 0)

        aux = data[0] * data[0] * data[2] * S_f / photon_gains
        # ax1.loglog(data_0_eV * gamma / (1 + z), aux * data[4],  label="$e^- Synchrotron$")
        # ax1.loglog(data_0_eV * gamma / (1 + z), aux * data[5],  label="$p^+ Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[6],  label="$\pi^+ Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[7],  label="$\pi^- Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[8],  label="$\mu^+_L Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[9],  label="$\mu^+_R Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[10], label="$\mu^-_L Synchrotron$")
        # ax1.loglog(data_0_eV, aux * data[11], label="$\mu^-_R Synchrotron$")
        # ax1.loglog(data_0_eV * gamma / (1 + z), aux * data[12], label="$IC\ up$")
        # ax1.loglog(data_0_eV, aux * data[13], label="$IC\ down$")
        # ax1.loglog(data_0_eV * gamma / (1 + z), aux * data[14], label="$\pi^0\ decay$")

        # photon_losses = [1e-308] * len(data[0])
        # photon_losses += np.sum(data[15:], 0)

        # aux = data[0] * data[0] * data[2] * S_f / photon_losses
        # ax1.loglog(data_0_eV, aux * data[-2], label="$PP$")
        # ax1.loglog(data_0_eV, aux * data[-1], label="$Escape$")

        print(config_file)
        print(photon_data_filename)
        print(lnlike(data))

        if first_plot:
            ax1.lines[0].set_label("Simulated SED")
            ax1.lines[1].set_label("Neutrino Flux")

            ax2.set_xlabel("Frequency (Hz)")
            # ax2.xaxis.set_label_coords(0.5,1.05)
            ax2.xaxis.set_label_coords(0.5,1.075)
            ax2.set_xlim(tic_function(np.array(ax1.get_xbound())))
            ax2.loglog(tic_function(data_0_eV * gamma / (1 + z)), data[0]*data[0]*data[2] * S_f, color='b')

            ax3.set_xlabel("Normalized Energy")
            ax3.set_xlim(np.array(ax1.get_xbound()) / 511000)
            ax3.loglog(data[0] * gamma / (1 + z), data[0]*data[0]*data[2] * S_f, color='b')
            ax3.spines['top'].set_position(('outward', 60))

            # indices = data_0_eV < 0.4
            # ax1.loglog(data_0_eV[indices] * gamma / (1 + z),
                       # data[0][indices]*data[0][indices]*data[2][indices] * S_f, color='black', linewidth=3, label="Trusted range")

            # indices = np.logical_and(1e5 < data_0_eV, data_0_eV < 3e7)
            # ax1.loglog(data_0_eV[indices] * gamma / (1 + z),
                       # data[0][indices]*data[0][indices]*data[2][indices] * S_f, color='black', linewidth=3)


            first_plot = False

    except IOError as e:
        break

    break
    i += step

ax1.errorbar(     radio_obs_data[0],      radio_obs_data[1], yerr=0.5*radio_obs_data[1], uplims=[True], fmt=' m', capsize=5, label="Experimental Data")
ax1.errorbar(   optical_obs_data[0],    optical_obs_data[1], yerr=optical_obs_errors   ,                fmt='om')
ax1.errorbar(    x_rays_obs_data[0],     x_rays_obs_data[1], yerr=x_rays_obs_errors    ,                fmt='om')
ax1.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], yerr=gamma_rays_obs_errors,                fmt='om')

ax1.errorbar(2.5e14, 5e-11, yerr=2e-11, uplims=[True], capsize=5)
ax1.text    (2.5e14, 7e-11,
        """Flux Corresponding to\n1 $\\nu_{\\mu}$ in IceCube\nper 1/2 year""", ha="center")

x_extra = np.logspace(np.log10(data_0_eV[0]), np.log10(data_0_eV[-1]), 1024)
p_extra_1 = 1.2
y_extra_1 = x_extra**(p_extra_1) * radio_obs_data[1][0] / (radio_obs_data[0][0] / gamma * (1+z))**(p_extra_1)

p_extra_2 = -0.5
y_extra_2 = x_extra**(p_extra_2) * optical_obs_data[1][0] / (optical_obs_data[0][0] / gamma * (1+z))**(p_extra_2)

p_extra_3 = 0.35
y_extra_3 = x_extra**(p_extra_3) * x_rays_obs_data[1][10] / (x_rays_obs_data[0][10] / gamma * (1+z))**(p_extra_3)

# ax1.loglog(x_extra,y_extra_1)
# ax1.loglog(x_extra,y_extra_2)
# ax1.loglog(x_extra,y_extra_3)

# ax2.plot([g_min*g_min * nu_0, g_min*g_min * nu_0],[1e-17,1e-9])
# ax2.plot([g_c  *g_c   * nu_0, g_c  *g_c   * nu_0],[1e-17,1e-9])
# ax2.plot([g_max*g_max * nu_0, g_max*g_max * nu_0],[1e-17,1e-9])

# ax2.plot([(g_min*g_min)**2 * nu_0, (g_min*g_min)**2 * nu_0],[1e-17,1e-9])
# ax2.plot([(g_min*g_c)**2 * nu_0, (g_min*g_c)**2 * nu_0],[1e-17,1e-9])
# ax2.plot([(g_c  *g_c)**2 * nu_0, (g_c  *g_c)**2 * nu_0],[1e-17,1e-9])
# ax2.plot([(g_max*g_c)**2 * nu_0, (g_max*g_c)**2 * nu_0],[1e-17,1e-9])
# ax2.plot([(g_min*g_max)**2 * nu_0, (g_min*g_max)**2 * nu_0],[1e-17,1e-9])

# ax3.plot([g_min, g_min],[1e-17,1e-9])
# ax3.plot([g_c, g_c]    ,[1e-17,1e-9])
# ax3.plot([g_max, g_max],[1e-17,1e-9])

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.tight_layout(0)
ax2.set_xlim(tic_function(np.array(ax1.get_xbound())))
ax3.set_xlim(np.array(ax1.get_xbound()) / 511000)
ax1.legend(loc='upper left')
# plt.savefig("test.png")
plt.show()
