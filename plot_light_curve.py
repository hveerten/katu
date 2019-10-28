import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integ
plt.close()

plt.rcdefaults()
plt.rcParams.update({'font.size': 20})
fig = plt.figure(figsize=(24,12), dpi=128)
ax1 = fig.add_subplot(111)
ax1.grid()

gamma  = 10
z      = 0.3365
factor = 2e-7
step   = 1

q = 4.8032045e-10
m = 9.109e-28
c = 3e10
h_ev = 4.135667662e-15

H_0 = 67.66
Omega_m = 0.3111
Omega_l = 0.6889

D_l = (1 + z) * c / (H_0 * 1e5) * integ.quad(lambda x: 1 / np.sqrt((1 + x)**3 * Omega_m + Omega_l), 0, z)[0] * 3.086e24

R = 9.62949420e+18
gamma = 5.34897856
S_f = m * c * c * c / np.pi * (R * gamma / D_l)**2  # Shell

factor = S_f

indices  = []
energies = []
data_times      = []
data_radio      = []
data_optical    = []
data_x_rays     = []
data_gamma_rays = []
data_neutrinos  = []
first_data = True
i = 0
while True:
    filename_photon   = "data/photon_data_{:04d}.tsv".format(i)
    filename_neutrino = "data/neutrino_data_{:04d}.tsv".format(i)

    try:
        data_photon   = np.loadtxt(filename_photon,unpack=True)
        data_neutrino = np.loadtxt(filename_neutrino,unpack=True)

        if first_data:
            data_0_eV = data_photon[0] * 511000

            indices.append((np.abs(data_0_eV - 10**-5)).argmin()) # Radio
            indices.append((np.abs(data_0_eV - 10**0 )).argmin()) # Optical
            indices.append((np.abs(data_0_eV - 10**3 )).argmin()) # X-Rays
            indices.append((np.abs(data_0_eV - 10**9 )).argmin()) # Gamma-Rays
            energies.append(data_photon[0][indices[0]])
            energies.append(data_photon[0][indices[1]])
            energies.append(data_photon[0][indices[2]])
            energies.append(data_photon[0][indices[3]])

            data_0_eV = data_neutrino[0] * 511000
            indices.append((np.abs(data_0_eV - 10**14)).argmin()) # neutrinos
            energies.append(data_neutrino[0][indices[4]])

            indices = np.array(indices)

            first_data = False

        with open(filename_photon) as f:
            time = float(f.readline().strip().split()[1])
            data_times.append(time)

        data_radio     .append(energies[0] * energies[0] * data_photon[2][indices[0]])
        data_optical   .append(energies[1] * energies[1] * data_photon[2][indices[1]])
        data_x_rays    .append(energies[2] * energies[2] * data_photon[2][indices[2]])
        data_gamma_rays.append(energies[3] * energies[3] * data_photon[2][indices[3]])
        data_neutrinos .append(energies[4]**4            * (data_neutrino[3][indices[4]] + data_neutrino[4][indices[4]]))
    except Exception as e:
        break

    i += step

ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Flux (AU)")

data_radio = np.array(data_radio)
data_optical = np.array(data_optical)
data_x_rays = np.array(data_x_rays)
data_gamma_rays = np.array(data_gamma_rays)
data_neutrinos = np.array(data_neutrinos)

ax1.semilogy(data_times, data_radio,      label='radio')
ax1.semilogy(data_times, data_optical,    label='optical')
ax1.semilogy(data_times, data_x_rays,     label='x_rays')
ax1.semilogy(data_times, data_gamma_rays, label='gamma_rays')
ax1.semilogy(data_times, data_neutrinos,  label='neutrinos')

t1 = 180000
t2 = t1 + 45000
t3 = t2 + 225050
t4 = t3 + 45000

ax1.plot([t1, t1], [1e-3, 1e3])
ax1.plot([t2, t2], [1e-3, 1e3])
ax1.plot([t3, t3], [1e-3, 1e3])
ax1.plot([t4, t4], [1e-3, 1e3])

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
# figManager.window.state('zoomed')
# figManager.full_screen_toggle()
# fig.tight_layout(0)
plt.tight_layout(1)
plt.legend()
# plt.show()
plt.savefig("./data/light_curve.png", dpi='figure')

# fig = plt.figure()
plt.clf()
ax1 = fig.add_subplot(111)
ax1.grid()
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Flux Normalized")

data_optical *= np.min(data_radio) / np.min(data_optical)
data_x_rays *= np.min(data_radio) / np.min(data_x_rays)
data_gamma_rays *= np.min(data_radio) / np.min(data_gamma_rays)
data_neutrinos *= np.min(data_radio) / np.min(data_neutrinos)

data_optical = (data_optical - np.min(data_optical)) / (np.max(data_optical) - np.min(data_optical)) * (np.max(data_radio) - np.min(data_radio)) + np.min(data_radio)
data_x_rays = (data_x_rays - np.min(data_x_rays)) / (np.max(data_x_rays) - np.min(data_x_rays)) * (np.max(data_radio) - np.min(data_radio)) + np.min(data_radio)
data_gamma_rays = (data_gamma_rays - np.min(data_gamma_rays)) / (np.max(data_gamma_rays) - np.min(data_gamma_rays)) * (np.max(data_radio) - np.min(data_radio)) + np.min(data_radio)
data_neutrinos = (data_neutrinos - np.min(data_neutrinos)) / (np.max(data_neutrinos) - np.min(data_neutrinos)) * (np.max(data_radio) - np.min(data_radio)) + np.min(data_radio)

ax1.semilogy(data_times, data_radio,      label='radio')
ax1.semilogy(data_times, data_optical,    label='optical')
ax1.semilogy(data_times, data_x_rays,     label='x_rays')
ax1.semilogy(data_times, data_gamma_rays, label='gamma_rays')
ax1.semilogy(data_times, data_neutrinos,  label='neutrinos')

ax1.plot([t1, t1], [3e-3, 6e-2])
ax1.plot([t2, t2], [3e-3, 6e-2])
ax1.plot([t3, t3], [3e-3, 6e-2])
ax1.plot([t4, t4], [3e-3, 6e-2])
ax1.set_ylim([3e-3,6e-2])

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
# figManager.window.state('zoomed')
# figManager.full_screen_toggle()
# fig.tight_layout(0)
plt.tight_layout(1)
plt.legend()
# plt.show()
plt.savefig("./data/light_curve_normalized.png", dpi='figure')
