import matplotlib.pyplot as plt
import numpy as np
plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)

gamma  = 10
z      = 1.34
factor = 2e-7
step   = 1

q = 4.8032045e-10
m = 9.109e-28
c = 3e10
h_ev = 4.135667662e-15

indices = []
data_times      = []
data_radio      = []
data_optical    = []
data_x_rays     = []
data_gamma_rays = []
first_data = True
i = 0
while True:
    filename = "data/photon_data_{:04d}.tsv".format(i)

    try:
        data = np.loadtxt(filename,unpack=True)

        if first_data:
            data_0_eV = data[0]   * 511000

            indices.append((np.abs(data_0_eV - 10**-4.5)).argmin()) # Radio
            indices.append((np.abs(data_0_eV - 10**0   )).argmin()) # Optical
            indices.append((np.abs(data_0_eV - 10**3   )).argmin()) # X-Rays
            indices.append((np.abs(data_0_eV - 10**9   )).argmin()) # Gamma-Rays

            indices = np.array(indices)

            first_data = False

        with open(filename) as f:
            time = float(f.readline().strip().split()[1])
            data_times.append(time)
        data_radio     .append(data[0][indices[0]] * data[1][indices[0]] * data[2][indices[0]])
        data_optical   .append(data[0][indices[1]] * data[1][indices[1]] * data[2][indices[1]])
        data_x_rays    .append(data[0][indices[2]] * data[1][indices[2]] * data[2][indices[2]])
        data_gamma_rays.append(data[0][indices[3]] * data[1][indices[3]] * data[2][indices[3]])
    except Exception as e:
        break

    i += step

data_radio = np.array(data_radio)

ax1.loglog(data_times, data_radio,      label='radio')
ax1.loglog(data_times, data_optical,    label='optical')
ax1.loglog(data_times, data_x_rays,     label='x_rays')
ax1.loglog(data_times, data_gamma_rays, label='gamma_rays')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
# figManager.window.state('zoomed')
# figManager.full_screen_toggle()
# fig.tight_layout(0)
plt.tight_layout(0)
plt.legend()
plt.show()
