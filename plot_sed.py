import matplotlib.pyplot as plt
import numpy as np
plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax3 = ax1.twiny()

B = 0.17
gamma  = 10
z      = 1.34
factor = 2e-7
step   = 1

g_min = 10**2.5
g_c   = 10**4.2
g_max = 10**5.5
q = 4.8032045e-10
m = 9.109e-28
c = 3e10

nu_0 = 3 / (4 * np.pi) * B *q / (m*c)

def tic_function(x):
    x2 = x / 4.13e-15
    return x2

ax1.set_ylim([1e-17, 1e-9])
ax1.grid()

first_plot = True
i = 0
while i<17:
    # filename = ""
    # label = ""
    # if i == 0:
        # filename = "./data/data_linn_logg_implicit_NO_limit/photon_data_{:04d}.tsv".format(16)
        # label = "NO limit"
    # elif i == 1:
        # filename = "./data/data_linn_logg_implicit_WITH_limit/photon_data_{:04d}.tsv".format(16)
        # label = "WITH limit"
    # elif i == 2:
        # filename = "./data/data_logn_logg_2_order/photon_data_{:04d}.tsv".format(16)
        # label = "Logn"
    # else:
        # i+=1
        # continue


    # filename = "./data/data_linn_logg_implicit_NO_limit/photon_data_{:04d}.tsv".format(i)
    # filename = "./data/data_linn_logg_implicit_WITH_limit/photon_data_{:04d}.tsv".format(i)
    filename = "./data/data_logn_logg_2_order/photon_data_{:04d}.tsv".format(i)

    try:
        print(filename)
        data = np.loadtxt(filename,unpack=True)
        data_0_eV = data[0] * 511000
        ax1.loglog(data_0_eV, data[0]*data[1]*data[2] * factor)

        if first_plot:
            ax2.set_xlim(tic_function(np.array(ax1.get_xlim())))
            ax2.loglog(tic_function(data_0_eV), data[0]*data[1]*data[2] * factor, 'r')

            ax3.set_xlim(np.array(ax1.get_xlim()) / 511000)
            ax3.loglog(data[0], data[0]*data[1]*data[2] * factor, 'g')
            ax3.spines['top'].set_position(('outward', 30))

            first_plot = False
    except Exception as e:
        break

    i += step

obs_data = np.loadtxt("data.dat", unpack=True)
obs_data[0] = (obs_data[0] - 150) / (592 - 150)
obs_data[1] = (obs_data[1] - 113) / (725 - 113)
obs_data_x = 10**(-6  * (1 - obs_data[0]) + obs_data[0] * 3)
obs_data_y = 10**(-10 * (1 - obs_data[1]) + obs_data[1] * -14)

ax1.loglog(obs_data_x / gamma * z, obs_data_y, 'o')
ax1.errorbar(obs_data_x / gamma * z, obs_data_y, yerr=0.5*obs_data_y, fmt='o')

x_extra = np.logspace(-7, 12, 1024)
p_extra_1 = 1.2
y_extra_1 = x_extra**(p_extra_1) * obs_data_y[0] / (obs_data_x[0] / gamma * z)**(p_extra_1)

p_extra_2 = -0.5
y_extra_2 = x_extra**(p_extra_2) * obs_data_y[27] / (obs_data_x[27] / gamma * z)**(p_extra_2)

p_extra_3 = 0.35
y_extra_3 = x_extra**(p_extra_3) * obs_data_y[28] / (obs_data_x[28] / gamma * z)**(p_extra_3)

ax1.loglog(x_extra,y_extra_1)
ax1.loglog(x_extra,y_extra_2)
ax1.loglog(x_extra,y_extra_3)

ax2.plot([g_min*g_min * nu_0, g_min*g_min * nu_0],[1e-17,1e-9])
ax2.plot([g_c  *g_c   * nu_0, g_c  *g_c   * nu_0],[1e-17,1e-9])
ax2.plot([g_max*g_max * nu_0, g_max*g_max * nu_0],[1e-17,1e-9])

ax2.plot([g_c**4 * nu_0, g_c**4 * nu_0],[1e-17,1e-9])
ax2.plot([(g_min*g_max)**2 * nu_0, (g_min*g_max)**2 * nu_0],[1e-17,1e-9])

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
# figManager.window.state('zoomed')
# figManager.full_screen_toggle()
# fig.tight_layout(0)
plt.tight_layout(0)
plt.show()
