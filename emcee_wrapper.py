import emcee

import subprocess
import numpy as np
from io import StringIO
from string import Template
import os

import matplotlib.pyplot as plt

import corner

z = 0.3365
ELECTRON_MASS = 9.109e-28
LIGHT_SPEED   = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED

def obs_to_data(obs_data):
    a = (obs_data[0] - 150) / (592 - 150)
    b = (obs_data[1] - 113) / (725 - 113)

    c = 10**(-6  * (1 - a) + a * 3)
    d = 10**(-10 * (1 - b) + b * -14)

    return np.array([c,d])

# Inlined experimental data
radio_obs = np.array([
    [200, 708],
    [207, 683],
    [216, 652],
    [224, 622],
    [227, 606],
    [232, 591],
    [232, 586],
    [238, 573]
]).transpose()
radio_obs_data = obs_to_data(radio_obs)

optical_obs = np.array([
    [460, 185],
    [459, 185],
    [463, 170],
    [463, 177],
    [467, 188],
    [467, 179],
    [472, 177],
    [479, 185],
    [482, 177],
    [484, 182]
]).transpose()
optical_obs_data = obs_to_data(optical_obs)
# Assume that the errors in optical data are small
optical_obs_errors = 0.1 * optical_obs_data[1]

x_rays_obs = np.array([
    [573, 350],
    [582, 347],
    [588, 361],
    [593, 382],
    [597, 382],
    [601, 391],
    [604, 396],
    [610, 414],
    [618, 420],
    [625, 430],
    [632, 429],
    [636, 415],
    [643, 412],
    [650, 398],
    [663, 400],
    [679, 380]
]).transpose()
x_rays_obs_data = obs_to_data(x_rays_obs)
# Assume that the errors in x-ray data are small
x_rays_obs_errors = 0.2 * x_rays_obs_data[1]

gamma_rays_obs = np.array([
    [ 851, 148],
    [ 854, 121],
    [ 875, 152],
    [ 900, 153],
    [ 924, 150],
    [ 948, 164],
    [ 973, 200],
    [ 981, 286],
    [ 998, 334],
    [1015, 399]
]).transpose()
gamma_rays_obs_data = obs_to_data(gamma_rays_obs)
# Assume that the upper errors in gamma ray data are small
# but the lower errors are awful
gamma_rays_obs_errors = [0.2 * gamma_rays_obs_data[1], 0.5 * gamma_rays_obs_data[1]]

config_template = Template("""
    [general]
    density = $dens
    magnetic_field = $mag
    dt = 0.1
    t_max = 1e7
    R = $radius

    [electrons]
    gamma_min = $g_min
    gamma_max = $g_max
    size      = 160
    distribution_type = "broken_power_law"
    break_point = $g_break
    first_slope  =  $p1
    second_slope =  $p2

    [protons]
    size      = 160
    distribution_type = "power_law"
    slope  =  2

    [photons]
    size = 128
    """)

# def theta_to_params(theta):
    # log_gamma_min, log_gamma_break, log_gamma_max,  \
    # first_slope, second_slope,                      \
    # log_R, B, log_S_f, Gamma, log_density = theta

    # return (10**log_gamma_min, 10**log_gamma_break, 10**log_gamma_max,
            # first_slope, second_slope,
            # 10**log_R, B, 10**log_S_f, Gamma, 10**log_density)

# def params_to_theta(params):
    # gamma_min, gamma_break, gamma_max,  \
    # first_slope, second_slope,          \
    # R, B, S_f, Gamma, density = params

    # return (np.log10(gamma_min), np.log10(gamma_break), np.log10(gamma_max),
            # first_slope, second_slope,
            # np.log10(R), B, np.log10(S_f), Gamma, np.log10(density))

def theta_to_params(theta):
    log_gamma_min, log_gamma_break, log_gamma_max,  \
    first_slope, second_slope,                      \
    log_R, B, log_S_f, Gamma, log_density_R = theta

    return (10**log_gamma_min, 10**log_gamma_break, 10**log_gamma_max,
            first_slope, second_slope,
            10**log_R, B, 10**log_S_f, Gamma, 10**(log_density_R - log_R))

def params_to_theta(params):
    gamma_min, gamma_break, gamma_max,  \
    first_slope, second_slope,          \
    R, B, S_f, Gamma, density = params

    return (np.log10(gamma_min), np.log10(gamma_break), np.log10(gamma_max),
            first_slope, second_slope,
            np.log10(R), B, np.log10(S_f), Gamma, np.log10(density * R))

def write_config_file(theta, config_filename):
    params = theta_to_params(theta)
    config = config_template.substitute(
                g_min = params[0],
                g_break = params[1],
                g_max = params[2],
                p1 = params[3],
                p2 = params[4],
                radius = params[5],
                mag = params[6],
                dens = params[9])

    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

def generate_sim_data(theta):
    config_filename = "config_emcee_{}.toml".format(os.getpid())

    write_config_file(theta, config_filename)
    process = subprocess.Popen(['./model',config_filename], stdout=subprocess.PIPE, universal_newlines=True)
    output_data, err = process.communicate()
    data = np.loadtxt(StringIO(output_data), unpack=True)
    os.remove(config_filename)

    return data

def lnlike(theta):
    gamma_min, gamma_break, gamma_max,  \
    first_slope, second_slope,          \
    R, B, S_f, Gamma, density = theta_to_params(theta)

    sim_data = generate_sim_data(theta)

    # sim_data[0] in epsilon -> eV
    # sim_data[1] in pop     -> erg
    sim_data[1] *= sim_data[0] * sim_data[0] * ELECTRON_ENERGY
    sim_data[0] *= 511000

    # Transform the simulated data from the fluid frame to the obs frame
    sim_data[0] *= Gamma / (1 + z)
    sim_data[1] *= S_f

    radio_sim_data      = np.interp(radio_obs_data[0],      sim_data[0], sim_data[1])
    optical_sim_data    = np.interp(optical_obs_data[0],    sim_data[0], sim_data[1])
    x_rays_sim_data     = np.interp(x_rays_obs_data[0],     sim_data[0], sim_data[1])
    gamma_rays_sim_data = np.interp(gamma_rays_obs_data[0], sim_data[0], sim_data[1])

    # The radio band is an upper limit
    if np.any(radio_sim_data > radio_obs_data): return -np.inf, sim_data

    optical_lnlike = -np.sum(0.5*np.power((optical_sim_data - optical_obs_data[1]) / (optical_obs_errors), 2)
                             + np.log(optical_obs_errors))
    x_rays_lnlike = -np.sum(0.5*np.power((x_rays_sim_data - x_rays_obs_data[1]) / (x_rays_obs_errors), 2)
                             + np.log(x_rays_obs_errors))

    # Gamma rays are trickier because the errors are asymmetric
    gamma_rays_obs_errors_temp = gamma_rays_obs_errors[0].copy()
    gamma_rays_obs_errors_temp[gamma_rays_sim_data < gamma_rays_obs_data[1]] = \
        gamma_rays_obs_errors[1][gamma_rays_sim_data < gamma_rays_obs_data[1]]

    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2)
                                + np.log(gamma_rays_obs_errors_temp))

    optical_lnlike = -np.sum(0.5*np.power((optical_sim_data - optical_obs_data[1]) / (optical_obs_errors), 2))
    x_rays_lnlike = -np.sum(0.5*np.power((x_rays_sim_data - x_rays_obs_data[1]) / (x_rays_obs_errors), 2))
    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2))

    return optical_lnlike + x_rays_lnlike + gamma_rays_lnlike, sim_data

def lnprior(theta):
    gamma_min, gamma_break, gamma_max,  \
    first_slope, second_slope,          \
    R, B, S_f, Gamma, density = theta_to_params(theta)

    if 10        < gamma_min    < gamma_max and \
       gamma_min < gamma_break  < gamma_max and \
       gamma_min < gamma_max    < 1e7 and       \
       1         < first_slope  < 5 and         \
       1         < second_slope < 5 and         \
                                                \
       1e13  <  R       < 1e18 and  \
       0.1   <  B       < 0.3 and   \
       1e-10 <  S_f     < 1e-5 and  \
       5     <  Gamma   < 20 and    \
       1e-24 <  density < 1e-19:
           return 0
    else:
        return -np.inf

def lnprob(theta):
    lp = lnprior(theta)
    if np.isinf(lp):
        return -np.inf, np.array([[0]*128,[0]*128])
    else:
        lnlike_data = lnlike(theta)
        return lp + lnlike_data[0], lnlike_data[1]

ndim = 10
walkers = 128

# limits = np.array([(1, 4),      # log10 g_min
                   # (2, 6),      # log10 g_break
                   # (3, 7),      # log10 g_max
                   # (1, 4),      # p1
                   # (2, 5),      # p2
                   # (13, 18),    # log10 R
                   # (0.1, 0.3),  # B
                   # (-10, -5),   # log10 S_f
                   # (5, 20),     # Gamma
                   # (-24, -19)]) # log10 d

# labels=["$\\log_{10}\\left(\\gamma_{min}\\right)$", \
        # "$\\log_{10}\\left(\\gamma_{break}\\right)$", \
        # "$\\log_{10}\\left(\\gamma_{max}\\right)$", \
        # "$p_1$", \
        # "$p_2$", \
        # "$\\log_{10}\\left(R\\right)$", \
        # "$B$", \
        # "$\\log_{10}\\left(S_f\\right)$", \
        # "$\\Gamma$", \
        # "$\\log_{10}\\left(\\rho\\right)$"]

# real_params = [2.5, 4.2,  5.5, \
               # 2.17, 4.15, \
               # np.log10(6.3e15), 0.17, np.log10(3.5e-7), 10, np.log10(4e-22)]

# Note that this are the limits for the corner plot!
# They are different from the ones in our priors
limits = np.array([(1, 4),      # log10 g_min
                   (3, 6),      # log10 g_break
                   (4, 7),      # log10 g_max
                   (1, 3),      # p1
                   (3, 5),      # p2
                   (14, 17),    # log10 R
                   (0.1, 0.25), # B
                   (-8, -5),    # log10 S_f
                   (5, 15),     # Gamma
                   (-7, -4)])   # log10 d R

labels=["$\\log_{10}\\left(\\gamma_{min}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{break}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{max}\\right)$", \
        "$p_1$", \
        "$p_2$", \
        "$\\log_{10}\\left(R\\right)$", \
        "$B$", \
        "$\\log_{10}\\left(S_f\\right)$", \
        "$\\Gamma$", \
        "$\\log_{10}\\left(\\rho \cdot R\\right)$"]

real_params = [10**2.5, \
               10**4.2, \
               10**5.5, \
               2.17,    \
               4.15,    \
               6.3e15,  \
               0.17,    \
               3.5e-7,  \
               10,      \
               4e-22]

real_theta = params_to_theta(real_params)
print(real_theta)

input_parameters_lnprob, input_parameters_data = lnprob(real_theta)

plt.figure(figsize=(24,12), dpi=128)
plt.ylim([1e-17, 1e-9])
plt.xlim([1e-6, 1e12])
plt.grid()
plt.loglog(input_parameters_data[0], input_parameters_data[1])
plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
# plt.loglog(optical_obs_data[0], optical_obs_data[1], 'o')
# plt.loglog(x_rays_obs_data[0], x_rays_obs_data[1], 'o')
# plt.loglog(gamma_rays_obs_data[0], gamma_rays_obs_data[1], 'o')
plt.errorbar(optical_obs_data[0], optical_obs_data[1], fmt='o', yerr=optical_obs_errors)
plt.errorbar(x_rays_obs_data[0], x_rays_obs_data[1], fmt='o', yerr=x_rays_obs_errors)
plt.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], fmt='o', yerr=gamma_rays_obs_errors)

# plt.legend(loc='best')
plt.savefig("input_parameters.png", dpi='figure')
plt.clf()
plt.close()
print(input_parameters_lnprob)

# pos = np.array([real_theta * (1 + 0.1*np.random.normal(size=len(real_theta)) * 4 / np.array(real_theta)) for i in range(walkers)])
pos = np.array([np.random.normal(loc=real_theta, scale=(limits[:,1] - limits[:,0]) / 8) for i in range(walkers)])
print(pos[0])

lnprobs = None
blobs = None
sampler = emcee.EnsembleSampler(walkers, ndim, lnprob, threads=4)

# plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
# plt.loglog(optical_obs_data[0], optical_obs_data[1], 'o')
# plt.loglog(x_rays_obs_data[0], x_rays_obs_data[1], 'o')
# plt.loglog(gamma_rays_obs_data[0], gamma_rays_obs_data[1], 'o')
# plt.show()

total_iterations = 100 * 1

# print("Percent\tlog10(g_min)\tlog10(g_break)\tlog10(g_max)\tp_1\t\tp_2\t\tlog10(R)\tB\t\tlog10(S_f)\tGamma\t\tlog10(rho)")
print("Percent\tlog10(g_min)\tlog10(g_break)\tlog10(g_max)\tp_1\t\tp_2\t\tlog10(R)\tB\t\tlog10(S_f)\tGamma\t\tlog10(rho.R)")
for i in range(100):
    pos, lnprobs, ignored, blobs = sampler.run_mcmc(pos, total_iterations / 100, lnprob0=lnprobs, blobs0=blobs)
    samples = sampler.chain[:,:,:].reshape((-1,ndim))

    blobs = np.array(blobs)

    lnprobs_sort_indices = lnprobs.argsort()
    sorted_lnprobs = lnprobs[lnprobs_sort_indices[::-1]]
    sorted_pos     = pos[lnprobs_sort_indices[::-1]]
    sorted_blobs   = blobs[lnprobs_sort_indices[::-1]]

    # log_g_min_mcmc, log_g_break_mcmc, log_g_max_mcmc, \
    # p1_mcmc, p2_mcmc, \
    # log_R_mcmc, B_mcmc, log_S_f_mcmc, Gamma_mcmc, log_density_mcmc = \
            # map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
                    # zip(*np.percentile(samples, [16,50,84], axis=0)))

    log_g_min_mcmc, log_g_break_mcmc, log_g_max_mcmc, \
    p1_mcmc, p2_mcmc, \
    log_R_mcmc, B_mcmc, log_S_f_mcmc, Gamma_mcmc, log_density_R_mcmc = \
            map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
                    zip(*np.percentile(samples, [16,50,84], axis=0)))

    print("{:d}%".format(i+1), end='\t')
    print("{:.3f} ± {:.3f}".format(log_g_min_mcmc[0],    (log_g_min_mcmc[1]    + log_g_min_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(log_g_break_mcmc[0],    (log_g_break_mcmc[1]    + log_g_break_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(log_g_max_mcmc[0],    (log_g_max_mcmc[1]    + log_g_max_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(p1_mcmc[0],    (p1_mcmc[1]    + p1_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(p2_mcmc[0],    (p2_mcmc[1]    + p2_mcmc[2])    / 2), end='\t')

    print("{:.3f} ± {:.3f}".format(log_R_mcmc[0],    (log_R_mcmc[1]    + log_R_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(B_mcmc[0],    (B_mcmc[1]    + B_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(log_S_f_mcmc[0],    (log_S_f_mcmc[1]    + log_S_f_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(Gamma_mcmc[0],    (Gamma_mcmc[1]    + Gamma_mcmc[2])    / 2), end='\t')
    print("{:.3f} ± {:.3f}".format(log_density_R_mcmc[0],    (log_density_R_mcmc[1]    + log_density_R_mcmc[2])    / 2), end='\t')
    # print("{:.3f} ± {:.3f}".format(log_density_mcmc[0],    (log_density_mcmc[1]    + log_density_mcmc[2])    / 2), end='\t')

    print("", end='\n')
    print("BEST:", end='\t')
    print("10^{:.3f}".format(sorted_pos[0][0]), end='\t')
    print("10^{:.3f}".format(sorted_pos[0][1]), end='\t')
    print("10^{:.3f}".format(sorted_pos[0][2]), end='\t')
    print("{:.3f}".format(sorted_pos[0][3]), end='\t\t')
    print("{:.3f}".format(sorted_pos[0][4]), end='\t\t')
    print("10^{:.3f}".format(sorted_pos[0][5]), end='\t')
    print("{:.3f}".format(sorted_pos[0][6]), end='\t\t')
    print("10^{:.3f}".format(sorted_pos[0][7]), end='\t')
    print("{:.3f}".format(sorted_pos[0][8]), end='\t\t')
    print("10^{:.3f}".format(sorted_pos[0][9]), end='\t')
    print("{:.3f}".format(sorted_lnprobs[0]))
    print("", end='\n')

    plt.figure(figsize=(24,12), dpi=128)
    plt.ylim([1e-17, 1e-9])
    plt.xlim([1e-6, 1e12])
    plt.grid()
    plt.loglog(sorted_blobs[0][0], sorted_blobs[0][1])
    plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
    # plt.loglog(optical_obs_data[0], optical_obs_data[1], 'o')
    # plt.loglog(x_rays_obs_data[0], x_rays_obs_data[1], 'o')
    # plt.loglog(gamma_rays_obs_data[0], gamma_rays_obs_data[1], 'o')
    plt.errorbar(optical_obs_data[0], optical_obs_data[1], fmt='o', yerr=optical_obs_errors)
    plt.errorbar(x_rays_obs_data[0], x_rays_obs_data[1], fmt='o', yerr=x_rays_obs_errors)
    plt.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], fmt='o', yerr=gamma_rays_obs_errors)

    # plt.legend(loc='best')
    plt.savefig("emcee_best_at_step_{:03d}.png".format(i), dpi='figure')
    plt.clf()
    plt.close()

    fig = corner.corner(samples, labels=labels, range=limits, truths=real_theta, quantiles=[0.16,0.5,0.84], show_titles=True)
    fig.savefig("triangle_at_step_{:03d}.png".format(i))
    plt.close(fig)

samples = sampler.chain[:,:,:].reshape((-1,ndim))
np.savetxt("data_mcmc.gz", samples)
fig = corner.corner(samples, labels=labels, range=limits, truths=real_theta, quantiles=[0.16,0.5,0.84], show_titles=True)
fig.savefig("triangle.png")
