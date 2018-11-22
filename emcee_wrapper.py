import emcee

import subprocess
import numpy as np
import scipy.integrate as integ
from io import StringIO
from string import Template

import os
import sys

import matplotlib.pyplot as plt

import corner

z = 0.3365
ELECTRON_MASS = 9.109e-28
LIGHT_SPEED   = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED

H_0 = 67.74
Omega_m = 0.3089
Omega_l = 0.6911

D_l = (1 + z) * LIGHT_SPEED / H_0 * integ.quad(lambda x: 1 / np.sqrt((1 + x)**3 * Omega_m + Omega_l), 0, z)[0] * 3.086e24

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
# Assume that the upper errors in gamma ray data are small (0.2), but
# we add also a high of an error due to possible absorption by the EGL (0.5)
# but the lower errors are awful
gamma_rays_obs_errors = [0.5 * gamma_rays_obs_data[1], (0.2 + 0.5) * gamma_rays_obs_data[1]]

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

def theta_to_params(theta):
    log_gamma_min, log_gamma_break, log_gamma_max,  \
    first_slope, second_slope,                      \
    log_Sf_R2, log_h, B, Gamma, log_density_h_g_min = theta

    return (10**log_gamma_min, 10**log_gamma_break, 10**log_gamma_max,
            first_slope, second_slope,
            10**log_Sf_R2, 10**log_h, B, Gamma, 10**(log_density_h_g_min - log_h - 2/3*log_gamma_min))

def params_to_theta(params):
    gamma_min, gamma_break, gamma_max,  \
    first_slope, second_slope,          \
    Sf_R2, h, B, Gamma, density = params

    return (np.log10(gamma_min), np.log10(gamma_break), np.log10(gamma_max),
            first_slope, second_slope,
            np.log10(Sf_R2), np.log10(h), B, Gamma, np.log10(density * h * np.power(gamma_min,2/3)))

def write_config_file(theta, config_filename):
    params = theta_to_params(theta)
    config = config_template.substitute(
                g_min = params[0],
                g_break = params[1],
                g_max = params[2],
                p1 = params[3],
                p2 = params[4],
                radius = params[6], # not radius, h
                mag = params[7],
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
    Sf_R2, h, B, Gamma, density = theta_to_params(theta)

    # S_f = ELECTRON_ENERGY / 2     * (R * Gamma * LIGHT_SPEED / D_l)**2  # Spherical
    S_f = Sf_R2 * (Gamma * LIGHT_SPEED / D_l)**2  # Shell

    sim_data = generate_sim_data(theta)

    # sim_data[0] in epsilon -> eV
    # sim_data[1] in pop     -> erg
    sim_data[1] *= sim_data[0] * sim_data[0]
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
    gamma_rays_obs_errors_temp[gamma_rays_sim_data > gamma_rays_obs_data[1]] = \
        gamma_rays_obs_errors[1][gamma_rays_sim_data > gamma_rays_obs_data[1]]

    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2)
                                + np.log(gamma_rays_obs_errors_temp))

    optical_lnlike    = -np.sum(0.5*np.power((optical_sim_data    - optical_obs_data[1])    / (optical_obs_errors), 2))
    x_rays_lnlike     = -np.sum(0.5*np.power((x_rays_sim_data     - x_rays_obs_data[1])     / (x_rays_obs_errors), 2))
    gamma_rays_lnlike = -np.sum(0.5*np.power((gamma_rays_sim_data - gamma_rays_obs_data[1]) / (gamma_rays_obs_errors_temp), 2))

    return optical_lnlike + x_rays_lnlike + gamma_rays_lnlike, sim_data

def lnprior(theta):
    gamma_min, gamma_break, gamma_max,  \
    first_slope, second_slope,          \
    Sf_R2, h, B, Gamma, density = theta_to_params(theta)

    if 1 < gamma_min    < gamma_break  < gamma_max < 1e8 and \
       1 < first_slope  < second_slope < 7 and               \
                                                              \
       1e25  <  Sf_R2   < 1e35 and  \
       1e12  <  h       < 1e18 and  \
       0.05  <  B       < 0.3 and   \
       5     <  Gamma   < 30 and    \
       1e-24 <  density < 1e-18:
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

# Note that this are the limits for the corner plot!
# They are different from the ones in our priors
# BUT! they are used to generate the initial guesses
limits = np.array([(1, 7),          # log10 g_min
                   (1, 7),          # log10 g_break
                   (1, 7),          # log10 g_max
                   (1, 5),          # p1
                   (1, 5),          # p2
                   (25, 35),        # log10 Sf·R^2
                   (13, 16),        # log10 h
                   (0.08, 0.25),    # B
                   (5, 20),         # Gamma
                   (-10, +2)])      # log10 (d·h·g_min^2/3)

labels=["$\\log_{10}\\left(\\gamma_{min}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{break}\\right)$", \
        "$\\log_{10}\\left(\\gamma_{max}\\right)$", \
        "$p_1$", \
        "$p_2$", \
        "$\\log_{10}\\left(S_f \cdot R^2\\right)$", \
        "$\\log_{10}\\left(h\\right)$", \
        "$B$", \
        "$\\Gamma$", \
        "$\\log_{10}\\left(\\rho \cdot h \cdot \\gamma_{min}^{2/3}\\right)$"]

real_params = [10**2.22, \
               10**4.14, \
               10**5.87, \
               1.85,    \
               4.32,    \
               10**29.68, \
               10**14.94,  \
               0.17,    \
               11.51,      \
               10**(-5.39 - 14.94)]

real_theta = params_to_theta(real_params)
print(real_theta)

input_parameters_lnprob, input_parameters_data = lnprob(real_theta)

plt.figure(figsize=(24,12), dpi=128)
plt.ylim([1e-17, 1e-9])
plt.xlim([1e-6, 1e12])
plt.grid()
plt.loglog(input_parameters_data[0], input_parameters_data[1])
plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
plt.errorbar(optical_obs_data[0], optical_obs_data[1], fmt='o', yerr=optical_obs_errors)
plt.errorbar(x_rays_obs_data[0], x_rays_obs_data[1], fmt='o', yerr=x_rays_obs_errors)
plt.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], fmt='o', yerr=gamma_rays_obs_errors)

plt.savefig("input_parameters.png", dpi='figure')
plt.clf()
plt.close()
print(input_parameters_lnprob)

pos = np.array([np.random.normal(loc=real_theta, scale=(limits[:,1] - limits[:,0]) / 4) for i in range(walkers)])

print(pos[0])
pos[0] = real_theta

lnprobs = None
blobs = None
rstate = None
sampler = emcee.EnsembleSampler(walkers, ndim, lnprob, threads=4)

try:
    def print_medians(i, samples, ndim):
        print("{:03d}".format(i+1), end='\t')
        for j in range(ndim):
            medians = np.percentile(samples[:,j] , [16,50,84])
            medians[0] = medians[1] - medians[0]
            medians[2] = medians[2] - medians[1]
            print("{:.3f} ± {:.3f}".format(medians[1], (medians[0] + medians[2]) / 2), end='\t')
        print("", end='\n')

    def print_best(i, best, best_lnprob):
        print("BEST:", end='\t')
        print("10^{:.3f}".format(best[0]), end='\t')
        print("10^{:.3f}".format(best[1]), end='\t')
        print("10^{:.3f}".format(best[2]), end='\t')
        print("{:.3f}".format(best[3]), end='\t\t')
        print("{:.3f}".format(best[4]), end='\t\t')
        print("10^{:.3f}".format(best[5]), end='\t')
        print("10^{:.3f}".format(best[6]), end='\t')
        print("{:.3f}".format(best[7]), end='\t\t')
        print("{:.3f}".format(best[8]), end='\t\t')
        print("10^{:.3f}".format(best[9]), end='\t')
        print("{:.3f}".format(best_lnprob))
        print("", end='\n')

    def create_fig(i, best_blob):
        plt.figure(figsize=(24,12), dpi=128)
        plt.ylim([1e-17, 1e-9])
        plt.xlim([1e-6, 1e12])
        plt.grid()
        plt.loglog(best_blob[0], best_blob[1])
        plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
        plt.errorbar(optical_obs_data[0], optical_obs_data[1], fmt='o', yerr=optical_obs_errors)
        plt.errorbar(x_rays_obs_data[0], x_rays_obs_data[1], fmt='o', yerr=x_rays_obs_errors)
        plt.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], fmt='o', yerr=gamma_rays_obs_errors)

        plt.savefig("emcee_best_at_step_{:03d}.png".format(i), dpi='figure')
        plt.clf()
        plt.close()

    def create_corner(i, samples):
        fig = corner.corner(samples, labels=labels, range=limits, truths=real_theta, quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
        fig.savefig("triangle_at_step_{:03d}.png".format(i))
        plt.close(fig)

    def save_data(i, data):
        np.save("data_mcmc_at_step_{:03d}".format(i), data)

    print("Percent\tlog10(g_min)\tlog10(g_break)\tlog10(g_max)\tp_1\t\tp_2\t\tlog10(Sf·R^2)\tlog10(h)\tB\t\tGamma\t\tlog10(rho·h·g_min^2/3)")
    total_iterations = 0

    # First, do 10 quick iterations and output general data to see if
    # everything is behaving correctly
    for i in range(10):
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 1, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        blobs = np.array(blobs)

        lnprobs_max_index = (lnprobs.argsort())[::-1]
        best_lnprob = lnprobs[lnprobs_max_index][0]
        best_pos    = pos[lnprobs_max_index][0]
        best_blob   = blobs[lnprobs_max_index][0]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_fig(i, best_blob)
        create_corner(i, samples)
        save_data(i, sampler.chain)

        total_iterations += 1

    # Now, do some longish iterations as a 'burn-in' thing
    for i in range(10, 29):
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        blobs = np.array(blobs)

        lnprobs_max_index = (lnprobs.argsort())[::-1]
        best_lnprob = lnprobs[lnprobs_max_index][0]
        best_pos    = pos[lnprobs_max_index][0]
        best_blob   = blobs[lnprobs_max_index][0]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_fig(i, best_blob)
        create_corner(i, samples)
        save_data(i, sampler.chain)

        total_iterations += 10
    i = 29

    # total_iterations ~ 200

    # Now, '''Real''' MCMC process. This runs until some condition of steady
    # state is True
    while True:
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        blobs = np.array(blobs)

        lnprobs_max_index = (lnprobs.argsort())[::-1]
        best_lnprob = lnprobs[lnprobs_max_index][0]
        best_pos    = pos[lnprobs_max_index][0]
        best_blob   = blobs[lnprobs_max_index][0]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_fig(i, best_blob)
        create_corner(i, samples)
        save_data(i, sampler.chain)

        total_iterations += 10
        i += 1

        out_conditions = 0
        for dim in range(ndim):
            old_data = sampler.chain[:, -100:-50, dim]
            new_data = sampler.chain[:, -50 :   , dim]

            old_medians = np.percentile(old_data, [16, 50, 84])
            new_medians = np.percentile(new_data, [16, 50, 84])

            if ((new_medians - old_medians) / new_medians < 1e-3).all():
                out_conditions += 1

        if out_conditions == ndim:
            break

    # SIGINT -> KeyboardInterrupt
    # and that means that we are killing the program with Ctrl+C.
    # Just dump the data we have so far and be done with it
except KeyboardInterrupt as e:
    samples = sampler.chain[:,:,:].reshape((-1,ndim))
    np.save("saved_data_mcmc", sampler.chain)
    fig = corner.corner(samples, labels=labels, range=limits, truths=real_theta, quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
    fig.savefig("saved_triangle.png")

    sys.exit(1)

samples = sampler.chain[:,:,:].reshape((-1,ndim))
np.save("final_data_mcmc", sampler.chain)
fig = corner.corner(samples, labels=labels, range=limits, truths=real_theta, quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
fig.savefig("final_triangle.png")
