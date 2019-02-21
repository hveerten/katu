import emcee

import subprocess
import numpy as np
import scipy.integrate as integ
from io import StringIO

import os
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner

from emcee_wrapper_utils.obs_data import     \
    radio_obs_data,      \
    optical_obs_data,    \
    x_rays_obs_data,     \
    gamma_rays_obs_data, \
    optical_obs_errors,  \
    x_rays_obs_errors,   \
    gamma_rays_obs_errors

z = 0.3365
ELECTRON_MASS   = 9.109e-28
LIGHT_SPEED     = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED

H_0 = 67.66
Omega_m = 0.3111
Omega_l = 0.6889

# H_0 is km/s / Mpc
D_l = (1 + z) * LIGHT_SPEED / (H_0 * 1e5) * integ.quad(lambda x: 1 / np.sqrt((1 + x)**3 * Omega_m + Omega_l), 0, z)[0] * 3.086e24

volume = "shell"
model  = "steady_state"

assert(volume == "shell"        or volume == "sphere")
assert(model  == "steady_state" or model  == "injection")

utils_module = "emcee_wrapper_utils.utils_{}_{}".format(volume, model)
utils = import_module(utils_module)

ndim          = getattr(utils, 'ndim')
labels        = getattr(utils, 'labels')
initial_theta = getattr(utils, 'initial_theta')

lnprior = getattr(utils, 'lnprior')

theta_to_config = getattr(utils, 'theta_to_config')
theta_to_params = getattr(utils, 'theta_to_params')
params_to_theta = getattr(utils, 'params_to_theta')

get_gamma = getattr(utils, 'get_gamma')
get_R     = getattr(utils, 'get_R')

S_f_prefactor = ELECTRON_ENERGY * LIGHT_SPEED / D_l**2
if   volume == "shell":  S_f_prefactor *= 1 / np.pi
elif volume == "sphere": S_f_prefactor *= 4 / 9

walkers = 128

def write_config_file(theta, config_filename):
    config = theta_to_config(theta)
    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

def generate_sim_data(theta):
    config_filename = "config_emcee_{}.toml".format(os.getpid())

    write_config_file(theta, config_filename)
    process = subprocess.Popen(['./model',config_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output_data, err = process.communicate()
    data             = np.loadtxt(StringIO(output_data), unpack=True)
    energy_injection = np.loadtxt(StringIO(err), unpack=True)
    os.remove(config_filename)

    return data, energy_injection

def lnlike(theta):
    Gamma = get_gamma(theta)
    R     = get_R(theta)

    S_f = S_f_prefactor (R * Gamma)**2

    sim_data, energy_injection = generate_sim_data(theta)

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

    # If one of the energy injections is NaN return 'impossible'
    if np.any(np.isnan(energy_injection)):
        return -np.inf, (sim_data, (0,0,0), (0,0))

    # The radio band is an upper limit
    if np.any(radio_sim_data > radio_obs_data[1]): return -np.inf, (sim_data, (0,0,0), (0,0))

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
           (sim_data,                                           \
            (optical_lnlike, x_rays_lnlike, gamma_rays_lnlike), \
            energy_injection)

def lnprob(theta):
    lp = lnprior(theta)
    if np.isinf(lp):
        return -np.inf, (np.array([[0]*128,[0]*128]), (0,0,0))
    else:
        lnlike_data = lnlike(theta)
        return lp + lnlike_data[0], lnlike_data[1]

initial_params = theta_to_params(initial_theta[1])
print(initial_theta[1])
print(initial_params)

# Assert that theta_to_params and params_to_theta are actually inverse functions
np.testing.assert_allclose(params_to_theta(initial_params),initial_theta[1])

input_parameters_lnprob, input_parameters_data = lnprob(initial_theta[1])

plt.figure(figsize=(24,12), dpi=128)
plt.ylim([1e-17, 1e-9])
plt.xlim([1e-6, 1e12])
plt.grid()
plt.loglog(input_parameters_data[0][0], input_parameters_data[0][1])
plt.loglog(radio_obs_data[0], radio_obs_data[1], 'o')
plt.errorbar(optical_obs_data[0], optical_obs_data[1], fmt='o', yerr=optical_obs_errors)
plt.errorbar(x_rays_obs_data[0], x_rays_obs_data[1], fmt='o', yerr=x_rays_obs_errors)
plt.errorbar(gamma_rays_obs_data[0], gamma_rays_obs_data[1], fmt='o', yerr=gamma_rays_obs_errors)

plt.savefig("input_parameters.png", dpi='figure')
plt.clf()
plt.close()
print(input_parameters_lnprob, input_parameters_data[1])

pos = np.array([np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 4) for i in range(walkers)])
pos[0] = initial_theta[1]

# Make sure that ALL walkers start from a valid position
for i in range(walkers):
    prob, _ = lnprob(pos[i])

    while np.isinf(prob):
        pos[i] = np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 4)
        prob, _ = lnprob(pos[i])

    print("\r[{}{}] {:.2g}%".format("=" * int(80 * i/walkers), ' ' * int(80 * (walkers - i)/walkers), 100 * i/walkers), end='')
print()

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

    def print_best(i, best, best_lnprob, best_lnprobs):
        print("BEST:", end='\t')
        for b in best:
            print("{:.3f}".format(b), end='\t')
        print("{:.3f}".format(best_lnprob), end='\t')
        print("({:.3f} + {:.3f} + {:.3f})".format(*best_lnprobs))
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
        tt = time.time()
        print("Creating Corner Plot", end='\t')
        fig = corner.corner(samples, labels=labels, truths=initial_theta[1], quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
        fig.savefig("triangle_at_step_{:03d}.png".format(i))
        plt.close(fig)
        print("Done\t({} seconds)".format(time.time() - tt))

    def save_data(i, data, lnprobs, blobs):
        tt = time.time()
        print("Saving Data", end='\t')
        energy_injection = [[blobs[i][w][2] for w in range(len(blobs[i]))] for i in range(len(blobs))]
        energy_injection = np.array(energy_injection)
        np.savez("data_mcmc_at_step_{:03d}".format(i), data=data, lnprobs=lnprobs, energy_injection=energy_injection)
        print("Done\t({} seconds)".format(time.time() - tt))

    print("Percent\tlog10(g_min)\tlog10(g_break)\tlog10(g_max)\tp_1\t\tp_2\t\tlog10(B)\tlog10(Gamma)\tlog10(h)\tlog10(R·g_break^-3/2)\tlog10(rho·h·g_min^2/3)")
    total_iterations = 0

    # First, do 10 quick iterations and output general data to see if
    # everything is behaving correctly
    for i in range(10):
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 1, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]
        best_blob   = blobs[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob, best_blob[1])
        create_fig(i, best_blob[0])
        create_corner(i, samples)
        save_data(i, sampler.chain, sampler.lnprobability, sampler.blobs)

        total_iterations += 1

    # Now, do some longish iterations as a 'burn-in' thing
    # and check that no walker has gotten stuck in a bad
    # position
    for i in range(10, 29):
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]
        best_blob   = blobs[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob, best_blob[1])
        create_fig(i, best_blob[0])
        create_corner(i, samples)
        save_data(i, sampler.chain, sampler.lnprobability, sampler.blobs)

        for j in range(walkers):
            if (sampler.chain[j, total_iterations:, :] == sampler.chain[j, total_iterations , :]).all() and \
                lnprobs[j] < 25 * best_lnprob:
                print("Walker {} seems to be stuck at {} with {}, restarting".format(j, pos[j], lnprobs[j]))
                lnprobs[j] = -np.inf
                while lnprobs[j] < 25 * best_lnprob or np.isinf(lnprobs[j]):
                    pos[j] = np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 2)
                    lnprobs[j], _ = lnprob(pos[j])

        total_iterations += 10

    i = 29

    # total_iterations ~ 200

    # Now, '''Real''' MCMC process. This runs until some condition of steady
    # state is True
    while True:
        pos, lnprobs, rstate, blobs = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs, blobs0=blobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]
        best_blob   = blobs[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob, best_blob[1])
        create_fig(i, best_blob[0])
        create_corner(i, samples)
        save_data(i, sampler.chain, sampler.lnprobability, sampler.blobs)

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
    np.savez("saved_data_mcmc", data=sampler.chain, lnprobs=sampler.lnprobability)
    fig = corner.corner(samples, labels=labels, truths=initial_theta[1], quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
    fig.savefig("saved_triangle.png")

    sys.exit(1)

samples = sampler.chain[:,:,:].reshape((-1,ndim))
np.savez("final_data_mcmc", data=sampler.chain, lnprobs=sampler.lnprobability)
fig = corner.corner(samples, labels=labels, truths=initial_theta[1], quantiles=[0.16,0.5,0.84], show_titles=True, quiet=True)
fig.savefig("final_triangle.png")
