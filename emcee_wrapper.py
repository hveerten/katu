import emcee

from importlib import import_module

import numpy as np
import scipy.integrate as integ
from io import StringIO

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import corner

from utils.general_utils import *

from utils.emcee_utils import *

Object = "3C279"
Object_data_file = "3C279_2014_A"

volume = "shell"
model  = "steady_state"

assert(volume == "shell"        or volume == "sphere")
assert(model  == "steady_state" or model  == "injection")

directory     = "{}_results".format(Object)
sub_directory = "emcee_{}_{}".format(volume, model)
path          = directory + '/' + sub_directory

try:
    os.makedirs(path)
except OSError: pass

# name of the output files
data_prefix   = "{}/data_".format(path)
corner_prefix = "{}/corner_".format(path)

# Import the data from the data_module
data_module = "Object_Data.{}".format(Object_data_file)
data_imports = import_module(data_module)

z        = getattr(data_imports, 'z')
obs_data = getattr(data_imports, 'data')

# Import everything from the utils_module
utils_module = "utils.utils_{}_{}".format(volume, model)
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

assert_consistency = getattr(utils, 'assert_consistency')

setup_functions(theta_to_config)
# End of importing things from utils_module

D_l = calculate_luminosity_distance(z)

S_f_prefactor = ELECTRON_ENERGY * LIGHT_SPEED / D_l**2 / (1 + z)
if   volume == "shell":  S_f_prefactor *= 1 / np.pi
elif volume == "sphere": S_f_prefactor *= 4 / 9

walkers = 128

def lnlike(theta):
    Gamma = get_gamma(theta)
    R     = get_R(theta)

    S_f = S_f_prefactor * (R * Gamma**2)**2

    try:
        sim_data = generate_sim_data(theta)
    except ValueError as e:
        return -np.inf

    # sim_data[0] in epsilon -> eV
    # sim_data[1] in pop     -> erg/cm^2/s
    sim_data[1] *= sim_data[0] * sim_data[0]
    sim_data[0] *= 511000

    # Transform the simulated data from the fluid frame to the obs frame
    sim_data[0] *= Gamma / (1 + z)
    sim_data[1] *= S_f

    return general_lnlike(sim_data, obs_data)

def lnprob(theta):
    lp = lnprior(theta)
    if np.isinf(lp):
        return -np.inf
    else:
        return lp + lnlike(theta)

# Checks for the initial parameters
assert_consistency(initial_theta[1])

initial_params = theta_to_params(initial_theta[1])

input_parameters_lnprob = lnprob(initial_theta[1])

print(initial_theta[1])
print(initial_params)
print(input_parameters_lnprob)

# Create the first walkers
pos = np.array([np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 4) for i in range(walkers)])
pos[0] = initial_theta[1]

# Make sure that ALL walkers start from a valid position
for i in range(walkers):
    prob = lnprob(pos[i])

    while np.isinf(prob):
        pos[i] = np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 4)
        prob = lnprob(pos[i])

    print("\r[{}{}] {:.2g}%".format("=" * int(80 * i/walkers), ' ' * int(80 * (walkers - i)/walkers), 100 * i/walkers), end='')
print()

lnprobs = None
rstate = None
sampler = emcee.EnsembleSampler(walkers, ndim, lnprob, threads=1)

print("Percent\t", "\t".join(labels))
total_iterations = 0

try:
    # First, do 10 quick iterations and output general data to see if
    # everything is behaving correctly
    for i in range(10):
        pos, lnprobs, rstate = sampler.run_mcmc(pos, 1, rstate0=rstate, lnprob0=lnprobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_corner("{}{}.png".format(corner_prefix,i), samples, labels, best_pos)
        save_data("{}{:03d}".format(data_prefix,i), sampler.chain, sampler.lnprobability)

        total_iterations += 1

    # Now, do some longish iterations as a 'burn-in' thing
    # and check that no walker has gotten stuck in a bad
    # position
    for i in range(10, 29):
        pos, lnprobs, rstate = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_corner("{}{}.png".format(corner_prefix,i), samples, labels, best_pos)
        save_data("{}{:03d}".format(data_prefix,i), sampler.chain, sampler.lnprobability)

        for j in range(walkers):
            if (sampler.chain[j, total_iterations:, :] == sampler.chain[j, total_iterations , :]).all() and \
                lnprobs[j] < 25 * best_lnprob:
                print("Walker {} seems to be stuck at {} with {}, restarting".format(j, pos[j], lnprobs[j]))
                lnprobs[j] = -np.inf
                while lnprobs[j] < 25 * best_lnprob or np.isinf(lnprobs[j]):
                    pos[j] = np.random.normal(loc=initial_theta[1], scale=(initial_theta[2] - initial_theta[0]) / 2)
                    lnprobs[j] = lnprob(pos[j])

        total_iterations += 10

    i = 29

    # total_iterations ~ 200

    # Now, '''Real''' MCMC process. This runs until some condition of steady
    # state is True
    while True:
        pos, lnprobs, rstate = sampler.run_mcmc(pos, 10, rstate0=rstate, lnprob0=lnprobs)
        samples = sampler.chain[:,:,:].reshape((-1,ndim))

        lnprobs_max_index = ((lnprobs.argsort())[::-1])[0]

        best_lnprob = lnprobs[lnprobs_max_index]
        best_pos    = pos[lnprobs_max_index]

        print_medians(i, samples, ndim)
        print_best(i, best_pos, best_lnprob)
        create_corner("{}{}.png".format(corner_prefix,i), samples, labels, best_pos)
        save_data("{}{:03d}".format(data_prefix,i), sampler.chain, sampler.lnprobability)

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
    pass

samples = sampler.chain[:,:,:].reshape((-1,ndim))

create_corner("{}final.png".format(corner_prefix), samples, labels, initial_theta[1])
save_data("{}final".format(data_prefix), sampler.chain, sampler.lnprobability)
