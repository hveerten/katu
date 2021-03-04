from importlib import import_module

import subprocess
import scipy.integrate as integ
import numpy as np
from io import StringIO

import os
import sys
import json

from utils.general_utils import *

Object = "3C279"
Object_data_file = "3C279_2014_A"

volume = "sphere"
model  = "injection"

assert(volume == "shell"        or volume == "sphere")
assert(model  == "steady_state" or model  == "injection")


directory     = "{}_results".format(Object)
sub_directory = "multinest_{}_{}".format(volume, model)
path          = directory + '/' + sub_directory

try:
    os.makedirs(path)
except OSError: pass

# name of the output files
prefix = "{}/data_".format(path)


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

lnprior = getattr(utils, 'lnprior')

theta_to_config = getattr(utils, 'theta_to_config')

get_gamma = getattr(utils, 'get_gamma')
get_R     = getattr(utils, 'get_R')

assert_consistency = getattr(utils, 'assert_consistency')

setup_functions(theta_to_config)
# End of importing things from utils_module

D_l = calculate_luminosity_distance(z)

S_f_prefactor = ELECTRON_ENERGY * LIGHT_SPEED / D_l**2 / (1 + z)
if   volume == "shell":  S_f_prefactor *= 1 / np.pi
elif volume == "sphere": S_f_prefactor *= 4 / 9


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

# Transform the hypercube to theta-space
def prior(hypercube):
    _log_g_min,  \
    _log_g_max,  \
    _p,          \
    _log_B,      \
    _log_Gamma,  \
    _log_R,      \
    _log_L = hypercube

    lerp = lambda x,a,b : (1-x)*a + x*b

    return np.array([lerp(_log_g_min, 1,              10),
                     lerp(_log_g_max, 1,              10),
                     lerp(_p,         1,              6),
                     lerp(_log_B,     np.log10(5e-3), np.log10(5)),
                     lerp(_log_Gamma, np.log10(2),    np.log10(192)),
                     lerp(_log_R,     10,             26),
                     lerp(_log_L,     38,             50)])

# run MultiNest
result = solve(
        LogLikelihood=lnprob, Prior=prior,
        n_dims=ndim, outputfiles_basename=prefix,
        verbose=True, n_iter_before_update=10)

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for name, col in zip(labels, result['samples'].transpose()):
        print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))

# make marginal plots by running:
# $ python multinest_marginals.py chains/3-
# For that, we need to store the parameter names:
with open('%sparams.json' % prefix, 'w') as f:
        json.dump(labels, f, indent=2)
