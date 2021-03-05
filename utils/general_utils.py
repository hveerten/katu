import subprocess
import os
from io import StringIO

import scipy.integrate as integ
import numpy as np

ELECTRON_MASS   = 9.109e-28
LIGHT_SPEED     = 2.99792458e10
ELECTRON_ENERGY = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED

# H_0 is km/s / Mpc
H_0 = 67.66
Omega_m = 0.3111
Omega_l = 0.6889

def setup_functions(theta_to_config):
    global _internal_theta_to_config
    _internal_theta_to_config = theta_to_config

def calculate_luminosity_distance(z):
    return (1 + z) * LIGHT_SPEED / (H_0 * 1e5) * \
            integ.quad(lambda x: 1 / np.sqrt((1 + x)**3 * Omega_m + Omega_l), 0, z)[0] * 3.086e24

def write_config_file(theta, config_filename):
    config = _internal_theta_to_config(theta)
    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

# Simple version of `generate_sim_data` which only returns
# the data that it gets from stdout. This can be used as a basis for
# customized get_sim_data functions
def generate_sim_data(theta):
    config_filename = "config_{}.toml".format(os.getpid())

    write_config_file(theta, config_filename)
    process = subprocess.Popen(['./example_wrapper',config_filename],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)

    output_data, err = process.communicate()
    data             = np.loadtxt(StringIO(output_data), unpack=True)

    os.remove(config_filename)

    return data

# sim_data is a 2,n ndarray
# obs_data is a 3,n ndarray
#
# No unit conversion done internally, so they must be compatible.
# The third component of obs_data is the error
def general_lnlike(sim_data, obs_data):

    interp_sim_data = np.exp(np.interp(np.log(obs_data[0]), np.log(sim_data[0]), np.log(sim_data[1])))

    lnlike = -np.sum(0.5*np.power((interp_sim_data - obs_data[1]) / (obs_data[2]), 2))

    return lnlike

# Returns -infinity if any point is above the upper limits given by obs_data
def general_upper_limit(sim_data, obs_data):

    interp_sim_data = np.exp(np.interp(np.log(obs_data[0]), np.log(sim_data[0]), np.log(sim_data[1])))

    if np.any(interp_sim_data > obs_data[1]):
        return -np.inf

    return 0
