import sys

from importlib import import_module

import numpy as np
import matplotlib.pyplot as plt
import os

import corner


from utils.analyzer_utils import  \
    plot_bin_evolution,     \
    plot_lnprob_evolution,  \
    plot_random_walkers,    \
    plot_histogram,         \
    plot_correlation,       \
    create_correlation,     \
    create_histogram,       \
    create_corner,          \
    create_report

def write_config_file(theta, config_filename):
    config = theta_to_config(theta)
    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

levels = 1.0 - np.exp(-0.5 * np.linspace(0.4, 3.5, 10) ** 2)

Object = "3C279"

volume = "sphere"
model  = "injection"

directory     = "{}_results".format(Object)
sub_directory = "emcee_{}_{}".format(volume, model)
path          = directory + '/' + sub_directory

prefix = "{}/data_".format(path)

metadata_filename = directory + '.' + sub_directory + '.metadata'

data_filename = "{}final.npz".format(prefix)
data = np.load(data_filename)

# Load things from Metadata
metadata = import_module(metadata_filename)

discarded_data     = getattr(metadata, 'discarded_data', 0)
labels             = getattr(metadata, 'labels')
theta_to_config    = getattr(metadata, 'theta_to_config')


# Load data
full_data = np.load(data_filename)
raw_data = full_data['data']
lnprobs  = full_data['lnprobs']

walkers, total_iterations, ndim = raw_data.shape

if volume == "sphere":      print("Processing data from a spherical model")
if volume == "shell":       print("Processing data from a shell model")
if model == "steady_state": print("Processing data from a steady state model")
if model == "injection":    print("Processing data from an injection model")

# Discard data
raw_data = raw_data[:, discarded_data:, :]
lnprobs  = lnprobs[:, discarded_data:]

print("Data after discarding:\t", raw_data.shape)
print()

# Reshape the data to be usable
raw_data_samples = raw_data.reshape((-1, ndim))
lnprobs_samples  = lnprobs.reshape(-1)

# Find best point
best_lnprob_index = np.where(lnprobs_samples == lnprobs_samples.max())
best_theta        = raw_data_samples[best_lnprob_index][0]

print("best lnprob:\t", lnprobs.max())
print("best theta:\t",  best_theta)

create_corner(raw_data_samples, labels, best_theta, levels, prefix + 'final.png')

print()
create_report(raw_data_samples, labels, best_theta)

write_config_file(best_theta, "{}best_parameters.toml".format(prefix))
