import sys

from importlib import import_module

import numpy as np
import matplotlib.pyplot as plt
import os

from utils.multinest_analyzer_utils import \
        create_report, \
        create_corner

def write_config_file(theta, config_filename):
    config = theta_to_config(theta)
    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

levels = 1.0 - np.exp(-0.5 * np.linspace(0.4, 3.5, 10) ** 2)

plt.rcdefaults()

Object = "3C279"

volume = "sphere"
model  = "injection"

directory     = "{}_results".format(Object)
sub_directory = "multinest_{}_{}".format(volume, model)
path          = directory + '/' + sub_directory

prefix = "{}/data_".format(path)

metadata_filename = directory + '.' + sub_directory + '.metadata'
metadata_imports = import_module(metadata_filename)

labels          = getattr(metadata_imports, 'labels')
theta_to_config = getattr(metadata_imports, 'theta_to_config')


data_filename = "{}.txt".format(prefix)
data = np.loadtxt(data_filename)

probability   =  data[:,0]
loglikelihood = -data[:,1]/2
parameters    =  data[:,2:]

n_dim = parameters.shape[1]

print()
best_parameters    = parameters   [loglikelihood == np.max(loglikelihood)][0]
best_loglikelihood = loglikelihood[loglikelihood == np.max(loglikelihood)]
print("Best Parameters:\t{}".format(best_parameters))
print("Best Loglikelihood:\t{}".format(best_loglikelihood))

write_config_file(best_parameters, "{}best_parameters.toml".format(prefix))
create_report(probability, parameters, labels, best_parameters)

mask = probability > 1e-30

create_corner(
        parameters[mask,:],
        weights[mask],
        labels,
        best_parameters,
        levels,
        prefix + 'corner.png')
