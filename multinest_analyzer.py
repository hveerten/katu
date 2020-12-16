from importlib import import_module

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import corner

import toml
import json

def write_config_file(theta, config_filename):
    config = theta_to_config(theta)
    fout = open(config_filename, "w")
    fout.write(config)
    fout.close()

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

median = np.zeros((n_dim))
quantile_16 = np.zeros((n_dim))
quantile_84 = np.zeros((n_dim))

for i in range(n_dim):
    a = np.c_[probability, parameters[:,i]]

    a = a[a[:,1].argsort()]

    a = a.T

    a[0] = a[0].cumsum()

    median[i]      = np.interp(0.50, a[0], a[1])
    quantile_16[i] = np.interp(0.16, a[0], a[1])
    quantile_84[i] = np.interp(0.84, a[0], a[1])

    print("{}:\t{:.2f}^{{+{:.2f}}}_{{-{:.2f}}}".format(
        labels[i],
        median[i],
        quantile_84[i] - median[i],
        median[i] - quantile_16[i]))

print()
best_parameters    = parameters[loglikelihood    == np.max(loglikelihood)][0]
best_loglikelihood = loglikelihood[loglikelihood == np.max(loglikelihood)]
print("Best Parameters:\t{}".format(best_parameters))
print("Best loglikelihood:\t{}".format(best_loglikelihood))

write_config_file(best_parameters, "{}best_parameters.toml".format(prefix))

mask = probability > 1e-30

corner.corner(parameters[mask,:], weights=probability[mask],
        labels=labels, show_titles=True,
        truths=best_parameters)
plt.savefig(prefix + 'corner.png')
plt.close()
