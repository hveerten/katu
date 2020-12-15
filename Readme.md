# Katu

## Introduction

Katu is a library aimed to simulate the interaction of particles in plasma.
It currently traces protons, neutrons, leptons, pions and neutrinos.

## Installation

Katu requires `gsl` with `cblas` support and `pthreads` as dependencies and
`cmake` as a make dependency.

Building the code is very simple and just requires creating a new directory:

`mkdir build`

configuring using cmake:

`cd build && cmake ..`

and making:

`make`

This will compile a static library in the folder `src`, four examples in
the folder `examples` and a util in the `tables` folder.

## Usage

Katu is mainly intended to be used as a library, but the examples can be used
right away to produce results. In other cases, you can use them as simple
templates from where to write your own programs.

Katu produces a lot of look-up tables (LUTs) and, even if I have tryed hard
to optimize its creation, in some cases this has been impossible. So, the first
step to launch any of the examples, or any program that relies on Katu is to
create the Bethe-Heitler LUT. For that, you can use the program that should
have been compiled in the `tables` folder.

Before launching any of the examples, you have to make sure that the table
that the `generate_table_BH_ka` program has created is called `BH_KA_table.tab`
and is situated in a folder called `tables` in the same directory as the
executable. Likewise, the file `default.toml` must be in the same
directory. Like this:

```
 .
 |- executable
 |- default.toml
 |
 |-tables
 |   |
 |   |- BH_KA_table.tab
```

Launch any of the examples by passing `default.toml` (Or any other kind of
configuration file, to be explained later) as first argument.

`example_simple` will simply simulate the blob until steady state and output
the photon population at the end. `example_data` will also create a folder
called `data` where it will dump information about every particle species
every 10 steps. `example_graphics` will also use gnuplot to plot the resulting
particle species every two steps.

`example_spheres` is a bit different as it requires two configuration files
to represent concentric spheres. This example will simulate two blobs, one
inside the other and will output the photon population of the outer one
at the end

## Configuration

Configuring Katu is done with the [Toml](https://toml.io/en/) language. An
example configuration can be found in the `default.toml` file. For the sake
of completeness, an explanation of all the options is given here.

The `[general]` table refers to general options for the simulation and the
plasma. For now, this table comprises both physical and numerical keys,
although that might change in the future. The possible keys are:

* `density`: The density of the background plasma
* `magnetic_field`: The magnitude of the background magnetic field
* `eta`: The ratio of protons to electrons in the background plasma
* `cfe_ratio`: The ratio of charged to free escape
* `t_acc`: The acceleration timescale
* `dt`: The initial value of the time step
* `dt_max`: The maximum value of the time step
* `t_max`: The maximum simulation time
* `tol`: The requested tolerance for the steady state

The `[volume]` table refers to options refered to the shape of the simulation
volume.

* `shape`: The shape of the volume, currently it supports `"sphere"` and `"disk"`
* `R`: The radius of the volume
* `h`: The height of the disk. (Only used for the `disk` shape)

The `external_injection` table refers to the injection of particles from
outside the simulated volume. The possible keys are:

* `luminosity`: The integrated luminosity of the injected particles
* `eta`: The ratio of protons to electrons to inject.

Note that the subtable `external_injection.photons` also has `luminosity`
as a key separated from the particle injection.

The distribution tables and subtables (`[protons]`, `[electrons]` and
`[photons]`) refer to the distribution of particles, their limits and number
of bins. The available keys are:

* `gamma_min`: The minimum value for gamma for protons and electrons
* `gamma_max`: The maximum value for gamma for protons and electrons
* `epsilon_min`: The minimum value for epsilon for photons
* `epsilon_max`: The maximum value for epsilon for photons
* `size`: The number of bins for the particle kind

Note that the distribution limits for the rest of the particle species are
derived from these.

The distribution itself is controlled by the key `distribution_type` and
supports the following options, along with the corresponding keys:

* `"maxwell_juttner"`: This option gives a thermal distribution of particles
    * `temperature`: The temperature of the distribution, in terms of the rest
                     mass energy
* `"black_body"`: This option gives a thermal distribution of photons
    * `temperature`: The temperature of the distribution, in terms of the rest
                     mass energy of an electron
* `"power_law"`: This option gives a simple power law
    * `slope`: The value of the exponent of the power law
* `"broken_power_law"`: This option gives a broken power law
    * `break_point`: The value of the energy where the the power law changes
                     from the first slope to the second
    * `first_slope`: The value of the first exponent of the power law
    * `second_slope`: The value of the second exponent of the power law
* `"connected_power_law"`: This option gives a smooth power law with two slopes
    * `connection_point`: The value of the energy where the the power law changes
                     from the first slope to the second
    * `first_slope`: The value of the first exponent of the power law
    * `second_slope`: The value of the second exponent of the power law
* `"power_law_with_exponential_cutoff"`: This option gives a power law with an
                                          exponential cutoff
    * `slope`: The value of the exponent of the power law
    * `break_point`: The normalization value for the energy in the exponential
                     cutoff
* `"hybrid"`: This option gives a thermal distribution of particles at low
              energies coupled with a power law at higher energies
    * `temperature`: The temperature of the distribution, in terms of the rest
                     mass energy
    * `slope`: The value of the exponent of the power law


## Wrappers

Katu comes with two wrappers to allow it to work with two common packages to
do bayesian analysis: `emcee` and `multinest`.

`emcee` is described in [`emcee: the MCMC Hammer`](https://arxiv.org/abs/1202.3665)
and implements Markov Chain Monte Carlo analysis with Affine Invariant ensembles.
The interface we employ for this is the python package [`emcee`](https://github.com/dfm/emcee) 
by Dan Foreman-Mackey.

`multinest` is described in [`Importance Nested Sampling and the MultiNest Algorithm`](http://arxiv.org/abs/1306.2144)
and implements Nested Sampling analysis. The interface we employ for this is
the python package [`pymultinest`](https://github.com/JohannesBuchner/PyMultiNest) by
Johannes Buchner.

In the case of using `emcee` the corresponding package has to be installed by doing

`pip install --user emcee`

In the case of using `multinest`, then multinest itself must be installed from
the distribution's repositories or from multinest's [source](https://github.com/rjw57/MultiNest)
and `pymultinest` must be installed:

`pip install --user pymultinest`

Note that in some cases, manual setting of the path to the multinest library
must be done.


The usage of the wrappers is a bit complex because the user is expected to
choose and write the model and the functions that go from the analyzed space
to a Katu configuration file. Some examples can be seen in the
`wrapper_utils/utils_{shell,sphere}_{steady_state,injection}.py` files.

In particular, the user must supply the variables:

* `ndim`: The number of dimensions of the space
* `labels`: A list with the names of the variables
* `lnprior`: A function that returns the a priori probability for the chosen point
* `get_gamma`: A function that returns the value of Gamma from a point in space
* `get_R`: A function that returns the value of R from a point in space

For `emcee`, the user must supply this aditional variable:

* `initial_theta`: A 3D array with `ndim` entries that sets the spread (0,2) and
    the central value (1) for the variables to start from.

For `multinest`, the user must supply the aditional function:

* `prior`: A function that transforms from the unit cube to the space coordinates.

## Publications

Non-eshaustive list of publications that use Katu in any form:

* A Bayesian Approach to Modelling Multi-Messenger Emission from Blazars using
    Lepto-Hadronic Kinetic Equations [ADS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.500.3613J/abstract) [arXiv](https://arxiv.org/abs/2006.01543v2)
