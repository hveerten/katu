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

This will compile a static library in the folder `src`, three examples in
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

`
 .
 |- executable
 |- default.toml
 |
 |-tables
 |   |
 |   |- BH_KA_table.tab`

Launch any of the examples by passing `default.toml` (Or any other kind of
configuration file, to be explained later) as first argument.

`example_simple` will simply simulate the blob until steady state and output
the photon population at the end. `example_data` will also create a folder
called `data` where it will dump information about every particle species
every 10 steps. `example_graphics` will also use gnuplot to plot the resulting
particle species every two steps.

## Configuration

Configuring Katu is done with the [Toml](https://toml.io/en/) language. An
example configuration can be found in the `default.toml` file.
