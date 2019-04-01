#ifndef POPULATION_H
#define POPULATION_H
#pragma once

#include "distribution.h"

enum particle_type
{
    photon,
    electron,
    positron,
    proton,
    neutron,
    positive_pion,
    negative_pion,
    neutral_pion,
    positive_left_muon,
    positive_right_muon,
    negative_left_muon,
    negative_right_muon,
    electron_neutrino,
    electron_antineutrino,
    muon_neutrino,
    muon_antineutrino
};

typedef struct
{
    double *population, *log_population;
    double *energy,     *log_energy;

    double *tentative_population;

    double mass;
    unsigned int size;
} population_t;

void init_population(population_t *p, enum particle_type pt,
        double min_energy, double max_energy, unsigned int size);

void free_population(population_t *p);

void generate_population(population_t *p, distribution_metadata_t *dm);

double calculate_population(population_t *p);
double calculate_energy(population_t *p);

/*
 *void init_population(state_t *st, enum particle_type pt,
 *        double min_energy, double max_energy, unsigned int size);
 *
 *void generate_population(state_t *st,
 *        enum particle_type pt, enum distribution_type dt,
 *        double *params);
 */

#endif /* end of include guard: POPULATION_H */
