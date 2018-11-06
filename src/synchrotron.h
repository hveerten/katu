#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H
#pragma once

#include "population.h"
/*#include "state.h"*/

struct state_t;
typedef struct synchrotron_t
{
    double photon_gains_factor;
    double photon_losses_factor;
    double particle_losses_factor;
    double nu_0;

    double dt; // Needed for implicit schemes

    double *temp1, *temp2;  // Needed for the derivatives!

    double *photon_gains;
    double *photon_losses;
    double *particle_losses;

    population_t *particles;
    population_t *photons;

    double *LUT_xCS;

#if USE_THREADS
    void *(*synchrotron_function)(void *args);
#else
    void (*synchrotron_function)(struct synchrotron_t *synchro);
#endif
} synchrotron_t;

void init_synchrotron(struct state_t *st, synchrotron_t *synchro, enum particle_type pt);

void synchrotron_update_dt(synchrotron_t *synchro, double dt);
void synchrotron_update_B(synchrotron_t *synchro,  double B);

#endif /* end of include guard: SYNCHROTRON_H */
