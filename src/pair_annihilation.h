#ifndef PAIR_ANNIHILATION_H
#define PAIR_ANNIHILATION_H
#pragma once

#include "state.h"

void pair_annihilation_process_photon_gains(state_t *st);
void pair_annihilation_process_lepton_losses(state_t *st);

void init_pair_annihilation_LUT_photon_gains(state_t *st);
void init_pair_annihilation_LUT_lepton_losses(state_t *st);

void calculate_pair_annihilation_LUT_photon_gains(state_t *st);
void calculate_pair_annihilation_LUT_lepton_losses(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *pair_annihilation_photon_gains_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    pair_annihilation_process_photon_gains(st);
    return NULL;
}

static inline void *pair_annihilation_lepton_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    pair_annihilation_process_lepton_losses(st);
    return NULL;
}

#endif


#endif /* end of include guard: PAIR_ANNIHILATION_H */
