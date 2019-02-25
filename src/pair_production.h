#ifndef PAIR_PRODUCTION_H
#define PAIR_PRODUCTION_H
#pragma once

#include "state.h"

void pair_production_process_photon_losses(state_t *st);

void init_pair_production_LUT_photon_losses(state_t *st);
void calculate_pair_production_LUT_photon_losses(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *pair_production_photon_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    pair_production_process_photon_losses(st);
    return NULL;
}

#endif


#endif /* end of include guard: PAIR_PRODUCTION_H */
