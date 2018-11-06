#ifndef MUON_DECAY_H
#define MUON_DECAY_H
#pragma once

#include "state.h"

void muon_decay(state_t *st);

void init_muon_decay_LUT_lifetime(state_t *st);
void init_muon_decay_LUT_neutrino_functions(state_t *st);

void calculate_muon_decay_LUT_lifetime(state_t *st);
void calculate_muon_decay_LUT_neutrino_functions(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *muon_decay_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    muon_decay(st);
    return NULL;
}

#endif

#endif /* end of include guard: MUON_DECAY_H */
