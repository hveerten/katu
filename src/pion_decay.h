#ifndef PION_DECAY_H
#define PION_DECAY_H
#pragma once

#include "state.h"

void charged_pion_decay(state_t *st);

void init_pion_decay_LUT_lifetime(state_t *st);
void init_pion_decay_LUT_muon_functions(state_t *st);

void calculate_pion_decay_LUT_lifetime(state_t *st);
void calculate_pion_decay_LUT_muon_functions(state_t *st);

void neutral_pion_decay(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *charged_pion_decay_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    charged_pion_decay(st);
    return NULL;
}

#endif

#endif /* end of include guard: PION_DECAY_H */
