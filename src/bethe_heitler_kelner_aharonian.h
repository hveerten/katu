#ifndef BETHE_HEITLER_KELNER_AHARONIAN_H
#define BETHE_HEITLER_KELNER_AHARONIAN_H
#pragma once

#include "state.h"

void bethe_heitler_ka_process_lepton_gains(state_t *st);

void init_bethe_heitler_ka_LUT_lepton_gains(state_t *st);

void calculate_bethe_heitler_ka_LUT_lepton_gains(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *bethe_heitler_ka_lepton_gains_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    bethe_heitler_ka_process_lepton_gains(st);
    return NULL;
}

#endif

#endif /* end of include guard: BETHE_HEITLER_KELNER_AHARONIAN_H */
