#ifndef BETHE_HEITLER_H
#define BETHE_HEITLER_H
#pragma once

#include "state.h"

void bethe_heitler_process_lepton_gains(state_t *st);

void init_bethe_heitler_LUT_lepton_gains(state_t *st);

void calculate_bethe_heitler_LUT_lepton_gains(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *bethe_heitler_lepton_gains_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    bethe_heitler_process_lepton_gains(st);
    return NULL;
}

#endif

#endif /* end of include guard: BETHE_HEITLER_H */
