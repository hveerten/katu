#ifndef PION_PRODUCTION_H
#define PION_PRODUCTION_H
#pragma once

#include "state.h"

void multi_resonance_pion_production(state_t *st);
void multi_resonance_pion_production_hadron_gains(state_t *st);
void multi_resonance_pion_production_hadron_losses(state_t *st);

void direct_pion_production(state_t *st);
void direct_pion_production_hadron_gains(state_t *st);
void direct_pion_production_hadron_losses(state_t *st);


void init_multi_resonances_LUT_pion_gains(state_t *st);
void init_multi_resonances_LUT_hadron_losses(state_t *st);
void init_multi_resonances_LUT_hadron_gains(state_t *st);

void init_direct_LUT_pion_gains(state_t *st);
void init_direct_pion_production_LUT_hadron_losses(state_t *st);
void init_direct_pion_production_LUT_hadron_gains(state_t *st);

void calculate_multi_resonances_LUT_pion_gains(state_t *st);
void calculate_multi_resonances_LUT_hadron_losses(state_t *st);
void calculate_multi_resonances_LUT_hadron_gains(state_t *st);

void calculate_direct_LUT_pion_gains(state_t *st);
void calculate_direct_pion_production_LUT_hadron_losses(state_t *st);
void calculate_direct_pion_production_LUT_hadron_gains(state_t *st);

#ifdef USE_THREADS
#include <stdlib.h>

static inline void *multi_resonance_pion_production_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    multi_resonance_pion_production(st);
    return NULL;
}

static inline void *multi_resonance_pion_production_hadron_gains_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    multi_resonance_pion_production_hadron_gains(st);
    return NULL;
}

static inline void *multi_resonance_pion_production_hadron_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    multi_resonance_pion_production_hadron_losses(st);
    return NULL;
}

static inline void *direct_pion_production_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    direct_pion_production(st);
    return NULL;
}

static inline void *direct_pion_production_hadron_gains_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    direct_pion_production_hadron_gains(st);
    return NULL;
}

static inline void *direct_pion_production_hadron_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    direct_pion_production_hadron_losses(st);
    return NULL;
}
#endif

#endif /* end of include guard: PION_PRODUCTION_H */
