#ifndef INVERSE_COMPTON_H
#define INVERSE_COMPTON_H
#pragma once

#include "state.h"

void inverse_compton_process_head_on(state_t *st);
void inverse_compton_process_head_on_upscattering(state_t *st);
void inverse_compton_process_head_on_downscattering(state_t *st);

void inverse_compton_process_photon_losses(state_t *st);
void inverse_compton_process_electron_losses(state_t *st);

void init_inverse_compton_LUT_g1_g2(state_t *st);
void init_inverse_compton_LUT_losses_reaction_rate(state_t *st);
void calculate_inverse_compton_LUT_g1_g2(state_t *st);
void calculate_inverse_compton_LUT_losses_reaction_rate(state_t *st);


#ifdef USE_THREADS
#include <stdlib.h>

static inline void *inverse_compton_head_on_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    inverse_compton_process_head_on(st);
    return NULL;
}

static inline void *inverse_compton_head_on_upscattering_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    inverse_compton_process_head_on_upscattering(st);
    return NULL;
}

static inline void *inverse_compton_head_on_downscattering_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    inverse_compton_process_head_on_downscattering(st);
    return NULL;
}

static inline void *inverse_compton_photon_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    inverse_compton_process_photon_losses(st);
    return NULL;
}

static inline void *inverse_compton_electron_losses_wrapper(void *args)
{
    state_t *st = (state_t *)args;
    inverse_compton_process_electron_losses(st);
    return NULL;
}
#endif

#endif /* end of include guard: INVERSE_COMPTON_H */
