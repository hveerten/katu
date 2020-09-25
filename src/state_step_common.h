#ifndef STATE_STEP_COMMON_H
#define STATE_STEP_COMMON_H
#pragma once

#include "state.h"

#include <stdbool.h>

#if LEPTON_STEADY_STATE == 0
#define APPLY_TO_ELECTRONS STEP_INTERNAL_FUNCTION(electrons)
#define APPLY_TO_POSITRONS STEP_INTERNAL_FUNCTION(positrons)
#else
#define APPLY_TO_ELECTRONS
#define APPLY_TO_POSITRONS
#endif

#if PROTON_STEADY_STATE == 0
#define APPLY_TO_PROTONS STEP_INTERNAL_FUNCTION(protons)
#else
#define APPLY_TO_PROTONS
#endif

#define APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES   \
    STEP_INTERNAL_FUNCTION(photons)                 \
    APPLY_TO_ELECTRONS \
    APPLY_TO_POSITRONS \
    APPLY_TO_PROTONS \
    STEP_INTERNAL_FUNCTION(neutrons)                \
                                                    \
    STEP_INTERNAL_FUNCTION(neutral_pions)           \
    STEP_INTERNAL_FUNCTION(positive_pions)          \
    STEP_INTERNAL_FUNCTION(negative_pions)          \
                                                    \
    STEP_INTERNAL_FUNCTION(positive_left_muons)     \
    STEP_INTERNAL_FUNCTION(positive_right_muons)    \
    STEP_INTERNAL_FUNCTION(negative_left_muons)     \
    STEP_INTERNAL_FUNCTION(negative_right_muons)    \
                                                    \
    STEP_INTERNAL_FUNCTION(electron_neutrinos)      \
    STEP_INTERNAL_FUNCTION(electron_antineutrinos)  \
    STEP_INTERNAL_FUNCTION(muon_neutrinos)          \
    STEP_INTERNAL_FUNCTION(muon_antineutrinos)

void step_calculate_processes(state_t *st);

void step_update_populations(state_t *st, double dt);

void step_experimental_update_populations(state_t *st, double dt);
void step_experimental_update_populations_injection(state_t *st, double dt);

bool step_check_tentative_populations(state_t *st);
void step_fix_tentative_zeros(state_t *st);
void step_accept_tentative_population(state_t *st);

void step_log_populations(state_t *st);

void step_propagate_new_dt(state_t *st, double dt);

void step_report_population_update(state_t *st, enum particle_type pt);

#endif /* end of include guard: STATE_STEP_COMMON_H */
