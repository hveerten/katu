#ifndef STATE_H
#define STATE_H
#pragma once

#include "acceleration.h"
#include "distribution.h"
#include "escape.h"
#include "population.h"
#include "synchrotron.h"

#include <stdbool.h>

typedef struct
{
    double *stage0_population;
    double *stage_delta[4];
} RK_information_t;

typedef struct state_t
{
    double t, dt, dt_max;

    double B;
    double R;
    double density;

    // For now, we won't use this
    /*double *big_chunk_of_memory;*/

    double photon_energy;

    population_t photons;
    population_t electrons;
    population_t protons;
    population_t neutrons;

    population_t neutral_pions;
    population_t positive_pions;
    population_t negative_pions;

    population_t positive_left_muons;
    population_t positive_right_muons;
    population_t negative_left_muons;
    population_t negative_right_muons;

    population_t electron_neutrinos;
    population_t electron_antineutrinos;
    population_t muon_neutrinos;
    population_t muon_antineutrinos;

    RK_information_t photons_RK_information;
    RK_information_t electrons_RK_information;
    RK_information_t protons_RK_information;
    RK_information_t neutrons_RK_information;

    RK_information_t neutral_pions_RK_information;
    RK_information_t positive_pions_RK_information;
    RK_information_t negative_pions_RK_information;

    RK_information_t positive_left_muons_RK_information;
    RK_information_t positive_right_muons_RK_information;
    RK_information_t negative_left_muons_RK_information;
    RK_information_t negative_right_muons_RK_information;

    RK_information_t electron_neutrinos_RK_information;
    RK_information_t electron_antineutrinos_RK_information;
    RK_information_t muon_neutrinos_RK_information;
    RK_information_t muon_antineutrinos_RK_information;

    double inverse_compton_electron_losses_factor;

    double *inverse_compton_photon_gains;
    double *inverse_compton_photon_gains_downscattering;
    double *inverse_compton_photon_gains_upscattering;
    double *inverse_compton_photon_losses;
    double *inverse_compton_electron_losses;

    double *inverse_compton_inner_downscattering;
    double *inverse_compton_inner_upscattering;

    double *inverse_compton_LUT_g1;
    double *inverse_compton_LUT_g2;
    double *inverse_compton_LUT_downscattering_e_min_index;
    double *inverse_compton_LUT_downscattering_e_max_index;
    double *inverse_compton_LUT_upscattering_e_min_index;
    double *inverse_compton_LUT_upscattering_e_max_index;
    double *inverse_compton_LUT_upscattering_g_min_index;

    double *inverse_compton_LUT_losses_reaction_rate;

    double *pair_production_losses;
    double *pair_production_LUT_R;
    double *pair_production_LUT_index_e_min;

    double *multi_resonances_neutral_pion_gains;
    double *multi_resonances_positive_pion_gains;
    double *multi_resonances_negative_pion_gains;

    double *multi_resonances_proton_losses;
    double *multi_resonances_proton_gains;
    double *multi_resonances_neutron_losses;
    double *multi_resonances_neutron_gains;

    double *multi_resonances_LUT_pion_gains;
    double *multi_resonances_LUT_pion_gains_e_min_index;

    double *multi_resonances_LUT_hadron_gains;
    double *multi_resonances_LUT_hadron_gains_e_min_index;
    double *multi_resonances_LUT_hadron_losses;
    double *multi_resonances_LUT_hadron_losses_e_min_index;

    double *direct_neutral_pion_gains;
    double *direct_positive_pion_gains;
    double *direct_negative_pion_gains;

    double *direct_pion_production_proton_gains;
    double *direct_pion_production_proton_losses;
    double *direct_pion_production_neutron_gains;
    double *direct_pion_production_neutron_losses;

    double *direct_LUT_pion_gains;
    double *direct_LUT_pion_gains_e_min_index;

    double *direct_pion_production_LUT_hadron_gains;
    double *direct_pion_production_LUT_hadron_gains_e_min_index;
    double *direct_pion_production_LUT_hadron_losses;
    double *direct_pion_production_LUT_hadron_losses_e_min_index;

    double *pion_decay_positive_left_muon_gains;
    double *pion_decay_positive_right_muon_gains;
    double *pion_decay_negative_left_muon_gains;
    double *pion_decay_negative_right_muon_gains;

    double *pion_decay_muon_neutrino_gains;
    double *pion_decay_muon_antineutrino_gains;

    double *pion_decay_positive_pion_losses;
    double *pion_decay_negative_pion_losses;

    double *pion_decay_LUT_lifetime;
    double *pion_decay_LUT_positive_right;
    double *pion_decay_LUT_positive_left;

    double *pion_decay_photon_gains;

    double *muon_decay_electron_neutrino_gains;
    double *muon_decay_electron_antineutrino_gains;
    double *muon_decay_muon_neutrino_gains;
    double *muon_decay_muon_antineutrino_gains;

    double *muon_decay_positive_left_muon_losses;
    double *muon_decay_positive_right_muon_losses;
    double *muon_decay_negative_left_muon_losses;
    double *muon_decay_negative_right_muon_losses;

    double *muon_decay_LUT_lifetime;
    double *muon_decay_LUT_positive_electron_minus_1;
    double *muon_decay_LUT_positive_electron_plus_1;
    double *muon_decay_LUT_positive_muon_minus_1;
    double *muon_decay_LUT_positive_muon_plus_1;

    acceleration_t electron_acceleration;
    acceleration_t proton_acceleration;

    escape_t photon_escape;
    escape_t electron_escape;

    synchrotron_t electron_synchrotron;
    synchrotron_t proton_synchrotron;
    synchrotron_t positive_pion_synchrotron;
    synchrotron_t negative_pion_synchrotron;
    synchrotron_t positive_left_muon_synchrotron;
    synchrotron_t positive_right_muon_synchrotron;
    synchrotron_t negative_left_muon_synchrotron;
    synchrotron_t negative_right_muon_synchrotron;
} state_t;

void init_state_synchrotron(state_t *st, double B);

void init_state_populations(state_t *st,
        double g_electron_min, double g_electron_max, unsigned int electron_size,
        double g_proton_min, double g_proton_max, unsigned int proton_size,
        unsigned int pion_size, unsigned int muon_size,
        unsigned int photon_size, unsigned int neutrino_size);

void init_state_aux_memory(state_t *st);

void state_init_LUTs(state_t *st);

void state_init_RK_information(state_t *st);

void step(state_t *st);
void step_tentative(state_t *st, bool try_new_step);

void state_print_data_to_file(state_t *st, enum particle_type pt, char *filename);

void state_save_state_to_file(state_t *st, char *filename);
void state_load_state_from_file(state_t *st, char *filename);

#endif /* end of include guard: STATE_H */
