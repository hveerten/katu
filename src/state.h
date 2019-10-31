#ifndef STATE_H
#define STATE_H
#pragma once

#include "config.h"

#include "acceleration.h"
#include "distribution.h"
#include "escape.h"
#include "decay_and_escape.h"
#include "population.h"
#include "synchrotron.h"

#include "utils/tpool.h"

#include <stdbool.h>

typedef struct
{
    double *stage0_population;
    double *stage_delta[4];
} RK_information_t;

typedef struct state_t
{
    double  t,  t_max;
    double dt, dt_max;

    double B;
    double R;
    double density;
    double volume;

    double eta;

    // For now, we won't use this
    /*double *big_chunk_of_memory;*/

    double photon_energy;

    population_t photons;
    population_t electrons;
    population_t positrons;
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
    RK_information_t positrons_RK_information;
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

    double inverse_compton_lepton_losses_factor;

    double *inverse_compton_photon_gains;
    double *inverse_compton_photon_gains_downscattering;
    double *inverse_compton_photon_gains_upscattering;
    double *inverse_compton_photon_losses;
    double *inverse_compton_electron_losses;
    double *inverse_compton_positron_losses;

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

    double *pair_production_photon_losses;
    double *pair_production_electron_gains;
    double *pair_production_positron_gains;

    double *pair_production_LUT_photon_losses_R;
    double *pair_production_LUT_photon_losses_index_e_min;
    double *pair_production_LUT_lepton_gains_R;
    double *pair_production_LUT_lepton_gains_index_e_min;

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

    double *bethe_heitler_LUT_proton_gamma;
    double *bethe_heitler_LUT_inelasticity;
    double *bethe_heitler_LUT_reaction_rate;

    /*double *bethe_heitler_photon_losses;*/
    /*double *bethe_heitler_proton_losses;*/
    double *bethe_heitler_electron_gains;
    double *bethe_heitler_positron_gains;

    double *pion_decay_positive_left_muon_gains;
    double *pion_decay_positive_right_muon_gains;
    double *pion_decay_negative_left_muon_gains;
    double *pion_decay_negative_right_muon_gains;

    double *pion_decay_muon_neutrino_gains;
    double *pion_decay_muon_antineutrino_gains;

    double *pion_decay_LUT_positive_right;
    double *pion_decay_LUT_positive_left;

    double *pion_decay_photon_gains;

    double *muon_decay_electron_neutrino_gains;
    double *muon_decay_electron_antineutrino_gains;
    double *muon_decay_muon_neutrino_gains;
    double *muon_decay_muon_antineutrino_gains;

    double *muon_decay_LUT_positive_electron_minus_1;
    double *muon_decay_LUT_positive_electron_plus_1;
    double *muon_decay_LUT_positive_muon_minus_1;
    double *muon_decay_LUT_positive_muon_plus_1;

    acceleration_t electron_acceleration;
    acceleration_t positron_acceleration;
    acceleration_t proton_acceleration;

    escape_t photon_escape;
    escape_t electron_escape;
    escape_t positron_escape;
    escape_t proton_escape;
    decay_and_escape_t neutron_decay_and_escape;

    decay_and_escape_t neutral_pion_decay_and_escape;
    decay_and_escape_t positive_pion_decay_and_escape;
    decay_and_escape_t negative_pion_decay_and_escape;

    decay_and_escape_t positive_left_muon_decay_and_escape;
    decay_and_escape_t positive_right_muon_decay_and_escape;
    decay_and_escape_t negative_left_muon_decay_and_escape;
    decay_and_escape_t negative_right_muon_decay_and_escape;

    escape_t electron_neutrino_escape;
    escape_t electron_antineutrino_escape;
    escape_t muon_neutrino_escape;
    escape_t muon_antineutrino_escape;

    synchrotron_t electron_synchrotron;
    synchrotron_t positron_synchrotron;
    synchrotron_t proton_synchrotron;
    synchrotron_t positive_pion_synchrotron;
    synchrotron_t negative_pion_synchrotron;
    synchrotron_t positive_left_muon_synchrotron;
    synchrotron_t positive_right_muon_synchrotron;
    synchrotron_t negative_left_muon_synchrotron;
    synchrotron_t negative_right_muon_synchrotron;

    struct
    {
        double *photons;
        double *electrons;
        double *positrons;
        double *protons;
        double *neutrons;

        double *neutral_pions;
        double *positive_pions;
        double *negative_pions;

        double *positive_left_muons;
        double *positive_right_muons;
        double *negative_left_muons;
        double *negative_right_muons;

        double *electron_neutrinos;
        double *electron_antineutrinos;
        double *muon_neutrinos;
        double *muon_antineutrinos;
    }  external_injection;

#ifdef USE_THREAD_POOL
    thread_pool_t thread_pool;
#endif

    void (*step_function)(struct state_t *st);
    void (*tentative_step_function)(struct state_t *st, bool try_new_step);
    void (*update_function)(struct state_t *st, double dt);
} state_t;

void state_init_from_config(state_t *st, config_t *cfg);

void init_state_synchrotron(state_t *st, double B);
void init_state_escape(state_t *st, double t, double cfe_ratio);
void init_state_decay_and_escape(state_t *st, double t);

void init_state_populations(state_t *st,
        double g_lepton_min, double g_lepton_max, unsigned int lepton_size,
        double g_proton_min, double g_proton_max, unsigned int proton_size,
        double e_proton_min, double e_proton_max, unsigned int photon_size,
        unsigned int pion_size, unsigned int muon_size,
        unsigned int neutrino_size);

void init_state_aux_memory(state_t *st);

void state_init_LUTs(state_t *st);

void state_init_RK_information(state_t *st);

void step(state_t *st);
void step_tentative(state_t *st, bool try_new_step);

void state_report_general_info(state_t *st);
void state_report_injection_info(state_t *st, config_t *cfg);

void state_print_data_to_file(state_t *st, enum particle_type pt, char *filename);

void state_save_state_to_file(state_t *st, char *filename);
void state_load_state_from_file(state_t *st, char *filename);

#endif /* end of include guard: STATE_H */
