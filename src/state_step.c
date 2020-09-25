#include "state.h"
#include "constants.h"

#include "state_step_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <float.h>
#include <stdbool.h>

#include <pthread.h>
#include <assert.h>

void step(state_t *st)
{
    step_calculate_processes(st);
    st->update_function(st, st->dt);

    if(step_check_tentative_populations(st))
        assert(0);
    step_fix_tentative_zeros(st);
    step_accept_tentative_population(st);

    step_log_populations(st);

    st->t += st->dt;
}

void step_tentative(state_t *st, bool try_new_step)
{
    double dt_old = st->dt;
    double dt_new = fmin(st->dt * 1.1, st->dt_max);

    if(try_new_step)
    {
        /*fprintf(stderr,"Testing: %lg\n", dt_new);*/
        step_propagate_new_dt(st, dt_new);
    }
    else
        dt_old = st->dt / 1.1;

    step_calculate_processes(st);
    st->update_function(st, st->dt);

    if(step_check_tentative_populations(st))
        goto abort;
    step_fix_tentative_zeros(st);

    goto end;

abort:
    /* We have aborted the new timestep.
     * So, set the old one and recalculate LUTs.
     * Finally, do a step with the old timestep */
    step_propagate_new_dt(st, dt_old);
    /*fprintf(stderr,"Back To: %lg\n", dt_old);*/

    step_tentative(st, 0);
    return;

end:;
    step_accept_tentative_population(st);
    step_log_populations(st);

    st->t += st->dt;
}

static void step_accept_tentative_population_RK(state_t *st)
{
    double *temp_population;
#define SWAP_POP(X) \
    temp_population = st->X.population; \
    st->X.population = st->X.tentative_population; \
    st->X.tentative_population = st->X##_RK_information.stage0_population; \
    st->X##_RK_information.stage0_population = temp_population;
#define STEP_INTERNAL_FUNCTION(X) SWAP_POP(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef SWAP_POP
}

static void step_calculate_deltas_RK(state_t *st, unsigned int stage)
{
    unsigned int i;

    for(i = 0; i < st->photons.size; i++)
    {
        st->photons_RK_information.stage_delta[stage][i] =
            (st->electron_synchrotron.photon_gains[i] +
             st->electron_synchrotron.photon_losses[i] +
             st->proton_synchrotron.photon_gains[i] +
             st->proton_synchrotron.photon_losses[i] +
             st->positive_pion_synchrotron.photon_gains[i] +
             st->positive_pion_synchrotron.photon_losses[i] +
             st->negative_pion_synchrotron.photon_gains[i] +
             st->negative_pion_synchrotron.photon_losses[i] +
             st->positive_left_muon_synchrotron.photon_gains[i] +
             st->positive_left_muon_synchrotron.photon_losses[i] +
             st->positive_right_muon_synchrotron.photon_gains[i] +
             st->positive_right_muon_synchrotron.photon_losses[i] +
             st->negative_left_muon_synchrotron.photon_gains[i] +
             st->negative_left_muon_synchrotron.photon_losses[i] +
             st->negative_right_muon_synchrotron.photon_gains[i] +
             st->negative_right_muon_synchrotron.photon_losses[i] +
             st->inverse_compton_photon_gains[i] + st->inverse_compton_photon_losses[i] +
             st->pair_production_photon_losses[i] +
             st->photon_escape.losses[i]);
    }

    for(i = 0; i < st->protons.size; i++)
    {
        st->protons_RK_information.stage_delta[stage][i] =
            (st->proton_synchrotron.particle_losses[i] +
             st->proton_acceleration.gains[i] +
             st->multi_resonances_proton_gains[i] +
             st->multi_resonances_proton_losses[i] +
             st->direct_pion_production_proton_gains[i] +
             st->direct_pion_production_proton_losses[i] +
             st->proton_escape.losses[i]);
    }

    for(i = 0; i < st->neutrons.size; i++)
    {
        st->neutrons_RK_information.stage_delta[stage][i] =
            (st->multi_resonances_neutron_gains[i] +
             st->multi_resonances_neutron_losses[i] +
             st->direct_pion_production_neutron_gains[i] +
             st->direct_pion_production_neutron_losses[i] +
             st->neutron_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->electrons.size; i++)
    {
        st->electrons_RK_information.stage_delta[stage][i] =
            (st->electron_synchrotron.particle_losses[i] +
             st->inverse_compton_electron_losses[i] +
             st->electron_acceleration.gains[i] +
             st->electron_escape.losses[i]);
    }

    for(i = 0; i < st->neutral_pions.size; i++)
    {
        st->neutral_pions_RK_information.stage_delta[stage][i] =
            (st->multi_resonances_neutral_pion_gains[i] +
             st->direct_neutral_pion_gains[i] +
             st->neutral_pion_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_pions.size; i++)
    {
        st->positive_pions_RK_information.stage_delta[stage][i] =
            (st->multi_resonances_positive_pion_gains[i] +
             st->direct_positive_pion_gains[i] +
             st->positive_pion_synchrotron.particle_losses[i] +
             st->positive_pion_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_pions.size; i++)
    {
        st->negative_pions_RK_information.stage_delta[stage][i] =
            (st->multi_resonances_negative_pion_gains[i] +
             st->direct_negative_pion_gains[i] +
             st->negative_pion_synchrotron.particle_losses[i] +
             st->negative_pion_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_left_muons.size; i++)
    {
        st->positive_left_muons_RK_information.stage_delta[stage][i] =
            (st->pion_decay_positive_left_muon_gains[i] +
             st->positive_left_muon_synchrotron.particle_losses[i] +
             st->positive_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_right_muons.size; i++)
    {
        st->positive_right_muons_RK_information.stage_delta[stage][i] =
            (st->pion_decay_positive_right_muon_gains[i] +
             st->positive_right_muon_synchrotron.particle_losses[i] +
             st->positive_right_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_left_muons.size; i++)
    {
        st->negative_left_muons_RK_information.stage_delta[stage][i] =
            (st->pion_decay_negative_left_muon_gains[i] +
             st->negative_left_muon_synchrotron.particle_losses[i] +
             st->negative_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_right_muons.size; i++)
    {
        st->negative_right_muons_RK_information.stage_delta[stage][i] =
            (st->pion_decay_negative_right_muon_gains[i] +
             st->negative_right_muon_synchrotron.particle_losses[i] +
             st->negative_right_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->electron_neutrinos.size; i++)
    {
        st->electron_neutrinos_RK_information.stage_delta[stage][i] =
            (st->muon_decay_electron_neutrino_gains[i] +
             st->electron_neutrino_escape.losses[i]);
    }

    for(i = 0; i < st->electron_antineutrinos.size; i++)
    {
        st->electron_antineutrinos_RK_information.stage_delta[stage][i] =
            (st->muon_decay_electron_antineutrino_gains[i] +
             st->electron_antineutrino_escape.losses[i]);
    }

    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        st->muon_neutrinos_RK_information.stage_delta[stage][i] =
            (st->pion_decay_muon_neutrino_gains[i] +
             st->muon_decay_muon_neutrino_gains[i] +
             st->muon_neutrino_escape.losses[i]);
    }

    for(i = 0; i < st->muon_antineutrinos.size; i++)
    {
        st->muon_antineutrinos_RK_information.stage_delta[stage][i] =
            (st->pion_decay_muon_antineutrino_gains[i] +
             st->muon_decay_muon_antineutrino_gains[i] +
             st->muon_antineutrino_escape.losses[i]);
    }
}

static void step_update_populations_heun(state_t *st, double dt)
{
    unsigned int i;

#define UPDATE_POPULATIONS_HEUN(X) \
    for(i = 0; i < st->X.size; i++) \
        st->X.tentative_population[i] = \
            st->X##_RK_information.stage0_population[i] + dt * \
                (st->X##_RK_information.stage_delta[0][i] + \
                 st->X##_RK_information.stage_delta[1][i]) / 2;
#define STEP_INTERNAL_FUNCTION(X) UPDATE_POPULATIONS_HEUN(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef UPDATE_POPULATIONS_HEUN
}

void step_heun(state_t *st)
{
    // First step
    step_calculate_processes(st);
    step_calculate_deltas_RK(st, 0);
    step_update_populations(st, st->dt);
    step_accept_tentative_population_RK(st);
    step_log_populations(st);

    // Second step
    step_calculate_processes(st);
    step_calculate_deltas_RK(st, 1);

    // Combination step
    step_update_populations_heun(st, st->dt);
    step_accept_tentative_population(st);
    step_log_populations(st);

    st->t += st->dt;
}

void step_heun_tentative(state_t *st, bool try_new_step)
{
    double dt_old = st->dt;
    double dt_new = fmin(st->dt * 1.1, st->dt_max);

    if(try_new_step)
    {
        /*fprintf(stderr,"Testing: %lg\n", dt_new);*/
        step_propagate_new_dt(st, dt_new);
    }
    else
        dt_old = st->dt / 1.1;

    // First step
    step_calculate_processes(st);
    st->update_function(st, st->dt);

    if(step_check_tentative_populations(st))
        goto abort1;
    step_fix_tentative_zeros(st);

    step_calculate_deltas_RK(st, 0);
    step_accept_tentative_population_RK(st);
    step_log_populations(st);

    // Second step
    step_calculate_processes(st);
    step_calculate_deltas_RK(st, 1);

    // Combination step
    step_update_populations_heun(st, st->dt);

    if(step_check_tentative_populations(st))
        goto abort;
    step_fix_tentative_zeros(st);

    step_accept_tentative_population(st);
    step_log_populations(st);

    st->t += st->dt;
    return;

abort:;
      /* We have aborted the new timestep NOT in the first stage.
       * This means that the original population is in the RK_information.
       * Retrieve it and put it in st->X.population and fall-through */

    double *temp_population;
#define SWAP_POP(X) \
    temp_population = st->X.population; \
    st->X.population = st->X##_RK_information.stage0_population; \
    st->X##_RK_information.stage0_population = temp_population;
#define STEP_INTERNAL_FUNCTION(X) SWAP_POP(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef SWAP_POP

abort1:
    /* We have aborted the new timestep in the first stage.
     * This means that the original population is still in
     * st->X.population, so we don't need to recover it from
     * the RK_information.
     * So, set the old time step and recalculate LUTs.
     * Finally, do a step with the old timestep */
    step_propagate_new_dt(st, dt_old);
    /*fprintf(stderr,"Back To: %lg\n", dt_old);*/

    step_heun_tentative(st, 0);
    return;
}
