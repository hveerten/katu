#include "state_step_common.h"

#include "acceleration.h"
#include "distribution.h"
#include "synchrotron.h"
#include "inverse_compton.h"
#include "pion_production.h"
#include "pion_decay.h"
#include "pair_production.h"
#include "muon_decay.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>

void step_calculate_processes(state_t *st)
{
#if USE_THREADS
    pthread_t electron_synchrotron_thread;
    pthread_t proton_synchrotron_thread;
    pthread_t positive_pion_synchrotron_thread;
    pthread_t negative_pion_synchrotron_thread;

    pthread_t positive_left_muon_synchrotron_thread;
    pthread_t positive_right_muon_synchrotron_thread;
    pthread_t negative_left_muon_synchrotron_thread;
    pthread_t negative_right_muon_synchrotron_thread;

    pthread_t inverse_compton_upscattering_thread;
    pthread_t inverse_compton_downscattering_thread;
    pthread_t inverse_compton_photon_losses_thread;
    pthread_t inverse_compton_electron_losses_thread;
    pthread_t muon_decay_thread;
    pthread_t pair_production_photon_losses_thread;
    pthread_t charged_pion_decay_thread;
    pthread_t multi_resonance_pion_production_thread;
    pthread_t multi_resonance_pion_production_hadron_losses_thread;
    pthread_t multi_resonance_pion_production_hadron_gains_thread;
    pthread_t direct_pion_production_thread;
    pthread_t direct_pion_production_hadron_losses_thread;
    pthread_t direct_pion_production_hadron_gains_thread;

    pthread_create(&electron_synchrotron_thread, NULL, st->electron_synchrotron.synchrotron_function, &st->electron_synchrotron);
    pthread_create(&proton_synchrotron_thread, NULL, st->proton_synchrotron.synchrotron_function, &st->proton_synchrotron);
    pthread_create(&positive_pion_synchrotron_thread, NULL, st->positive_pion_synchrotron.synchrotron_function, &st->positive_pion_synchrotron);
    pthread_create(&negative_pion_synchrotron_thread, NULL, st->negative_pion_synchrotron.synchrotron_function, &st->negative_pion_synchrotron);
    pthread_create(&positive_left_muon_synchrotron_thread, NULL,  st->positive_left_muon_synchrotron.synchrotron_function,  &st->positive_left_muon_synchrotron);
    pthread_create(&positive_right_muon_synchrotron_thread, NULL, st->positive_right_muon_synchrotron.synchrotron_function, &st->positive_right_muon_synchrotron);
    pthread_create(&negative_left_muon_synchrotron_thread, NULL,  st->negative_left_muon_synchrotron.synchrotron_function,  &st->negative_left_muon_synchrotron);
    pthread_create(&negative_right_muon_synchrotron_thread, NULL, st->negative_right_muon_synchrotron.synchrotron_function, &st->negative_right_muon_synchrotron);
    pthread_create(&inverse_compton_upscattering_thread, NULL, inverse_compton_head_on_upscattering_wrapper, st);
    pthread_create(&inverse_compton_downscattering_thread, NULL, inverse_compton_head_on_downscattering_wrapper, st);
    pthread_create(&inverse_compton_photon_losses_thread, NULL, inverse_compton_photon_losses_wrapper, st);
    pthread_create(&inverse_compton_electron_losses_thread, NULL, inverse_compton_electron_losses_wrapper, st);
    pthread_create(&muon_decay_thread, NULL, muon_decay_wrapper, st);
    pthread_create(&pair_production_photon_losses_thread, NULL, pair_production_photon_losses_wrapper, st);
    pthread_create(&charged_pion_decay_thread, NULL, charged_pion_decay_wrapper, st);

    pthread_create(&multi_resonance_pion_production_thread, NULL, multi_resonance_pion_production_wrapper, st);
    pthread_create(&multi_resonance_pion_production_hadron_gains_thread, NULL,  multi_resonance_pion_production_hadron_gains_wrapper, st);
    pthread_create(&multi_resonance_pion_production_hadron_losses_thread, NULL, multi_resonance_pion_production_hadron_losses_wrapper, st);
    pthread_create(&direct_pion_production_thread, NULL, direct_pion_production_wrapper, st);
    pthread_create(&direct_pion_production_hadron_gains_thread, NULL,  direct_pion_production_hadron_gains_wrapper, st);
    pthread_create(&direct_pion_production_hadron_losses_thread, NULL, direct_pion_production_hadron_losses_wrapper, st);

    pthread_join(electron_synchrotron_thread, NULL);
    pthread_join(proton_synchrotron_thread, NULL);
    pthread_join(positive_pion_synchrotron_thread, NULL);
    pthread_join(negative_pion_synchrotron_thread, NULL);
    pthread_join(positive_left_muon_synchrotron_thread, NULL);
    pthread_join(positive_right_muon_synchrotron_thread, NULL);
    pthread_join(negative_left_muon_synchrotron_thread, NULL);
    pthread_join(negative_right_muon_synchrotron_thread, NULL);
    pthread_join(inverse_compton_upscattering_thread, NULL);
    pthread_join(inverse_compton_downscattering_thread, NULL);
    pthread_join(inverse_compton_photon_losses_thread, NULL);
    pthread_join(inverse_compton_electron_losses_thread, NULL);
    pthread_join(muon_decay_thread, NULL);
    pthread_join(pair_production_photon_losses_thread, NULL);
    pthread_join(charged_pion_decay_thread, NULL);

    pthread_join(multi_resonance_pion_production_thread, NULL);
    pthread_join(multi_resonance_pion_production_hadron_gains_thread, NULL);
    pthread_join(multi_resonance_pion_production_hadron_losses_thread, NULL);
    pthread_join(direct_pion_production_thread, NULL);
    pthread_join(direct_pion_production_hadron_gains_thread, NULL);
    pthread_join(direct_pion_production_hadron_losses_thread, NULL);

    st->electron_acceleration.acceleration_function(&st->electron_acceleration);
    /*st->proton_acceleration.acceleration_function(&st->proton_acceleration);*/
#else
    st->electron_synchrotron.synchrotron_function(&st->electron_synchrotron);
    st->proton_synchrotron.synchrotron_function(&st->proton_synchrotron);
    st->positive_pion_synchrotron.synchrotron_function(&st->positive_pion_synchrotron);
    st->negative_pion_synchrotron.synchrotron_function(&st->negative_pion_synchrotron);
    st->positive_left_muon_synchrotron.synchrotron_function(&st->positive_left_muon_synchrotron);
    st->positive_right_muon_synchrotron.synchrotron_function(&st->positive_right_muon_synchrotron);
    st->negative_left_muon_synchrotron.synchrotron_function(&st->negative_left_muon_synchrotron);
    st->negative_right_muon_synchrotron.synchrotron_function(&st->negative_right_muon_synchrotron);
    inverse_compton_process_head_on_upscattering(st);
    inverse_compton_process_head_on_downscattering(st);
    inverse_compton_process_photon_losses(st);
    inverse_compton_process_electron_losses(st);
    pair_production_process_photon_losses(st);
    multi_resonance_pion_production(st);
    charged_pion_decay(st);
    muon_decay(st);

    multi_resonance_pion_production_hadron_gains(st);
    multi_resonance_pion_production_hadron_losses(st);
    direct_pion_production_hadron_gains(st);
    direct_pion_production_hadron_losses(st);

    direct_pion_production(st);

    st->electron_acceleration.acceleration_function(&st->electron_acceleration);
    st->proton_acceleration.acceleration_function(&st->proton_acceleration);
#endif

    st->photon_escape.escape_function(&st->photon_escape);
    st->electron_escape.escape_function(&st->electron_escape);
    st->proton_escape.escape_function(&st->proton_escape);
    st->neutron_decay_and_escape.losses_function(&st->neutron_decay_and_escape);

    st->neutral_pion_decay_and_escape.losses_function(&st->neutral_pion_decay_and_escape);
    st->positive_pion_decay_and_escape.losses_function(&st->positive_pion_decay_and_escape);
    st->negative_pion_decay_and_escape.losses_function(&st->negative_pion_decay_and_escape);

    st->positive_left_muon_decay_and_escape.losses_function(&st->positive_left_muon_decay_and_escape);
    st->positive_right_muon_decay_and_escape.losses_function(&st->positive_right_muon_decay_and_escape);
    st->negative_left_muon_decay_and_escape.losses_function(&st->negative_left_muon_decay_and_escape);
    st->negative_right_muon_decay_and_escape.losses_function(&st->negative_right_muon_decay_and_escape);

    st->electron_neutrino_escape.escape_function(&st->electron_neutrino_escape);
    st->electron_antineutrino_escape.escape_function(&st->electron_antineutrino_escape);
    st->muon_neutrino_escape.escape_function(&st->muon_neutrino_escape);
    st->muon_antineutrino_escape.escape_function(&st->muon_antineutrino_escape);

    unsigned int i;
    for(i = 0; i < st->photons.size; i++)
        st->inverse_compton_photon_gains[i] =
            st->inverse_compton_photon_gains_upscattering[i] +
            st->inverse_compton_photon_gains_downscattering[i];
}

void step_update_populations(state_t *st, double dt)
{
    unsigned int i;

    for(i = 0; i < st->photons.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->photons.population[i]);*/

        st->photons.tentative_population[i] =
            st->photons.population[i] + dt *
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
             st->pair_production_losses[i] +
             st->photon_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->photons.tentative_population[i],
         *        st->electron_synchrotron.photon_gains[i],
         *        st->electron_synchrotron.photon_losses[i],
         *        st->proton_synchrotron.photon_gains[i],
         *        st->proton_synchrotron.photon_losses[i],
         *        st->positive_pion_synchrotron.photon_gains[i],
         *        st->negative_pion_synchrotron.photon_losses[i],
         *        st->inverse_compton_gains[i],
         *        st->inverse_compton_loses[i],
         *        st->pair_production_losses[i],
         *        st->photon_escape.losses[i]);
         */
    }

    for(i = 0; i < st->protons.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->protons.population[i]);*/

        st->protons.tentative_population[i] =
            st->protons.population[i] + dt *
            (st->proton_synchrotron.particle_losses[i] +
             st->proton_acceleration.gains[i] +
             st->multi_resonances_proton_gains[i] +
             st->multi_resonances_proton_losses[i] +
             st->direct_pion_production_proton_gains[i] +
             st->direct_pion_production_proton_losses[i] +
             st->proton_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->protons.tentative_population[i],
         *        st->proton_synchrotron.particle_losses[i],
         *        st->proton_acceleration.gains[i],
         *        st->multi_resonances_proton_gains[i],
         *        st->multi_resonances_proton_losses[i]);
         */
    }

    for(i = 0; i < st->neutrons.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->neutrons.population[i]);*/

        st->neutrons.tentative_population[i] =
            st->neutrons.population[i] + dt *
            (st->multi_resonances_neutron_gains[i] +
             st->multi_resonances_neutron_losses[i] +
             st->direct_pion_production_neutron_gains[i] +
             st->direct_pion_production_neutron_losses[i] +
             st->neutron_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->neutrons.tentative_population[i],
         *        st->multi_resonances_neutron_gains[i],
         *        st->multi_resonances_neutron_losses[i],
         *        st->direct_pion_production_neutron_gains[i],
         *        st->direct_pion_production_neutron_losses[i]);
         */
    }

    for(i = 0; i < st->electrons.size; i++)
    {
        /*fprintf(stderr,"%i:\t%lg\t->", i, st->electrons.population[i]);*/

        st->electrons.tentative_population[i] =
            st->electrons.population[i] + dt *
            (st->electron_synchrotron.particle_losses[i] +
             st->inverse_compton_electron_losses[i] +
             st->electron_acceleration.gains[i] +
             st->electron_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\n",
         *        st->electrons.tentative_population[i],
         *        st->electron_synchrotron.particle_losses[i],
         *        st->electron_acceleration.gains[i]);
         */
    }

    for(i = 0; i < st->neutral_pions.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->neutral_pions.population[i]);*/

        st->neutral_pions.tentative_population[i] =
            st->neutral_pions.population[i] + dt *
            (st->multi_resonances_neutral_pion_gains[i] +
             st->direct_neutral_pion_gains[i] +
             st->neutral_pion_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->neutral_pions.tentative_population[i],
         *        st->multi_resonances_neutral_pion_gains[i]);
         */
    }

    for(i = 0; i < st->positive_pions.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->positive_pions.population[i]);*/

        st->positive_pions.tentative_population[i] =
            st->positive_pions.population[i] + dt *
            (st->multi_resonances_positive_pion_gains[i] +
             st->direct_positive_pion_gains[i] +
             st->positive_pion_synchrotron.particle_losses[i] +
             st->positive_pion_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->positive_pions.tentative_population[i],
         *        st->multi_resonances_positive_pion_gains[i],
         *        st->direct_positive_pion_gains[i],
         *        st->positive_pion_synchrotron.particle_losses[i],
         *        st->pion_decay_positive_pion_losses[i]);
         */
    }

    for(i = 0; i < st->negative_pions.size; i++)
    {
        st->negative_pions.tentative_population[i] =
            st->negative_pions.population[i] + dt *
            (st->multi_resonances_negative_pion_gains[i] +
             st->direct_negative_pion_gains[i] +
             st->negative_pion_synchrotron.particle_losses[i] +
             st->negative_pion_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_left_muons.size; i++)
    {
        st->positive_left_muons.tentative_population[i] =
            st->positive_left_muons.population[i] + dt *
            (st->pion_decay_positive_left_muon_gains[i] +
             st->positive_left_muon_synchrotron.particle_losses[i] +
             st->positive_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_right_muons.size; i++)
    {
        st->positive_right_muons.tentative_population[i] =
            st->positive_right_muons.population[i] + dt *
            (st->pion_decay_positive_right_muon_gains[i] +
             st->positive_right_muon_synchrotron.particle_losses[i] +
             st->positive_right_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_left_muons.size; i++)
    {
        st->negative_left_muons.tentative_population[i] =
            st->negative_left_muons.population[i] + dt *
            (st->pion_decay_negative_left_muon_gains[i] +
             st->negative_left_muon_synchrotron.particle_losses[i] +
             st->negative_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_right_muons.size; i++)
    {
        st->negative_right_muons.tentative_population[i] =
            st->negative_right_muons.population[i] + dt *
            (st->pion_decay_negative_right_muon_gains[i] +
             st->negative_right_muon_synchrotron.particle_losses[i] +
             st->negative_right_muon_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"%u:\t%lg\t%lg\t%lg\n", i,
         *        st->pion_decay_negative_right_muon_gains[i],
         *        st->negative_right_muon_synchrotron.particle_losses[i],
         *        st->muon_decay_negative_right_muon_losses[i]);
         */
    }

    for(i = 0; i < st->electron_neutrinos.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->electron_neutrinos.population[i]);*/

        st->electron_neutrinos.tentative_population[i] =
            st->electron_neutrinos.population[i] + dt *
            (st->muon_decay_electron_neutrino_gains[i] +
             st->electron_neutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->electron_neutrinos.tentative_population[i],
         *        st->muon_decay_electron_neutrino_gains[i]);
         */
    }

    for(i = 0; i < st->electron_antineutrinos.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->electron_antineutrinos.population[i]);*/

        st->electron_antineutrinos.tentative_population[i] =
            st->electron_antineutrinos.population[i] + dt *
            (st->muon_decay_electron_antineutrino_gains[i] +
             st->electron_antineutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->electron_antineutrinos.tentative_population[i],
         *        st->muon_decay_electron_antineutrino_gains[i]);
         */
    }

    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        st->muon_neutrinos.tentative_population[i] =
            st->muon_neutrinos.population[i] + dt *
            (st->pion_decay_muon_neutrino_gains[i] +
             st->muon_decay_muon_neutrino_gains[i] +
             st->muon_neutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"%u:\t%lg\t%lg\n", i,
         *        st->pion_decay_muon_neutrino_gains[i],
         *        st->muon_decay_muon_neutrino_gains[i]);
         */
    }

    for(i = 0; i < st->muon_antineutrinos.size; i++)
    {
        st->muon_antineutrinos.tentative_population[i] =
            st->muon_antineutrinos.population[i] + dt *
            (st->pion_decay_muon_antineutrino_gains[i] +
             st->muon_decay_muon_antineutrino_gains[i] +
             st->muon_antineutrino_escape.losses[i]);
    }
}

void step_experimental_update_populations(state_t *st, double dt)
{
    unsigned int i;
    double g_turnover;
    double dlng;

    for(i = 0; i < st->photons.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->photons.population[i]);*/

        double phot_gains = st->electron_synchrotron.photon_gains[i] +
                            st->proton_synchrotron.photon_gains[i] +
                            st->positive_pion_synchrotron.photon_gains[i] +
                            st->negative_pion_synchrotron.photon_gains[i] +
                            st->positive_left_muon_synchrotron.photon_gains[i] +
                            st->positive_right_muon_synchrotron.photon_gains[i] +
                            st->negative_left_muon_synchrotron.photon_gains[i] +
                            st->negative_right_muon_synchrotron.photon_gains[i] +
                            st->inverse_compton_photon_gains[i];

        double phot_losses = st->electron_synchrotron.photon_losses[i] +
                             st->proton_synchrotron.photon_losses[i] +
                             st->positive_pion_synchrotron.photon_losses[i] +
                             st->negative_pion_synchrotron.photon_losses[i] +
                             st->positive_left_muon_synchrotron.photon_losses[i] +
                             st->positive_right_muon_synchrotron.photon_losses[i] +
                             st->negative_left_muon_synchrotron.photon_losses[i] +
                             st->negative_right_muon_synchrotron.photon_losses[i] +
                             st->inverse_compton_photon_losses[i] +
                             st->pair_production_losses[i];

        double L = phot_losses / st->photons.population[i] - 1 / st->photon_escape.t;
        double tau = L * dt;
        double aux = expm1(tau) / L;

        double new_pop = st->photons.population[i] * exp(tau) + phot_gains * aux;

        st->photons.tentative_population[i] = new_pop;

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->photons.tentative_population[i],
         *        st->electron_synchrotron.photon_gains[i],
         *        st->electron_synchrotron.photon_losses[i],
         *        st->proton_synchrotron.photon_gains[i],
         *        st->proton_synchrotron.photon_losses[i],
         *        st->positive_pion_synchrotron.photon_gains[i],
         *        st->negative_pion_synchrotron.photon_losses[i],
         *        st->inverse_compton_gains[i],
         *        st->inverse_compton_loses[i],
         *        st->pair_production_losses[i],
         *        st->photon_escape.losses[i]);
         */
    }

#if PROTON_STEADY_STATE == 0
    g_turnover = 1 / (st->proton_acceleration.t * st->proton_synchrotron.particle_losses_factor);
    dlng = st->protons.log_energy[1] - st->protons.log_energy[0];
    for(i = 0; i < st->protons.size && st->protons.energy[i] < g_turnover; i++)
    {
        double proton_gains =
             st->multi_resonances_proton_gains[i] +
             st->direct_pion_production_proton_gains[i];

        double proton_losses =
             st->multi_resonances_proton_losses[i] +
             st->direct_pion_production_proton_losses[i];

        double L = proton_losses / st->protons.population[i] +
                   2 * st->protons.energy[i] * st->proton_synchrotron.particle_losses_factor +
                   - 1 / st->proton_acceleration.t - 1 / st->proton_escape.t;
        double tau = L * dt;

        double aux1 = (st->proton_synchrotron.particle_losses_factor * st->protons.energy[i] - 1 / st->proton_acceleration.t) / dlng;
        double aux2 = expm1(tau) / L;

        double new_pop;
        if(i != 0)
        {
            new_pop = (st->protons.population[i] * exp(tau) + aux2 * (proton_gains - aux1 * st->protons.tentative_population[i - 1])) /
                      (1 - aux1 * aux2);
        }
        else
            new_pop = st->protons.population[i];

        st->protons.tentative_population[i] = new_pop;
    }
    for(i = st->protons.size - 1; i < st->protons.size && g_turnover < st->protons.energy[i]; i--)
    {
        double proton_gains =
             st->multi_resonances_proton_gains[i] +
             st->direct_pion_production_proton_gains[i];

        double proton_losses =
             st->multi_resonances_proton_losses[i] +
             st->direct_pion_production_proton_losses[i];

        double L = proton_losses / st->protons.population[i] + 2 * st->protons.energy[i] * st->proton_synchrotron.particle_losses_factor - 1 / st->proton_acceleration.t;
        double tau = L * dt;

        double aux1 = (st->proton_synchrotron.particle_losses_factor * st->protons.energy[i] - 1 / st->proton_acceleration.t) / dlng;
        double aux2 = expm1(tau) / L;

        double new_pop;
        if(i != st->protons.size - 1)
        {
            new_pop = (st->protons.population[i] * exp(tau) + aux2 * (proton_gains + aux1 * st->protons.tentative_population[i + 1])) /
                      (1 + aux1 * aux2);
        }
        else
            new_pop = st->protons.population[i];

        st->protons.tentative_population[i] = new_pop;
    }
#endif

    for(i = 0; i < st->neutrons.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->neutrons.population[i]);*/

        st->neutrons.tentative_population[i] =
            st->neutrons.population[i] + dt *
            (st->multi_resonances_neutron_gains[i] +
             st->multi_resonances_neutron_losses[i] +
             st->direct_pion_production_neutron_gains[i] +
             st->direct_pion_production_neutron_losses[i] +
             st->neutron_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->neutrons.tentative_population[i],
         *        st->multi_resonances_neutron_gains[i],
         *        st->multi_resonances_neutron_losses[i],
         *        st->direct_pion_production_neutron_gains[i],
         *        st->direct_pion_production_neutron_losses[i]);
         */
    }

#if ELECTRON_STEADY_STATE == 0
    double electron_aux0[st->electrons.size];
    double electron_aux1[st->electrons.size];
    for(i = 0; i < st->electrons.size; i++)
    {
        double g     = st->electrons.energy[i];
        double S     = st->electron_synchrotron.particle_losses_factor;
        double IC    = st->inverse_compton_electron_losses_factor;
        double t_acc = st->electron_acceleration.t;

        electron_aux0[i] = 2 * g * (S + IC) - 1 / t_acc;
        electron_aux1[i] =     g * (S + IC) - 1 / t_acc;
    }

    g_turnover = 1 / (st->electron_acceleration.t * (st->electron_synchrotron.particle_losses_factor + st->inverse_compton_electron_losses_factor));
    double g_max = st->electrons.energy[0] * exp(st->t / st->electron_acceleration.t) / (st->electrons.energy[0] / g_turnover * expm1(st->t / st->electron_acceleration.t) + 1);
    dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

#if 0
    st->electrons.tentative_population[0] = st->electrons.population[0];
    for(i = 1; i < st->electrons.size - 1 && st->electrons.energy[i] < g_turnover; i++)
    /*for(i = 1; i < st->electrons.size - 1 && st->electrons.energy[i] < g_turnover && st->electrons.energy[i] < g_max; i++)*/
    {
        /*double electron_gains  = 1e-2 * pow(st->electrons.energy[i], -2.3);*/
        double electron_gains  = 0;
        double electron_losses = 0;

        double L = electron_losses + electron_aux0[i] - 1/ st->electron_escape.t;
        double tau = L * dt;

        double aux1 = electron_aux1[i] / dlng;
        double aux2 = expm1(tau) / L;

        double new_pop = (st->electrons.population[i] * exp(tau) + aux2 * (electron_gains - aux1 * st->electrons.tentative_population[i - 1])) /
                         (1 - aux1 * aux2);

        st->electrons.tentative_population[i] = new_pop;
    }
#endif

#if 1
    st->electrons.tentative_population[0] = st->electrons.population[0];
    {
        double electron_losses = 0;

        double aux1 = electron_aux1[1];
        double aux2 = electron_aux0[1] - 1/ st->electron_escape.t;

        double log_new_pop = log(st->electrons.tentative_population[0]);

        log_new_pop = (st->electrons.log_population[1] + st->dt * (aux2 - aux1 * log_new_pop / dlng)) /
                        (1 - aux1 * st->dt / dlng);

        st->electrons.tentative_population[1] = exp(log_new_pop);
    }
    for(i = 2; i < st->electrons.size - 1 && st->electrons.energy[i] < g_turnover; i++)
    {
        double electron_losses = 0;

        double aux1 = electron_aux1[i];
        double aux2 = electron_aux0[i] - 1/ st->electron_escape.t;

        double log_new_pop = log(st->electrons.tentative_population[i - 1]);
        double log_old_pop = log(st->electrons.tentative_population[i - 2]);

        log_new_pop = (st->electrons.log_population[i] + st->dt * (aux2 - aux1 * (4*log_new_pop - log_old_pop) / dlng/2)) /
                        (1 - 3*aux1 * st->dt / dlng/2);

        st->electrons.tentative_population[i] = exp(log_new_pop);

        /*fprintf(stderr,"%u:\t%lg\t%lg\n", i, log_new_pop, log(st->electrons.tentative_population[i - 2]));*/
    }
    st->electrons.tentative_population[0] = st->electrons.population[0];
#endif


    // ELECTRON COOLING BELOW
    /*for(i = st->electrons.size - 1; i != 0 && g_turnover < st->electrons.energy[i]; i--)*/
    if(0)
    {
        /*double electron_gains  = 1e-1 * pow(st->electrons.energy[i], -2.3);*/
        double electron_gains  = 0;
        double electron_losses = 0;

        double L = electron_losses + \
                   2 * st->electrons.energy[i] * (st->electron_synchrotron.particle_losses_factor + \
                                                  st->inverse_compton_electron_losses_factor) - \
                   1 / st->electron_acceleration.t - 1/ st->electron_escape.t;
        double tau = L * dt;

        double aux1 = ((st->electron_synchrotron.particle_losses_factor + \
                        st->inverse_compton_electron_losses_factor) * st->electrons.energy[i] - 1 / st->electron_acceleration.t) / dlng;
        double aux2 = expm1(tau) / L;

        double new_pop;
        if(i != st->electrons.size - 1)
        {
            new_pop = (st->electrons.population[i] * exp(tau) + aux2 * (electron_gains + aux1 * st->electrons.tentative_population[i + 1])) /
                      (1 + aux1 * aux2);
        }
        else
            new_pop = st->electrons.population[i];

        st->electrons.tentative_population[i] = new_pop;
    }
    /*for(i = st->electrons.size - 1; i != 0 && g_turnover < st->electrons.energy[i]; i--)*/
    if(0)
    {
        double electron_gains  = 0;
        double electron_losses = 0;

        double aux1 = ((st->electron_synchrotron.particle_losses_factor + \
                        st->inverse_compton_electron_losses_factor) * st->electrons.energy[i] - 1 / st->electron_acceleration.t);
        double aux2 = 2 * st->electrons.energy[i] * (st->electron_synchrotron.particle_losses_factor + \
                                                     st->inverse_compton_electron_losses_factor) - \
                        1 / st->electron_acceleration.t - 1/ st->electron_escape.t;
        double L = electron_losses + aux2 + aux1 * (st->electrons.log_population[i + 1] - st->electrons.log_population[i]) / dlng;
        double tau = L * dt;

        double new_pop;
        if(i != st->electrons.size - 1)
        {
            new_pop = st->electrons.population[i] * exp(tau);
        }
        else
            new_pop = st->electrons.population[i];

        st->electrons.tentative_population[i] = new_pop;
    }
    for(i = st->electrons.size - 1; i != 0 && g_turnover < st->electrons.energy[i]; i--)
    /*if(0)*/
    {
        double electron_losses = 0;

        double aux1 = ((st->electron_synchrotron.particle_losses_factor + \
                        st->inverse_compton_electron_losses_factor) * st->electrons.energy[i] - 1 / st->electron_acceleration.t);
        double aux2 = 2 * st->electrons.energy[i] * (st->electron_synchrotron.particle_losses_factor + \
                                                     st->inverse_compton_electron_losses_factor) - \
                        1 / st->electron_acceleration.t - 1/ st->electron_escape.t;

        double log_new_pop = log(st->electrons.tentative_population[i + 1]);
        if(i != st->electrons.size - 1)
        {
            log_new_pop = (st->electrons.log_population[i] + st->dt * (aux2 + aux1 * log_new_pop / dlng)) /
                          (1 + aux1 * st->dt / dlng);
        }
        else
            log_new_pop = st->electrons.log_population[i];

        st->electrons.tentative_population[i] = exp(log_new_pop);
    }
#endif

    for(i = 0; i < st->neutral_pions.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->neutral_pions.population[i]);*/

        st->neutral_pions.tentative_population[i] =
            st->neutral_pions.population[i] + dt *
            (st->multi_resonances_neutral_pion_gains[i] +
             st->direct_neutral_pion_gains[i] +
             st->neutral_pion_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->neutral_pions.tentative_population[i],
         *        st->multi_resonances_neutral_pion_gains[i]);
         */
    }

    for(i = 0; i < st->positive_pions.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->positive_pions.population[i]);*/

        st->positive_pions.tentative_population[i] =
            st->positive_pions.population[i] + dt *
            (st->multi_resonances_positive_pion_gains[i] +
             st->direct_positive_pion_gains[i] +
             st->positive_pion_synchrotron.particle_losses[i] +
             st->positive_pion_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\t%lg\t%lg\t%lg\n",
         *        st->positive_pions.tentative_population[i],
         *        st->multi_resonances_positive_pion_gains[i],
         *        st->direct_positive_pion_gains[i],
         *        st->positive_pion_synchrotron.particle_losses[i],
         *        st->pion_decay_positive_pion_losses[i]);
         */
    }

    for(i = 0; i < st->negative_pions.size; i++)
    {
        st->negative_pions.tentative_population[i] =
            st->negative_pions.population[i] + dt *
            (st->multi_resonances_negative_pion_gains[i] +
             st->direct_negative_pion_gains[i] +
             st->negative_pion_synchrotron.particle_losses[i] +
             st->negative_pion_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_left_muons.size; i++)
    {
        st->positive_left_muons.tentative_population[i] =
            st->positive_left_muons.population[i] + dt *
            (st->pion_decay_positive_left_muon_gains[i] +
             st->positive_left_muon_synchrotron.particle_losses[i] +
             st->positive_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->positive_right_muons.size; i++)
    {
        st->positive_right_muons.tentative_population[i] =
            st->positive_right_muons.population[i] + dt *
            (st->pion_decay_positive_right_muon_gains[i] +
             st->positive_right_muon_synchrotron.particle_losses[i] +
             st->positive_right_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_left_muons.size; i++)
    {
        st->negative_left_muons.tentative_population[i] =
            st->negative_left_muons.population[i] + dt *
            (st->pion_decay_negative_left_muon_gains[i] +
             st->negative_left_muon_synchrotron.particle_losses[i] +
             st->negative_left_muon_decay_and_escape.losses[i]);
    }

    for(i = 0; i < st->negative_right_muons.size; i++)
    {
        st->negative_right_muons.tentative_population[i] =
            st->negative_right_muons.population[i] + dt *
            (st->pion_decay_negative_right_muon_gains[i] +
             st->negative_right_muon_synchrotron.particle_losses[i] +
             st->negative_right_muon_decay_and_escape.losses[i]);

        /*
         *fprintf(stderr,"%u:\t%lg\t%lg\t%lg\n", i,
         *        st->pion_decay_negative_right_muon_gains[i],
         *        st->negative_right_muon_synchrotron.particle_losses[i],
         *        st->muon_decay_negative_right_muon_losses[i]);
         */
    }

    for(i = 0; i < st->electron_neutrinos.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->electron_neutrinos.population[i]);*/

        st->electron_neutrinos.tentative_population[i] =
            st->electron_neutrinos.population[i] + dt *
            (st->muon_decay_electron_neutrino_gains[i] +
             st->electron_neutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->electron_neutrinos.tentative_population[i],
         *        st->muon_decay_electron_neutrino_gains[i]);
         */
    }

    for(i = 0; i < st->electron_antineutrinos.size; i++)
    {
        /*fprintf(stderr,"%u:\t%lg\t->", i, st->electron_antineutrinos.population[i]);*/

        st->electron_antineutrinos.tentative_population[i] =
            st->electron_antineutrinos.population[i] + dt *
            (st->muon_decay_electron_antineutrino_gains[i] +
             st->electron_antineutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"\t%lg\t%lg\n",
         *        st->electron_antineutrinos.tentative_population[i],
         *        st->muon_decay_electron_antineutrino_gains[i]);
         */
    }

    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        st->muon_neutrinos.tentative_population[i] =
            st->muon_neutrinos.population[i] + dt *
            (st->pion_decay_muon_neutrino_gains[i] +
             st->muon_decay_muon_neutrino_gains[i] +
             st->muon_neutrino_escape.losses[i]);

        /*
         *fprintf(stderr,"%u:\t%lg\t%lg\n", i,
         *        st->pion_decay_muon_neutrino_gains[i],
         *        st->muon_decay_muon_neutrino_gains[i]);
         */
    }

    for(i = 0; i < st->muon_antineutrinos.size; i++)
    {
        st->muon_antineutrinos.tentative_population[i] =
            st->muon_antineutrinos.population[i] + dt *
            (st->pion_decay_muon_antineutrino_gains[i] +
             st->muon_decay_muon_antineutrino_gains[i] +
             st->muon_antineutrino_escape.losses[i]);
    }
}

/* Check that the new populations are sensible.
 *
 * Hence, what we do is first check for infinities, and in those cases
 * crash.
 * Then we check for negative numbers. In case we find those, it means
 * that we have overshoot the 0 and that the timestep was to big, so
 * return that we have problems.
 */
bool step_check_tentative_populations(state_t *st)
{
    unsigned int i;
#ifndef NDEBUG
#define CHECK_POPULATIONS_BODY(X) \
    for(i = 0; i < st->X.size; i++) \
    { \
        if(!isfinite(st->X.tentative_population[i])) \
        { \
            fprintf(stderr,#X" not finite:\t%u\t%lg\n",i,st->X.tentative_population[i]); \
            assert(0);\
        } \
        if(st->X.tentative_population[i] < -DBL_MIN) \
        {\
            fprintf(stderr,#X" negative:\t%u\t%lg\n",i,st->X.tentative_population[i]); \
            return true;\
        }\
    }
#else
#define CHECK_POPULATIONS_BODY(X) \
    for(i = 0; i < st->X.size; i++) \
    { \
        if(st->X.tentative_population[i] < -DBL_MIN) \
            return true;\
    }
#endif
#define STEP_INTERNAL_FUNCTION(X) CHECK_POPULATIONS_BODY(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef CHECK_POPULATIONS_BODY

    return false;
}

/* Raise 'zeros' to DBL_MIN */
void step_fix_tentative_zeros(state_t *st)
{
    unsigned int i;
#define FIX_ZEROS_BODY(X) \
    for(i = 0; i < st->X.size; i++) \
        st->X.tentative_population[i] = fmax(DBL_MIN, st->X.tentative_population[i]);
#define STEP_INTERNAL_FUNCTION(X) FIX_ZEROS_BODY(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef FIX_ZEROS_BODY
}

void step_accept_tentative_population(state_t *st)
{
    double *temp_population;
#define SWAP_POP(X) \
    temp_population = st->X.population; \
    st->X.population = st->X.tentative_population; \
    st->X.tentative_population = temp_population;
#define STEP_INTERNAL_FUNCTION(X) SWAP_POP(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef SWAP_POP
}

void step_log_populations(state_t *st)
{
    unsigned int i;

#define LOG_POPULATION(X) \
    for(i = 0; i < st->X.size; i++) \
        st->X.log_population[i] = log(st->X.population[i]);
#define STEP_INTERNAL_FUNCTION(X) LOG_POPULATION(X)
    APPLY_STEP_INTERNAL_FUNCTION_TO_PARTICLES;
#undef STEP_INTERNAL_FUNCTION
#undef LOG_POPULATION
}

void step_propagate_new_dt(state_t *st, double dt)
{
    st->dt = dt;

    update_acceleration(st, &st->proton_acceleration,   st->proton_acceleration.t);
    update_acceleration(st, &st->electron_acceleration, st->electron_acceleration.t);

    update_escape(st, &st->photon_escape,   st->photon_escape.t);
    update_escape(st, &st->electron_escape, st->electron_escape.t);
    update_escape(st, &st->proton_escape,   st->proton_escape.t);
    update_decay_and_escape(st, &st->neutron_decay_and_escape,  st->neutron_decay_and_escape.escape_lifetime);

    update_decay_and_escape(st, &st->neutral_pion_decay_and_escape,  st->neutral_pion_decay_and_escape.escape_lifetime);
    update_decay_and_escape(st, &st->positive_pion_decay_and_escape, st->positive_pion_decay_and_escape.escape_lifetime);
    update_decay_and_escape(st, &st->negative_pion_decay_and_escape, st->negative_pion_decay_and_escape.escape_lifetime);

    update_decay_and_escape(st, &st->positive_left_muon_decay_and_escape,  st->positive_left_muon_decay_and_escape.escape_lifetime);
    update_decay_and_escape(st, &st->positive_right_muon_decay_and_escape, st->positive_right_muon_decay_and_escape.escape_lifetime);
    update_decay_and_escape(st, &st->negative_left_muon_decay_and_escape,  st->negative_left_muon_decay_and_escape.escape_lifetime);
    update_decay_and_escape(st, &st->negative_right_muon_decay_and_escape, st->negative_right_muon_decay_and_escape.escape_lifetime);

    update_escape(st, &st->electron_neutrino_escape,     st->electron_neutrino_escape.t);
    update_escape(st, &st->electron_antineutrino_escape, st->electron_antineutrino_escape.t);
    update_escape(st, &st->muon_neutrino_escape,         st->muon_neutrino_escape.t);
    update_escape(st, &st->muon_antineutrino_escape,     st->muon_antineutrino_escape.t);

    synchrotron_update_dt(&st->electron_synchrotron, dt);
    synchrotron_update_dt(&st->proton_synchrotron, dt);
    synchrotron_update_dt(&st->positive_pion_synchrotron, dt);
    synchrotron_update_dt(&st->negative_pion_synchrotron, dt);
    synchrotron_update_dt(&st->negative_left_muon_synchrotron, dt);
    synchrotron_update_dt(&st->negative_right_muon_synchrotron, dt);
    synchrotron_update_dt(&st->positive_left_muon_synchrotron, dt);
    synchrotron_update_dt(&st->positive_right_muon_synchrotron, dt);
}
