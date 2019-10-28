#include "state.h"
#include "constants.h"
#include "config.h"
#include "acceleration.h"
#include "distribution.h"
#include "synchrotron.h"
#include "inverse_compton.h"
#include "pion_production.h"
#include "pion_decay.h"
#include "pair_production.h"
#include "muon_decay.h"
#include "bethe_heitler.h"

#include "state_step.h"
#include "state_step_common.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <assert.h>

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

static void state_init_injection(state_t *st, config_t *cfg)
{
    unsigned int i;

    population_t electron_aux_pop;
    population_t proton_aux_pop;

    init_population(&electron_aux_pop, electron, cfg->ei.electron_distribution.min, cfg->ei.electron_distribution.max, st->electrons.size);
    init_population(&proton_aux_pop,   proton,   cfg->ei.proton_distribution.min,   cfg->ei.proton_distribution.max,   st->protons.size);

    generate_population(&electron_aux_pop, &cfg->ei.electron_distribution);
    generate_population(&proton_aux_pop,   &cfg->ei.proton_distribution);

    double electron_energy_average = distribution_average(&cfg->ei.electron_distribution);
    double proton_energy_average   = distribution_average(&cfg->ei.proton_distribution);

    double Q_e = cfg->ei.luminosity / st->volume /
                    (electron_energy_average * ELECTRON_ENERGY +
                     proton_energy_average   * PROTON_ENERGY * cfg->ei.eta);

    gsl_spline       *electron_spline      = gsl_spline_alloc(gsl_interp_steffen, electron_aux_pop.size);
    gsl_interp_accel *electron_accelerator = gsl_interp_accel_alloc();
    gsl_spline_init(electron_spline, electron_aux_pop.log_energy, electron_aux_pop.log_population, electron_aux_pop.size);

    gsl_spline       *proton_spline      = gsl_spline_alloc(gsl_interp_steffen, proton_aux_pop.size);
    gsl_interp_accel *proton_accelerator = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, proton_aux_pop.log_energy, proton_aux_pop.log_population, proton_aux_pop.size);

    for(i = 0; i < st->electrons.size; i++)
    {
        double g_electron = st->electrons.energy[i];

        if(g_electron < electron_aux_pop.energy[0] ||
           g_electron > electron_aux_pop.energy[electron_aux_pop.size - 1])
            st->external_injection.electrons[i] = 0.0;
        else
            st->external_injection.electrons[i] = Q_e * exp(gsl_spline_eval(electron_spline, log(g_electron), electron_accelerator));
    }

    for(i = 0; i < st->protons.size; i++)
    {
        double g_proton = st->protons.energy[i];

        if(g_proton < proton_aux_pop.energy[0] ||
           g_proton > proton_aux_pop.energy[proton_aux_pop.size - 1])
            st->external_injection.protons[i] = 0.0;
        else
            st->external_injection.protons[i] =
                cfg->ei.eta * Q_e * exp(gsl_spline_eval(proton_spline, log(g_proton), proton_accelerator));
    }

    st->update_function = &step_experimental_update_populations_injection;

    gsl_spline_free(electron_spline);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(electron_accelerator);
    gsl_interp_accel_free(proton_accelerator);
    free_population(&electron_aux_pop);
    free_population(&proton_aux_pop);
}

void state_init_from_config(state_t *st, config_t *cfg)
{
    unsigned int i;
    unsigned int size = 256;
    unsigned int pion_size = size;
    unsigned int muon_size = size;
    unsigned int neutrino_size = size;

    st->t       = 0;
    st->t_max   = cfg->t_max;
    st->dt      = cfg->dt;
    st->dt_max  = cfg->dt_max;
    st->R       = cfg->R;
    st->density = cfg->density;
    st->eta     = cfg->eta;

    switch(cfg->v)
    {
        case sphere: st->volume = 4. / 3 * M_PI * pow(cfg->R, 3); break;
        case shell:  st->volume = cfg->h * M_PI * pow(cfg->R, 2); break;
        default: break;
    }

    init_state_populations(st,
            cfg->electron_distribution.min, cfg->electron_distribution.max, (unsigned int) cfg->electron_size,
            cfg->proton_distribution.min,   cfg->proton_distribution.max,   (unsigned int) cfg->proton_size,
            cfg->photon_distribution.min,   cfg->photon_distribution.max,   (unsigned int) cfg->photon_size,
            pion_size, muon_size, neutrino_size);
    init_state_synchrotron(st, cfg->magnetic_field);
    init_state_aux_memory(st);

    double t_acc = 1e100;
    init_acceleration(st, &st->electron_acceleration, electron, t_acc);
    init_acceleration(st, &st->positron_acceleration, positron, t_acc);
    init_acceleration(st, &st->proton_acceleration,   proton,   t_acc);

    double t_esc = 1e100;
    switch(cfg->v)
    {
        case sphere: t_esc = 3.   / 4 * cfg->R / LIGHT_SPEED; break;
        case shell:  t_esc = M_PI / 4 * cfg->h / LIGHT_SPEED; break;
        default: break;
    }

    init_state_escape(st, t_esc, cfg->cfe_ratio);
    init_state_decay_and_escape(st, t_esc);

    state_init_LUTs(st);

    state_init_RK_information(st);

    generate_population(&st->electrons, &cfg->electron_distribution);
    generate_population(&st->protons,   &cfg->proton_distribution);

    double particle_density = cfg->density / (cfg->eta * PROTON_MASS + ELECTRON_MASS);

    for(i = 0; i < st->electrons.size; i++)
    {
        st->electrons.population[i] *= particle_density;
        st->electrons.log_population[i] = log(st->electrons.population[i]);
    }

    for(i = 0; i < st->protons.size; i++)
    {
        st->protons.population[i] *= particle_density * cfg->eta;
        st->protons.log_population[i] = log(st->protons.population[i]);
    }

    st->electrons.population[st->electrons.size - 1] = DBL_MIN;
    st->electrons.log_population[st->electrons.size - 1] = log(DBL_MIN);
    st->protons.population[st->protons.size - 1] = DBL_MIN;
    st->protons.log_population[st->protons.size - 1] = log(DBL_MIN);

    st->step_function           = &step;
    st->tentative_step_function = &step_tentative;
    st->update_function = &step_experimental_update_populations;

    if(cfg->ei.luminosity != 0)
        state_init_injection(st, cfg);

#ifdef USE_THREAD_POOL
    thread_pool_init(&st->thread_pool, 5);
#endif
}

void init_state_synchrotron(state_t *st, double B)
{
    st->B = B;

    init_synchrotron(st, &st->electron_synchrotron,            electron);
    init_synchrotron(st, &st->positron_synchrotron,            positron);
    init_synchrotron(st, &st->proton_synchrotron,              proton);
    init_synchrotron(st, &st->positive_pion_synchrotron,       positive_pion);
    init_synchrotron(st, &st->negative_pion_synchrotron,       negative_pion);
    init_synchrotron(st, &st->positive_left_muon_synchrotron,  positive_left_muon);
    init_synchrotron(st, &st->positive_right_muon_synchrotron, positive_right_muon);
    init_synchrotron(st, &st->negative_left_muon_synchrotron,  negative_left_muon);
    init_synchrotron(st, &st->negative_right_muon_synchrotron, negative_right_muon);
}

void init_state_escape(state_t *st, double t, double cfe_ratio)
{
    init_escape(st, &st->photon_escape,   photon,   t);
    init_escape(st, &st->electron_escape, electron, t * cfe_ratio);
    init_escape(st, &st->positron_escape, positron, t * cfe_ratio);
    init_escape(st, &st->proton_escape,   proton,   t * cfe_ratio);

    init_escape(st, &st->electron_neutrino_escape,     electron_neutrino,     t);
    init_escape(st, &st->electron_antineutrino_escape, electron_antineutrino, t);
    init_escape(st, &st->muon_neutrino_escape,         muon_neutrino,         t);
    init_escape(st, &st->muon_antineutrino_escape,     muon_antineutrino,     t);
}

void init_state_decay_and_escape(state_t *st, double t)
{
    init_decay_and_escape(st, &st->neutron_decay_and_escape,  neutron,  t);

    init_decay_and_escape(st, &st->neutral_pion_decay_and_escape,  neutral_pion,  t);
    init_decay_and_escape(st, &st->positive_pion_decay_and_escape, positive_pion, t);
    init_decay_and_escape(st, &st->negative_pion_decay_and_escape, negative_pion, t);

    init_decay_and_escape(st, &st->positive_left_muon_decay_and_escape,  positive_left_muon,  t);
    init_decay_and_escape(st, &st->positive_right_muon_decay_and_escape, positive_right_muon, t);
    init_decay_and_escape(st, &st->negative_left_muon_decay_and_escape,  negative_left_muon,  t);
    init_decay_and_escape(st, &st->negative_right_muon_decay_and_escape, negative_right_muon, t);
}

void init_state_populations(state_t *st,
        double g_lepton_min, double g_lepton_max, unsigned int lepton_size,
        double g_proton_min, double g_proton_max, unsigned int proton_size,
        double e_photon_min, double e_photon_max, unsigned int photon_size,
        unsigned int pion_size, unsigned int muon_size,
        unsigned int neutrino_size)
{
    /*double proton_nu_0   = 3 / (4 * M_PI) * ELECTRON_CHARGE * st->B / (PROTON_MASS * LIGHT_SPEED);*/
    /*double electron_nu_0 = 3 / (4 * M_PI) * ELECTRON_CHARGE * st->B / (ELECTRON_MASS * LIGHT_SPEED);*/

    /* The minimum energy of the photons is set by the lower cutoff in
     * synchrotron production */
    /*
     *double e_photon_min =
     *    0.53 * PLANCK_CONSTANT / ELECTRON_ENERGY *
     *    fmin(proton_nu_0   / g_proton_max,
     *         electron_nu_0 / g_electron_max);
     */

    /* The maximum photon energy is set by the comptonized energy of the
     * highest synchrotron photons */
    /*
     *double e_photon_max =
     *    pow(g_electron_max, 2) *
     *    fmax(proton_nu_0   * pow(g_proton_max, 2),
     *         electron_nu_0 * pow(g_electron_max, 2));
     */

    /*double e_photon_min = pow(10, -12);*/
    /*double e_photon_max = g_electron_max;*/

    /* The limits for pions are derived from the limits for protons
     * using the formulas for pion production.
     * NOTE: the lower limit for pions can be anything bigger than 1
     * because they can lose energy by synchrotron */
    double g_pion_min = fmax(g_proton_min * 0.2 * PROTON_MASS / PION_MASS, pow(10, 0.5));
    double g_pion_max =      g_proton_max * 0.6 * PROTON_MASS / PION_MASS;

    /* The limits for muons are derived from those for pions by
     * the decay formulas.
     * NOTE: The lower limit again can be anything bigger than 1 */
    double g_muon_min = fmax(g_pion_min * PION_MASS / MUON_MASS, pow(10, 0.5));
    double g_muon_max =      g_pion_max * PION_MASS / MUON_MASS;

    /* The limits for neutrinos are derived from those of pions and muons
     * by the decay formulas.
     * NOTE: This is true for the max energy, as we can not obtain a
     * neutrino with an energy higher than that it is available.
     * This is not true for the lower energies though, so I have set
     * as minimum energy that of pions or muons. */
    double g_neutrino_min = fmin(g_pion_min * PION_MASS / ELECTRON_MASS,
                                 g_muon_min * MUON_MASS / ELECTRON_MASS);
                                 /*g_muon_min * MUON_MASS / ELECTRON_MASS * (1 - sqrt(1 - 1/(g_muon_min * g_muon_min))) / 2);*/
    double g_neutrino_max = fmax(g_pion_max * PION_MASS / ELECTRON_MASS * (1 - pow(MUON_MASS / PION_MASS, 2)),
                                 g_muon_max * MUON_MASS / ELECTRON_MASS);

    /*
     *e_photon_min   = pow(10, -17);
     *e_photon_max   = pow(10, 5.3);
     *g_electron_min = pow(10, 3.3);
     *g_electron_max = pow(10, 5.3);
     *g_proton_min   = pow(10, 1);
     *g_proton_max   = pow(10, 5.3);
     *g_pion_min     = pow(10, 0.5);
     *g_pion_max     = pow(10, 10);
     *g_muon_min     = pow(10, 0.5);
     *g_muon_max     = pow(10, 10);
     *g_neutrino_min = pow(10, -10);
     *g_neutrino_max = pow(10, 12);
     */

    init_population(&st->photons,   photon,   e_photon_min, e_photon_max, photon_size);

    init_population(&st->electrons, electron, g_lepton_min, g_lepton_max, lepton_size);
    init_population(&st->positrons, positron, g_lepton_min, g_lepton_max, lepton_size);

    init_population(&st->protons,   proton,   g_proton_min, g_proton_max, proton_size);
    init_population(&st->neutrons,  neutron,  g_proton_min, g_proton_max, proton_size);

    init_population(&st->neutral_pions,  neutral_pion,  g_pion_min, g_pion_max, pion_size);
    init_population(&st->positive_pions, positive_pion, g_pion_min, g_pion_max, pion_size);
    init_population(&st->negative_pions, negative_pion, g_pion_min, g_pion_max, pion_size);

    init_population(&st->positive_left_muons,  positive_left_muon,  g_muon_min, g_muon_max, muon_size);
    init_population(&st->positive_right_muons, positive_right_muon, g_muon_min, g_muon_max, muon_size);
    init_population(&st->negative_left_muons,  negative_left_muon,  g_muon_min, g_muon_max, muon_size);
    init_population(&st->negative_right_muons, negative_right_muon, g_muon_min, g_muon_max, muon_size);

    init_population(&st->electron_neutrinos,     electron_neutrino, g_neutrino_min, g_neutrino_max, neutrino_size);
    init_population(&st->electron_antineutrinos, electron_neutrino, g_neutrino_min, g_neutrino_max, neutrino_size);
    init_population(&st->muon_neutrinos,         muon_neutrino,     g_neutrino_min, g_neutrino_max, neutrino_size);
    init_population(&st->muon_antineutrinos,     muon_neutrino,     g_neutrino_min, g_neutrino_max, neutrino_size);
}

void init_state_aux_memory(state_t *st)
{
    MEM_PREPARE(st->inverse_compton_photon_gains,                st->photons.size);
    MEM_PREPARE(st->inverse_compton_photon_gains_upscattering,   st->photons.size);
    MEM_PREPARE(st->inverse_compton_photon_gains_downscattering, st->photons.size);
    MEM_PREPARE(st->inverse_compton_photon_losses,               st->photons.size);

    MEM_PREPARE(st->inverse_compton_electron_losses,      st->electrons.size);
    MEM_PREPARE(st->inverse_compton_positron_losses,      st->positrons.size);
    MEM_PREPARE(st->inverse_compton_inner_upscattering,   st->electrons.size);
    MEM_PREPARE(st->inverse_compton_inner_downscattering, st->electrons.size);

    MEM_PREPARE(st->pair_production_photon_losses,  st->photons.size);
    MEM_PREPARE(st->pair_production_electron_gains, st->electrons.size);
    MEM_PREPARE(st->pair_production_positron_gains, st->electrons.size);

    /*MEM_PREPARE(st->bethe_heitler_photon_losses, st->electrons.size);*/
    /*MEM_PREPARE(st->bethe_heitler_proton_losses, st->electrons.size);*/
    MEM_PREPARE(st->bethe_heitler_electron_gains, st->electrons.size);
    MEM_PREPARE(st->bethe_heitler_positron_gains, st->electrons.size);

    MEM_PREPARE(st->multi_resonances_neutral_pion_gains , st->neutral_pions.size);
    MEM_PREPARE(st->multi_resonances_positive_pion_gains, st->positive_pions.size);
    MEM_PREPARE(st->multi_resonances_negative_pion_gains, st->negative_pions.size);

    MEM_PREPARE(st->multi_resonances_proton_gains,   st->protons.size);
    MEM_PREPARE(st->multi_resonances_proton_losses,  st->protons.size);
    MEM_PREPARE(st->multi_resonances_neutron_gains,  st->neutrons.size);
    MEM_PREPARE(st->multi_resonances_neutron_losses, st->neutrons.size);

    MEM_PREPARE(st->direct_neutral_pion_gains , st->neutral_pions.size);
    MEM_PREPARE(st->direct_positive_pion_gains, st->positive_pions.size);
    MEM_PREPARE(st->direct_negative_pion_gains, st->negative_pions.size);

    MEM_PREPARE(st->direct_pion_production_proton_gains,   st->protons.size);
    MEM_PREPARE(st->direct_pion_production_proton_losses,  st->protons.size);
    MEM_PREPARE(st->direct_pion_production_neutron_gains,  st->neutrons.size);
    MEM_PREPARE(st->direct_pion_production_neutron_losses, st->neutrons.size);

    MEM_PREPARE(st->pion_decay_positive_left_muon_gains,  st->positive_left_muons.size);
    MEM_PREPARE(st->pion_decay_positive_right_muon_gains, st->positive_right_muons.size);
    MEM_PREPARE(st->pion_decay_negative_left_muon_gains,  st->negative_left_muons.size);
    MEM_PREPARE(st->pion_decay_negative_right_muon_gains, st->negative_right_muons.size);

    MEM_PREPARE(st->pion_decay_photon_gains, st->photons.size);

    MEM_PREPARE(st->pion_decay_muon_neutrino_gains,     st->muon_neutrinos.size);
    MEM_PREPARE(st->pion_decay_muon_antineutrino_gains, st->muon_antineutrinos.size);

    MEM_PREPARE(st->muon_decay_electron_neutrino_gains,     st->electron_neutrinos.size);
    MEM_PREPARE(st->muon_decay_electron_antineutrino_gains, st->electron_antineutrinos.size);
    MEM_PREPARE(st->muon_decay_muon_neutrino_gains,         st->muon_neutrinos.size);
    MEM_PREPARE(st->muon_decay_muon_antineutrino_gains,     st->muon_antineutrinos.size);

    // External Injection
    MEM_PREPARE(st->external_injection.photons,   st->photons.size);
    MEM_PREPARE(st->external_injection.electrons, st->electrons.size);
    MEM_PREPARE(st->external_injection.positrons, st->positrons.size);
    MEM_PREPARE(st->external_injection.protons,   st->protons.size);
    MEM_PREPARE(st->external_injection.neutrons,  st->neutrons.size);

    MEM_PREPARE(st->external_injection.positive_pions, st->positive_pions.size);
    MEM_PREPARE(st->external_injection.neutral_pions,  st->neutral_pions.size);
    MEM_PREPARE(st->external_injection.negative_pions, st->negative_pions.size);

    MEM_PREPARE(st->external_injection.negative_left_muons,  st->negative_left_muons.size);
    MEM_PREPARE(st->external_injection.negative_right_muons, st->negative_right_muons.size);
    MEM_PREPARE(st->external_injection.positive_left_muons,  st->positive_left_muons.size);
    MEM_PREPARE(st->external_injection.positive_right_muons, st->positive_right_muons.size);

    MEM_PREPARE(st->external_injection.electron_neutrinos,     st->electron_neutrinos.size);
    MEM_PREPARE(st->external_injection.electron_antineutrinos, st->electron_antineutrinos.size);
    MEM_PREPARE(st->external_injection.muon_neutrinos,         st->muon_neutrinos.size);
    MEM_PREPARE(st->external_injection.muon_antineutrinos,     st->muon_antineutrinos.size);
}

void state_init_LUTs(state_t *st)
{
    /* We don't need to init or calculate the
     * synchrotron LUTs because that is done when
     * initializing the synchrotron objects */
    init_inverse_compton_LUT_g1_g2(st);
    init_inverse_compton_LUT_losses_reaction_rate(st);
    init_pair_production_LUT_photon_losses(st);
    init_pair_production_LUT_lepton_gains(st);
    init_muon_decay_LUT_neutrino_functions(st);
    init_pion_decay_LUT_muon_functions(st);
    init_multi_resonances_LUT_pion_gains(st);
    init_multi_resonances_LUT_hadron_gains(st);
    init_multi_resonances_LUT_hadron_losses(st);
    init_direct_LUT_pion_gains(st);
    init_direct_pion_production_LUT_hadron_losses(st);
    init_direct_pion_production_LUT_hadron_gains(st);
    init_bethe_heitler_LUT_lepton_gains(st);

    calculate_inverse_compton_LUT_g1_g2(st);
    calculate_inverse_compton_LUT_losses_reaction_rate(st);
    calculate_pair_production_LUT_photon_losses(st);
    calculate_pair_production_LUT_lepton_gains(st);
    calculate_muon_decay_LUT_neutrino_functions(st);
    calculate_pion_decay_LUT_muon_functions(st);
    calculate_multi_resonances_LUT_pion_gains(st);
    calculate_multi_resonances_LUT_hadron_gains(st);
    calculate_multi_resonances_LUT_hadron_losses(st);
    calculate_direct_LUT_pion_gains(st);
    calculate_direct_pion_production_LUT_hadron_losses(st);
    calculate_direct_pion_production_LUT_hadron_gains(st);
    calculate_bethe_heitler_LUT_lepton_gains(st);
}

void state_init_RK_information(state_t *st)
{
#define MEM_PREPARE_EXTRA(X) \
    MEM_PREPARE(st->X##_RK_information.stage0_population, st->X.size); \
    for(unsigned int i = 0; i < 4; i++) \
        MEM_PREPARE(st->X##_RK_information.stage_delta[i], \
                    st->X.size);

    MEM_PREPARE_EXTRA(photons)
    MEM_PREPARE_EXTRA(electrons)
    MEM_PREPARE_EXTRA(positrons)
    MEM_PREPARE_EXTRA(protons)
    MEM_PREPARE_EXTRA(neutrons)

    MEM_PREPARE_EXTRA(neutral_pions)
    MEM_PREPARE_EXTRA(positive_pions)
    MEM_PREPARE_EXTRA(negative_pions)

    MEM_PREPARE_EXTRA(positive_left_muons)
    MEM_PREPARE_EXTRA(positive_right_muons)
    MEM_PREPARE_EXTRA(negative_left_muons)
    MEM_PREPARE_EXTRA(negative_right_muons)

    MEM_PREPARE_EXTRA(electron_neutrinos)
    MEM_PREPARE_EXTRA(electron_antineutrinos)
    MEM_PREPARE_EXTRA(muon_neutrinos)
    MEM_PREPARE_EXTRA(muon_antineutrinos)

#undef MEM_PREPARE_EXTRA
}

void state_report_general_info(state_t *st)
{
    fprintf(stderr,"Density:\t%lg\n", st->density);
    fprintf(stderr,"Proton To Electron Ratio:\t%lg\n", st->eta);
    fprintf(stderr,"B:\t%lg\n", st->B);
    fprintf(stderr,"dt:\t%lg\n", st->dt);
    fprintf(stderr,"Max time:\t%lg\n", st->t_max);
#define FORMAT_ARGUMENTS(X) st->X.energy[0], log10(st->X.energy[0]), st->X.energy[st->X.size - 1], log10(st->X.energy[st->X.size - 1]), st->X.size
    fprintf(stderr, "Photons:  \t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(photons));
    fprintf(stderr, "Protons:  \t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(protons));
    fprintf(stderr, "Neutrons: \t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(neutrons));
    fprintf(stderr, "Electrons:\t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(electrons));
    fprintf(stderr, "Pions:    \t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(neutral_pions));
    fprintf(stderr, "Muons:    \t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(negative_left_muons));
    fprintf(stderr, "Neutrinos:\t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\t(%u)\n", FORMAT_ARGUMENTS(electron_neutrinos));
#undef FORMAT_ARGUMENTS
    fprintf(stderr, "# Proton   Population:\t%lg\n", calculate_population(&st->protons));
    fprintf(stderr, "# Electron Population:\t%lg\n", calculate_population(&st->electrons));

    fprintf(stderr,"\n");
    fprintf(stderr,"Timescales\t     Escape\tAcceleration\tSynchrotron(g_min)\tSynchrotron(g_max)\t   Decay(1)\tDecay(g_min)\tDecay(g_max)\n");
    fprintf(stderr,"Photons:\t%11lg\n", st->photon_escape.t);
    fprintf(stderr,"Electrons:\t%11lg\t%11lg\t\t%11lg\t\t%11lg\n",
            st->electron_escape.t,
            st->electron_acceleration.t,
            1 / (st->electron_synchrotron.particle_losses_factor * st->electrons.energy[0]),
            1 / (st->electron_synchrotron.particle_losses_factor * st->electrons.energy[st->electrons.size - 1]));
    fprintf(stderr,"Protons:\t%11lg\t%11lg\t\t%11lg\t\t%11lg\n",
            st->proton_escape.t,
            st->proton_acceleration.t,
            1 / (st->proton_synchrotron.particle_losses_factor * st->protons.energy[0]),
            1 / (st->proton_synchrotron.particle_losses_factor * st->protons.energy[st->protons.size - 1]));
    fprintf(stderr,"Neutrons:\t%11lg\t\t\t\t\t\t\t\t\t%11lg\t%11lg\t%11lg\n",
            st->neutron_decay_and_escape.escape_lifetime,
            st->neutron_decay_and_escape.decay_lifetime,
            1/st->neutron_decay_and_escape.t[0],
            1/st->neutron_decay_and_escape.t[st->neutrons.size - 1]);
    fprintf(stderr,"Neutral Pions:\t%11lg\t\t\t\t\t\t\t\t\t%11lg\t%11lg\t%11lg\n",
            st->neutral_pion_decay_and_escape.escape_lifetime,
            st->neutral_pion_decay_and_escape.decay_lifetime,
            1/st->neutral_pion_decay_and_escape.t[0],
            1/st->neutral_pion_decay_and_escape.t[st->neutral_pions.size - 1]);
    fprintf(stderr,"Charged Pions:\t%11lg\t\t\t\t%11lg\t\t%11lg\t%11lg\t%11lg\t%11lg\n",
            st->positive_pion_decay_and_escape.escape_lifetime,
            1 / (st->positive_pion_synchrotron.particle_losses_factor * st->positive_pions.energy[0]),
            1 / (st->positive_pion_synchrotron.particle_losses_factor * st->positive_pions.energy[st->positive_pions.size - 1]),
            st->positive_pion_decay_and_escape.decay_lifetime,
            1/st->positive_pion_decay_and_escape.t[0],
            1/st->positive_pion_decay_and_escape.t[st->positive_pions.size - 1]);
    fprintf(stderr,"Muons:  \t%11lg\t\t\t\t%11lg\t\t%11lg\t%11lg\t%11lg\t%11lg\n",
            st->positive_left_muon_decay_and_escape.escape_lifetime,
            1 / (st->positive_left_muon_synchrotron.particle_losses_factor * st->positive_left_muons.energy[0]),
            1 / (st->positive_left_muon_synchrotron.particle_losses_factor * st->positive_left_muons.energy[st->positive_left_muons.size - 1]),
            st->positive_left_muon_decay_and_escape.decay_lifetime,
            1/st->positive_left_muon_decay_and_escape.t[0],
            1/st->positive_left_muon_decay_and_escape.t[st->positive_left_muons.size - 1]);
    fprintf(stderr,"\n\n");
}

// If, at some point, we add the ei data to the state, remove the
// config parameter
void state_report_injection_info(state_t *st, config_t *cfg)
{
    double electron_energy_average = distribution_average(&cfg->ei.electron_distribution);
    double proton_energy_average   = distribution_average(&cfg->ei.proton_distribution);

    double Q_e = cfg->ei.luminosity / st->volume /
                    (electron_energy_average * ELECTRON_ENERGY +
                     proton_energy_average   * PROTON_ENERGY * cfg->ei.eta);

    fprintf(stderr,"Luminosity:\t%lg\n", cfg->ei.luminosity);
    fprintf(stderr,"Injection density:\t%lg\n", Q_e);
    fprintf(stderr,"Proton To Electron Ratio:\t%lg\n", cfg->ei.eta);
    fprintf(stderr,"Photon Luminosity:\t%lg\n", cfg->ei.photon_luminosity);
#define FORMAT_ARGUMENTS(X) \
    "\t%11lg (10^%+4.2lf)\t%11lg (10^%+4.2lf)\n",  \
    cfg->ei.X##_distribution.min, log10(cfg->ei.X##_distribution.min), \
    cfg->ei.X##_distribution.max, log10(cfg->ei.X##_distribution.max)
    fprintf(stderr, "Photons:  " FORMAT_ARGUMENTS(photon));
    fprintf(stderr, "Protons:  " FORMAT_ARGUMENTS(proton));
    fprintf(stderr, "Electrons:" FORMAT_ARGUMENTS(electron));
#undef FORMAT_ARGUMENTS
    fprintf(stderr,"\n\n");
}

void state_print_data_to_file(state_t *st, enum particle_type pt, char *filename)
{
    unsigned int i;

    FILE *temp_file = fopen(filename, "w");

    fprintf(temp_file, "# t: %lg\n", st->t);

    switch(pt)
    {
        case photon:
        {
            fprintf(temp_file,
                    "#Energy\tEnergy_erg\tPopulation\t" \
                    "Injection\t" \
                    "e_sync_gains\tp_sync_gains\tpp_sync_gains\tnp_sync_gains\t" \
                    "plm_sync_gains\tprm_sync_gains\tnlm_sync_gains\tnrm_sync_gains\t" \
                    "IC_gains_up\tIC_gains_down\tPi0_decay\t" \
                    "e_sync_losses\tp_sync_losses\tpp_sync_losses\tnp_sync_losses\t" \
                    "plm_sync_losses\tprm_sync_losses\tnlm_sync_losses\tnrm_sync_losses\t" \
                    "IC_losses\tPP_losses\tesc_losses\n");
            for(i = 0; i < st->photons.size; i++)
            {
                fprintf(temp_file,"%lg\t%lg\t%lg\t",
                        st->photons.energy[i],
                        st->photons.energy[i] * ELECTRON_ENERGY,
                        st->photons.population[i]);

                fprintf(temp_file,
                        "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t",
                        st->external_injection.photons[i],
                        st->electron_synchrotron.photon_gains[i],
                        st->proton_synchrotron.photon_gains[i],
                        st->positive_pion_synchrotron.photon_gains[i],
                        st->negative_pion_synchrotron.photon_gains[i],
                        st->positive_left_muon_synchrotron.photon_gains[i],
                        st->positive_right_muon_synchrotron.photon_gains[i],
                        st->negative_left_muon_synchrotron.photon_gains[i],
                        st->negative_right_muon_synchrotron.photon_gains[i],
                        st->inverse_compton_photon_gains_upscattering[i],
                        st->inverse_compton_photon_gains_downscattering[i],
                        st->pion_decay_photon_gains[i]);

                fprintf(temp_file,
                        "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        st->electron_synchrotron.photon_losses[i],
                        st->proton_synchrotron.photon_losses[i],
                        st->positive_pion_synchrotron.photon_losses[i],
                        st->negative_pion_synchrotron.photon_losses[i],
                        st->positive_left_muon_synchrotron.photon_losses[i],
                        st->positive_right_muon_synchrotron.photon_losses[i],
                        st->negative_left_muon_synchrotron.photon_losses[i],
                        st->negative_right_muon_synchrotron.photon_losses[i],
                        st->inverse_compton_photon_losses[i],
                        st->pair_production_photon_losses[i],
                        st->photon_escape.losses[i]);
            }
            break;
        }

        case proton:
        case neutron:
        {
            fprintf(temp_file, "# Proton Population:\t%lg\n",  calculate_population(&st->protons));
            fprintf(temp_file, "# Neutron Population:\t%lg\n", calculate_population(&st->neutrons));

            fprintf(temp_file,
                    "#Energy\tP_Population\tN_Population\t" \
                    "MR_P_gains\tD_P_gains\tMR_N_gains\tD_N_gains\t" \
                    "MR_P_losses\tD_P_losses\tMR_N_losses\tD_N_losses\t" \
                    "P_escape\tN_escape\tN_decay\n");
            for(i = 0; i < st->protons.size; i++)
            {
                decay_and_escape_t ndae = st->neutron_decay_and_escape;

                fprintf(temp_file,"%lg\t%lg\t%lg\t",
                        st->protons.energy[i],
                        st->protons.population[i],
                        st->neutrons.population[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st->multi_resonances_proton_gains[i],
                        st->direct_pion_production_proton_gains[i],
                        st->multi_resonances_neutron_gains[i],
                        st->direct_pion_production_neutron_gains[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        st->multi_resonances_proton_losses[i],
                        st->direct_pion_production_proton_losses[i],
                        st->multi_resonances_neutron_losses[i],
                        st->direct_pion_production_neutron_losses[i],
                        st->proton_escape.losses[i],
                        ndae.losses[i] / (ndae.t[i] * ndae.escape_lifetime),
                        ndae.losses[i] / (ndae.t[i] * ndae.decay_lifetime * st->neutrons.energy[i]));
            }
            break;
        }

        case electron:
        case positron:
        {
            fprintf(temp_file, "# Electron Population:\t%lg\n", calculate_population(&st->electrons));

            fprintf(temp_file,
                    "#Energy\tPopulation\t" \
                    "External_Injection\t" \
                    "Acc_gains\t" \
                    "Pair_production\t" \
                    "Bethe_Heitler\t" \
                    "Sync_losses\tIC_Losses\tEscape\n");
            for(i = 0; i < st->electrons.size; i++)
            {
                fprintf(temp_file,"%lg\t%lg\t",
                        st->electrons.energy[i], st->electrons.population[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st->external_injection.electrons[i],
                        st->electron_acceleration.gains[i],
                        st->pair_production_electron_gains[i],
                        st->bethe_heitler_electron_gains[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\n",
                        st->electron_synchrotron.particle_losses[i],
                        st->inverse_compton_electron_losses[i],
                        st->electron_escape.losses[i]);
            }
            break;
        }

        case neutral_pion:
        case positive_pion:
        case negative_pion:
        {
            // Remember that we assume that the mass of all the pions is the same
            fprintf(temp_file,
                    "#Energy\t0Pion_Population\t+Pion_Population\t-Pion_Population\t" \
                    "0Pion_MR_gains\t0Pion_D_gains\t" \
                    "+Pion_MR_gains\t+Pion_D_gains\t" \
                    "-Pion_MR_gains\t-Pion_D_gains\t" \
                    "+Pion_sync_losses\t+Pion_decay\t+Pion_escape\t" \
                    "-Pion_sync_losses\t-Pion_decay\t-Pion_escape\t" \
                    "0Pion_decay\t0Pion_escape\n");
            for(i = 0; i < st->neutral_pions.size; i++)
            {
                decay_and_escape_t pos_dae = st->positive_pion_decay_and_escape;
                decay_and_escape_t neg_dae = st->negative_pion_decay_and_escape;
                decay_and_escape_t neu_dae = st->neutral_pion_decay_and_escape;

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st->neutral_pions.energy[i],
                        st->neutral_pions.population[i],
                        st->positive_pions.population[i],
                        st->negative_pions.population[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t",
                        st->multi_resonances_neutral_pion_gains[i],
                        st->direct_neutral_pion_gains[i],
                        st->multi_resonances_positive_pion_gains[i],
                        st->direct_positive_pion_gains[i],
                        st->multi_resonances_negative_pion_gains[i],
                        st->direct_negative_pion_gains[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        st->positive_pion_synchrotron.particle_losses[i],
                        pos_dae.losses[i] / (pos_dae.t[i] * pos_dae.decay_lifetime * st->positive_pions.energy[i]),
                        pos_dae.losses[i] / (pos_dae.t[i] * pos_dae.escape_lifetime),
                        st->negative_pion_synchrotron.particle_losses[i],
                        neg_dae.losses[i] / (neg_dae.t[i] * neg_dae.decay_lifetime * st->negative_pions.energy[i]),
                        neg_dae.losses[i] / (neg_dae.t[i] * neg_dae.escape_lifetime),
                        neu_dae.losses[i] / (neu_dae.t[i] * neu_dae.decay_lifetime * st->neutral_pions.energy[i]),
                        neu_dae.losses[i] / (neu_dae.t[i] * neu_dae.escape_lifetime));

            }
            break;
        }

        case positive_left_muon:
        case positive_right_muon:
        case negative_left_muon:
        case negative_right_muon:
        {
            // Remember that the mass of all the muons is the same
            fprintf(temp_file,
                    "#Energy\t"
                    "+lMuon_Population\t+rMuon_Population\t" \
                    "-lMuon_Population\t-rMuon_Population\t" \
                    "+lPD_gains\t+rPD_gains\t" \
                    "-lPD_gains\t-rPD_gains\t" \
                    "+l_sync_losses\t+r_sync_losses\t" \
                    "-l_sync_losses\t-r_sync_losses\t" \
                    "+l_decay\t+r_decay\t" \
                    "-l_decay\t-r_decay\t" \
                    "+l_escape\t+r_escape\t" \
                    "-l_escape\t-r_escape\n");
            for(i = 0; i < st->positive_left_muons.size; i++)
            {
                decay_and_escape_t pos_l_dae = st->positive_left_muon_decay_and_escape;
                decay_and_escape_t pos_r_dae = st->positive_right_muon_decay_and_escape;
                decay_and_escape_t neg_l_dae = st->negative_left_muon_decay_and_escape;
                decay_and_escape_t neg_r_dae = st->negative_right_muon_decay_and_escape;

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t",
                        st->positive_left_muons.energy[i],
                        st->positive_left_muons.population[i],
                        st->positive_right_muons.population[i],
                        st->negative_left_muons.population[i],
                        st->negative_right_muons.population[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st->pion_decay_positive_left_muon_gains[i],
                        st->pion_decay_positive_right_muon_gains[i],
                        st->pion_decay_negative_left_muon_gains[i],
                        st->pion_decay_negative_right_muon_gains[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st->positive_left_muon_synchrotron.particle_losses[i],
                        st->positive_right_muon_synchrotron.particle_losses[i],
                        st->negative_left_muon_synchrotron.particle_losses[i],
                        st->negative_right_muon_synchrotron.particle_losses[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        pos_l_dae.losses[i] / (pos_l_dae.t[i] * pos_l_dae.decay_lifetime * st->positive_left_muons.energy[i]),
                        pos_r_dae.losses[i] / (pos_r_dae.t[i] * pos_r_dae.decay_lifetime * st->positive_right_muons.energy[i]),
                        neg_l_dae.losses[i] / (neg_l_dae.t[i] * neg_l_dae.decay_lifetime * st->negative_left_muons.energy[i]),
                        neg_r_dae.losses[i] / (neg_r_dae.t[i] * neg_r_dae.decay_lifetime * st->negative_right_muons.energy[i]),
                        pos_l_dae.losses[i] / (pos_l_dae.t[i] * pos_l_dae.escape_lifetime),
                        pos_r_dae.losses[i] / (pos_r_dae.t[i] * pos_r_dae.escape_lifetime),
                        neg_l_dae.losses[i] / (neg_l_dae.t[i] * neg_l_dae.escape_lifetime),
                        neg_r_dae.losses[i] / (neg_r_dae.t[i] * neg_r_dae.escape_lifetime));
            }
            break;
        }

        case electron_neutrino:
        case electron_antineutrino:
        case muon_neutrino:
        case muon_antineutrino:
        {
            // Remember that we normalize the energy of the neutrinos with the mass of the electron
            fprintf(temp_file,
                    "#Energy\t"
                    "e_neutrino_Population\te_antineutrino_Population\t" \
                    "mu_neutrino_Population\tmu_antineutrino_Population\t" \
                    "e_neutrino_gains\te_antineutrino_gains\t" ///Remember, only from muons
                    "mu_neutrino_PD_gains\tmu_neutrino_MuD_gains\t"
                    "mu_antineutrino_PD_gains\tmu_antineutrino_MuD_gains\t" \
                    "e_neutrino_escape\te_antineutrino_escape\t" \
                    "mu_neutrino_escape\tmu_antineutrino_escape\n");
            for(i = 0; i < st->electron_neutrinos.size; i++)
            {
                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t",
                        st->electron_neutrinos.energy[i],
                        st->electron_neutrinos.population[i],
                        st->electron_antineutrinos.population[i],
                        st->muon_neutrinos.population[i],
                        st->muon_antineutrinos.population[i]);

                fprintf(temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        st->muon_decay_electron_neutrino_gains[i],
                        st->muon_decay_electron_antineutrino_gains[i],
                        st->pion_decay_muon_neutrino_gains[i],
                        st->muon_decay_muon_neutrino_gains[i],
                        st->pion_decay_muon_antineutrino_gains[i],
                        st->muon_decay_muon_antineutrino_gains[i],
                        st->electron_neutrino_escape.losses[i],
                        st->electron_antineutrino_escape.losses[i],
                        st->muon_neutrino_escape.losses[i],
                        st->muon_antineutrino_escape.losses[i]);
            }
            break;
        }

        default:
            assert(0);
    }

    fclose(temp_file);
}

void state_save_state_to_file(state_t *st, char *filename)
{
    FILE *fout = fopen(filename, "w");

    fwrite(&st->B,       sizeof(double), 1, fout);
    fwrite(&st->R,       sizeof(double), 1, fout);
    fwrite(&st->density, sizeof(double), 1, fout);
    fwrite(&st->dt,      sizeof(double), 1, fout);
    fwrite(&st->dt_max,  sizeof(double), 1, fout);

#define WRITE_POPULATION(X) \
    fwrite(&st->X.size,       sizeof(unsigned int), 1,          fout); \
    fwrite( st->X.energy,     sizeof(double),       st->X.size, fout); \
    fwrite( st->X.population, sizeof(double),       st->X.size, fout);

    WRITE_POPULATION(photons)
    WRITE_POPULATION(electrons)
    WRITE_POPULATION(protons)
    WRITE_POPULATION(neutrons)

    WRITE_POPULATION(neutral_pions)
    WRITE_POPULATION(positive_pions)
    WRITE_POPULATION(negative_pions)

    WRITE_POPULATION(positive_left_muons)
    WRITE_POPULATION(positive_right_muons)
    WRITE_POPULATION(negative_left_muons)
    WRITE_POPULATION(negative_right_muons)

    WRITE_POPULATION(electron_neutrinos)
    WRITE_POPULATION(electron_antineutrinos)
    WRITE_POPULATION(muon_neutrinos)
    WRITE_POPULATION(muon_antineutrinos)
#undef WRITE_POPULATION

    fwrite(&st->proton_acceleration.t,   sizeof(double), 1, fout);
    fwrite(&st->electron_acceleration.t, sizeof(double), 1, fout);
    fwrite(&st->electron_escape.t,       sizeof(double), 1, fout);

    fclose(fout);
}

void state_load_state_from_file(state_t *st, char *filename)
{
    unsigned int i;
    FILE *fin = fopen(filename, "r");

    if(fin == NULL)
    {
        fprintf(stderr,"Error opening state file: %s\n", filename);
        assert(0);
    }

    fread(&st->B,       sizeof(double), 1, fin);
    fread(&st->R,       sizeof(double), 1, fin);
    fread(&st->density, sizeof(double), 1, fin);
    fread(&st->dt,      sizeof(double), 1, fin);
    fread(&st->dt_max,  sizeof(double), 1, fin);

#define READ_POPULATION(X) \
    fread(&st->X.size,       sizeof(unsigned int), 1,          fin); \
    posix_memalign((void **) &st->X.energy,               32, sizeof(double) * st->X.size); \
    posix_memalign((void **) &st->X.log_energy,           32, sizeof(double) * st->X.size); \
    posix_memalign((void **) &st->X.population,           32, sizeof(double) * st->X.size); \
    posix_memalign((void **) &st->X.log_population,       32, sizeof(double) * st->X.size); \
    posix_memalign((void **) &st->X.tentative_population, 32, sizeof(double) * st->X.size); \
    fread( st->X.energy,     sizeof(double),       st->X.size, fin); \
    fread( st->X.population, sizeof(double),       st->X.size, fin); \
    for(i = 0; i < st->X.size; i++) st->X.log_energy[i] = log(st->X.energy[i]); \
    for(i = 0; i < st->X.size; i++) st->X.log_population[i] = log(st->X.population[i]);

    READ_POPULATION(photons)
    READ_POPULATION(electrons)
    READ_POPULATION(protons)
    READ_POPULATION(neutrons)

    READ_POPULATION(neutral_pions)
    READ_POPULATION(positive_pions)
    READ_POPULATION(negative_pions)

    READ_POPULATION(positive_left_muons)
    READ_POPULATION(positive_right_muons)
    READ_POPULATION(negative_left_muons)
    READ_POPULATION(negative_right_muons)

    READ_POPULATION(electron_neutrinos)
    READ_POPULATION(electron_antineutrinos)
    READ_POPULATION(muon_neutrinos)
    READ_POPULATION(muon_antineutrinos)
#undef READ_POPULATION

    st->photons.mass   = ELECTRON_MASS;
    st->electrons.mass = ELECTRON_MASS;
    st->protons.mass   = PROTON_MASS;
    st->neutrons.mass  = PROTON_MASS;

    st->neutral_pions.mass  = PION_MASS;
    st->positive_pions.mass = PION_MASS;
    st->negative_pions.mass = PION_MASS;

    st->positive_left_muons.mass  = MUON_MASS;
    st->positive_right_muons.mass = MUON_MASS;
    st->negative_left_muons.mass  = MUON_MASS;
    st->negative_right_muons.mass = MUON_MASS;

    st->electron_neutrinos.mass     = ELECTRON_MASS;
    st->electron_antineutrinos.mass = ELECTRON_MASS;
    st->muon_neutrinos.mass         = ELECTRON_MASS;
    st->muon_antineutrinos.mass     = ELECTRON_MASS;

    init_state_synchrotron(st, st->B);
    init_state_aux_memory(st);
    state_init_RK_information(st);

    double proton_acceleration_timescale;
    double electron_acceleration_timescale;
    double electron_escape_timescale;
    fread(&proton_acceleration_timescale,   sizeof(double), 1, fin);
    fread(&electron_acceleration_timescale, sizeof(double), 1, fin);
    fread(&electron_escape_timescale,       sizeof(double), 1, fin);

    init_acceleration(st, &st->proton_acceleration,   proton,   proton_acceleration_timescale);
    init_acceleration(st, &st->electron_acceleration, electron, electron_acceleration_timescale);

    init_state_escape(st, electron_escape_timescale, 1);

    state_init_LUTs(st);

    fclose(fin);
}
