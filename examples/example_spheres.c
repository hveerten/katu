#include "state.h"
#include "state_step.h"
#include "synchrotron.h"
#include "inverse_compton.h"
#include "muon_decay.h"
#include "pair_production.h"
#include "pion_decay.h"
#include "pion_production.h"
#include "constants.h"

#include "config.h"

#include "state_step_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <float.h>

#include <assert.h>

#include <unistd.h>

static bool state_check_steady_state(state_t *st, unsigned int step, double tol)
{
    unsigned int i;
    bool out = false;

    double photon_statistic   = 0;
    double electron_statistic = 0;
    /*double neutrino_statistic = 0;*/
    double neutrino_statistic1 = 0;
    double neutrino_statistic2 = 0;
    double neutrino_statistic3 = 0;
    double neutrino_statistic4 = 0;

    for(i = 0; i < st->photons.size; i++)
    {
        double n  = st->photons.population[i];
        double tn = st->photons.tentative_population[i];

        double Q = st->external_injection.photons[i] +
                   st->electron_synchrotron.photon_gains[i] +
                   st->positron_synchrotron.photon_gains[i] +
                   st->proton_synchrotron.photon_gains[i] +
                   st->positive_pion_synchrotron.photon_gains[i] +
                   st->negative_pion_synchrotron.photon_gains[i] +
                   st->positive_left_muon_synchrotron.photon_gains[i] +
                   st->positive_right_muon_synchrotron.photon_gains[i] +
                   st->negative_left_muon_synchrotron.photon_gains[i] +
                   st->negative_right_muon_synchrotron.photon_gains[i] +
                   st->inverse_compton_photon_gains[i] +
                   st->pair_annihilation_photon_gains[i];

        Q += st->pion_decay_photon_gains[i];

        double L = (st->electron_synchrotron.photon_losses[i] +
                    st->positron_synchrotron.photon_losses[i] +
                    st->proton_synchrotron.photon_losses[i] +
                    st->positive_pion_synchrotron.photon_losses[i] +
                    st->negative_pion_synchrotron.photon_losses[i] +
                    st->positive_left_muon_synchrotron.photon_losses[i] +
                    st->positive_right_muon_synchrotron.photon_losses[i] +
                    st->negative_left_muon_synchrotron.photon_losses[i] +
                    st->negative_right_muon_synchrotron.photon_losses[i] +
                    st->inverse_compton_photon_losses[i] +
                    st->pair_production_photon_losses[i]) / tn +
                   - 1 / st->photon_escape.t;

        double prd3 =
            (st->external_injection.photons[i] + st->pion_decay_photon_gains[i]
             - tn / st->photon_escape.t) +

             (st->electron_synchrotron.photon_gains[i] + st->electron_synchrotron.photon_losses[i]) +
             (st->positron_synchrotron.photon_gains[i] + st->positron_synchrotron.photon_losses[i]) +
             (st->proton_synchrotron.photon_gains[i]   + st->proton_synchrotron.photon_losses[i]) +
             (st->positive_pion_synchrotron.photon_gains[i] + st->positive_pion_synchrotron.photon_losses[i]) +
             (st->negative_pion_synchrotron.photon_gains[i] + st->negative_pion_synchrotron.photon_losses[i]) +
             (st->positive_left_muon_synchrotron.photon_gains[i]  + st->positive_left_muon_synchrotron.photon_losses[i]) +
             (st->positive_right_muon_synchrotron.photon_gains[i] + st->positive_right_muon_synchrotron.photon_losses[i]) +
             (st->negative_left_muon_synchrotron.photon_gains[i]  + st->negative_left_muon_synchrotron.photon_losses[i]) +
             (st->negative_right_muon_synchrotron.photon_gains[i] + st->negative_right_muon_synchrotron.photon_losses[i]) +

             (st->inverse_compton_photon_gains[i] + st->inverse_compton_photon_losses[i]) +
             (st->pair_annihilation_photon_gains[i] + st->pair_production_photon_losses[i]);

        prd3 /= tn;


        double aux0 = expm1(L * st->dt);

        double prd1 = (n - tn) / tn;
        double prd2 = (aux0 * (tn + Q/L)) / tn;

        if(tn < 1e-300) continue;
        double dn_ndt = L + Q / tn;

        photon_statistic += pow(dn_ndt, 2);
    }

    for(i = 0; i < st->electrons.size; i++)
    {
        double n  = st->electrons.population[i];
        double tn = st->electrons.tentative_population[i];

        double Dn_n = (tn - n) / n;

        electron_statistic += pow(Dn_n, 2);
    }
    electron_statistic /= (st->dt * st->dt);

    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        double n_ma  = st->muon_antineutrinos.population[i];
        double ln_ma = st->muon_antineutrinos.log_population[i];

        double tn_ma  = st->muon_antineutrinos.tentative_population[i];
        double ltn_ma = log(st->muon_antineutrinos.tentative_population[i]);

        double Q_ma = st->external_injection.muon_antineutrinos[i] +
                      st->pion_decay_muon_antineutrino_gains[i] +
                      st->muon_decay_muon_antineutrino_gains[i];

        double L = -1 / st->electron_neutrino_escape.t;

        double aux0 = expm1(L * st->dt);
        double aux1 = aux0 / L;

        double exp_delta = (aux0 * (n_ma + Q_ma / L) ) / (st->dt * n_ma);

        double Dn_ndt = ( tn_ma -  n_ma) / (n_ma * st->dt);
        double Dln_dt = (ltn_ma - ln_ma) / st->dt;

        double dn_ndt = L + Q_ma / n_ma;

        if(n_ma < 1e-300) continue;

        neutrino_statistic1 += pow(exp_delta, 2);
        neutrino_statistic2 += pow(Dn_ndt, 2);
        neutrino_statistic3 += pow(Dln_dt, 2);
        neutrino_statistic4 += pow(dn_ndt, 2);
    }

    if(sqrt(photon_statistic)   < tol &&
       sqrt(electron_statistic) < tol &&
       sqrt(neutrino_statistic4) < tol)
        out = true;

    fprintf(stderr,"%04u\t(%lg)\tDONE\tdt %lg\tphoton stat %lg\telectron stat %lg\tneutrino stat %lg\n",
            step / 20, st->t, st->dt, sqrt(photon_statistic), sqrt(electron_statistic), sqrt(neutrino_statistic4));

    return out;
}

void step_tentative_2(state_t *st1, state_t *st2, bool try_new_step)
{
    double dt_old = st1->dt;
    double dt_new = fmin(st1->dt * 1.1, st1->dt_max);

    if(try_new_step)
    {
        step_propagate_new_dt(st1, dt_new);
        step_propagate_new_dt(st2, dt_new);
    }
    else
        dt_old = st1->dt / 1.1;

    step_calculate_processes(st1);
    st1->update_function(st1, st1->dt);

    if(step_check_tentative_populations(st1))
        goto abort;

    step_calculate_processes(st2);
    st2->update_function(st2, st2->dt);

    if(step_check_tentative_populations(st2))
        goto abort;

    step_fix_tentative_zeros(st1);
    step_fix_tentative_zeros(st2);

    goto end;

abort:
    /* We have aborted the new timestep.
     * So, set the old one and recalculate LUTs.
     * Finally, do a step with the old timestep */
    step_propagate_new_dt(st1, dt_old);
    step_propagate_new_dt(st2, dt_old);
    /*fprintf(stderr,"Back To: %lg\n", dt_old);*/

    step_tentative_2(st1, st2, 0);
    return;

end:;
    step_accept_tentative_population(st1);
    step_accept_tentative_population(st2);
    step_log_populations(st1);
    step_log_populations(st2);

    st1->t += st1->dt;
    st2->t += st2->dt;
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    unsigned int i, j;

    state_t st_inner;
    state_t st_outer;
    config_t cfg_inner;
    config_t cfg_outer;

    if(argc != 3)
    {
        fprintf(stderr,"Please, pass as first and second arguments valid config files\n");
        fprintf(stderr,"Usage: %s <config_file_inner> <config_file_outer\n", argv[0]);

        return 1;
    }

    config_read_file(&cfg_inner, "default.toml");
    config_read_file(&cfg_inner, argv[1]);
    config_read_file(&cfg_outer, "default.toml");
    config_read_file(&cfg_outer, argv[2]);

    // Make sure that both volumes are spheres and that the inner one is
    // smaller than the outer one
    assert(cfg_inner.v == sphere);
    assert(cfg_outer.v == sphere);
    assert(cfg_inner.R <= cfg_outer.R);

    // Make sure that all the bin sizes are equal
    cfg_inner.photon_size = fmax(cfg_inner.photon_size, cfg_outer.photon_size);
    cfg_outer.photon_size = cfg_inner.photon_size;
    cfg_inner.electron_size = fmax(cfg_inner.electron_size, cfg_outer.electron_size);
    cfg_outer.electron_size = cfg_inner.electron_size;
    cfg_inner.proton_size = fmax(cfg_inner.proton_size, cfg_outer.proton_size);
    cfg_outer.proton_size = cfg_inner.proton_size;

    // Make sure that all the times are equal
    cfg_inner.dt = fmin(cfg_inner.dt, cfg_outer.dt);
    cfg_outer.dt = cfg_inner.dt;
    cfg_inner.t_max = fmin(cfg_inner.t_max, cfg_outer.t_max);
    cfg_outer.t_max = cfg_inner.t_max;

    state_init_from_config(&st_inner, &cfg_inner);
    state_init_from_config(&st_outer, &cfg_outer);

    // dt_max is clamped by the escape timescale so we need to equalize it
    // AFTER initing the states
    cfg_inner.dt_max = fmin(cfg_inner.dt_max, cfg_outer.dt_max);
    cfg_outer.dt_max = cfg_inner.dt_max;

    fprintf(stderr, "Inner Sphere Info\n");
    state_report_general_info(&st_inner);
    if(cfg_inner.ei.luminosity != 0)
        state_report_injection_info(&st_inner, &cfg_inner);

    fprintf(stderr, "Outer Sphere Info\n");
    state_report_general_info(&st_outer);
    if(cfg_outer.ei.luminosity != 0)
        state_report_injection_info(&st_outer, &cfg_outer);

    unsigned int i_max = 1e7;

    bool out = false;

    i = 0;
    while(true)
    {
        double surface_ratio = pow(st_inner.R / st_outer.R, 2);

        st_outer.tentative_step_function(&st_outer, i % 2 == 0);
        st_inner.tentative_step_function(&st_inner, i % 2 == 0);

        // Put particles from inside out and vice-versa
#define INTER_VOLUME_INJECTION(X)                               \
        for(j = 0; j < st_inner.X##s.size; j++)                 \
        {                                                       \
            st_inner.inter_volume_injection.X##s[j] =           \
                st_outer.X##s.population[j] * surface_ratio /   \
                st_outer.X##_escape.t * st_inner.dt;            \
            st_outer.inter_volume_injection.X##s[j] =           \
                st_inner.X##s.population[j]  /                  \
                st_inner.X##_escape.t * st_inner.dt;            \
        }
        INTER_VOLUME_INJECTION(photon);
        INTER_VOLUME_INJECTION(electron);
        INTER_VOLUME_INJECTION(proton);
#undef INTER_VOLUME_INJECTION

        if(i % 20 == 0)
        {
            out  = state_check_steady_state(&st_outer, i, st_outer.tol);
            out &= state_check_steady_state(&st_inner, i, st_inner.tol);
        }

        if(st_inner.t > st_inner.t_max)
            out = true;

        if(out || i == i_max) break;

        i++;
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"%04u\t(%lg)\tDONE\tdt %lg\n", i / 100, st_outer.t, st_outer.dt);

    for(i = 0; i < st_outer.photons.size; i++)
    {
        printf("%lg\t%lg\n", st_outer.photons.energy[i], st_outer.photons.population[i]);
    }

#ifdef USE_THREAD_POOL
    thread_pool_clear(&st_inner.thread_pool);
    thread_pool_clear(&st_outer.thread_pool);
#endif

    return 0;
}
