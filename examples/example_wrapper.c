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

    return out;
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    unsigned int i;

    state_t st;
    config_t cfg;

    if(argc != 2)
    {
        fprintf(stderr,"Please, pass as first argument a valid config file\n");
        fprintf(stderr,"Usage: %s <config_file>\n", argv[0]);

        return 1;
    }

    config_read_file(&cfg, "default.toml");
    config_read_file(&cfg, argv[1]);

    state_init_from_config(&st, &cfg);

    double tol = st.tol;
    unsigned int i_max = 1e7;

    bool out = false;

    i = 0;
    while(true)
    {
        st.tentative_step_function(&st, i % 2 == 0);

        if(i % 20 == 0)
            out = state_check_steady_state(&st, i, tol);

        if(st.t > st.t_max)
            out = true;

        if(out || i == i_max) break;

        i++;
    }

    for(i = 0; i < st.photons.size; i++)
    {
        printf("%lg\t%lg\n", st.photons.energy[i], st.photons.population[i]);
    }

#ifdef USE_THREAD_POOL
    thread_pool_clear(&st.thread_pool);
#endif

    return 0;
}
