#include "state.h"
#include "fast_start.h"
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

#include <sys/stat.h>

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

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    unsigned int i;

    state_t st;
    config_t cfg;

    config_read_file(&cfg, "default.toml");
    config_read_file(&cfg, argv[1]);

    if(1)
        state_init_from_config(&st, &cfg);
    else
        state_load_state_from_file(&st, "./data/state.st");

    state_report_general_info(&st);
    if(cfg.ei.luminosity != 0)
        state_report_injection_info(&st, &cfg);

    double tol = st.tol;
    unsigned int i_max = 1e7;

    bool out = false;

    mkdir("./data/", 0744);
    st.t = 0; i = 0;
    while(true)
    {
        st.tentative_step_function(&st, i % 2 == 0);

        if(i % 20 == 0)
            out = state_check_steady_state(&st, i, tol);

        if((i % 20) == 0 || out == true)
        {
            char photon_temp_file_filename[128];
            char lepton_temp_file_filename[128];
            char positron_temp_file_filename[128];
            char hadron_temp_file_filename[128];
            char pion_temp_file_filename[128];
            char muon_temp_file_filename[128];
            char neutrino_temp_file_filename[128];

            sprintf(photon_temp_file_filename,   "data/photon_data_%04u.tsv",   i/20);
            sprintf(lepton_temp_file_filename,   "data/lepton_data_%04u.tsv",   i/20);
            sprintf(positron_temp_file_filename, "data/positron_data_%04u.tsv", i/20);
            sprintf(hadron_temp_file_filename,   "data/hadron_data_%04u.tsv",   i/20);
            sprintf(pion_temp_file_filename,     "data/pion_data_%04u.tsv",     i/20);
            sprintf(muon_temp_file_filename,     "data/muon_data_%04u.tsv",     i/20);
            sprintf(neutrino_temp_file_filename, "data/neutrino_data_%04u.tsv", i/20);

            state_print_data_to_file(&st, photon,             photon_temp_file_filename);
            state_print_data_to_file(&st, electron,           lepton_temp_file_filename);
            state_print_data_to_file(&st, positron,           positron_temp_file_filename);
            state_print_data_to_file(&st, proton,             hadron_temp_file_filename);
            state_print_data_to_file(&st, positive_pion,      pion_temp_file_filename);
            state_print_data_to_file(&st, positive_left_muon, muon_temp_file_filename);
            state_print_data_to_file(&st, electron_neutrino,  neutrino_temp_file_filename);
        }
        if(out || i == i_max) break;

        i++;
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"%04u\t(%lg)\tDONE\tdt %lg\n", i / 100, st.t, st.dt);

    for(i = 0; i < st.photons.size; i++)
    {
        printf("%lg\t%lg\n", st.photons.energy[i], st.photons.population[i]);
    }

#ifdef USE_THREAD_POOL
    thread_pool_clear(&st.thread_pool);
#endif

    return 0;
}
