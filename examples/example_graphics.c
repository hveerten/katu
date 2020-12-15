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

    if(argc != 2)
    {
        fprintf(stderr,"Please, pass as first argument a valid config file\n");
        fprintf(stderr,"Usage: %s <config_file>\n", argv[0]);

        return 1;
    }

    config_read_file(&cfg, "default.toml");
    config_read_file(&cfg, argv[1]);

    state_init_from_config(&st, &cfg);

    state_report_general_info(&st);
    if(cfg.ei.luminosity != 0)
        state_report_injection_info(&st, &cfg);

#define GRAPHICAL_OUTPUT 1
#if GRAPHICAL_OUTPUT
    FILE *gnuplot_pipe;
    FILE *photon_temp_file;
    FILE *lepton_temp_file;
    FILE *positron_temp_file;
    FILE *hadron_temp_file;
    FILE *pion_temp_file;
    FILE *muon_temp_file;
    FILE *neutrino_temp_file;
    char photon_temp_file_filename[L_tmpnam]   = "/tmp/MK2_XXXXXX";
    char lepton_temp_file_filename[L_tmpnam]   = "/tmp/MK2_XXXXXX";
    char positron_temp_file_filename[L_tmpnam]   = "/tmp/MK2_XXXXXX";
    char hadron_temp_file_filename[L_tmpnam]   = "/tmp/MK2_XXXXXX";
    char pion_temp_file_filename[L_tmpnam]     = "/tmp/MK2_XXXXXX";
    char muon_temp_file_filename[L_tmpnam]     = "/tmp/MK2_XXXXXX";
    char neutrino_temp_file_filename[L_tmpnam] = "/tmp/MK2_XXXXXX";
    photon_temp_file   = fdopen(mkstemp(photon_temp_file_filename), "w");
    lepton_temp_file   = fdopen(mkstemp(lepton_temp_file_filename), "w");
    positron_temp_file = fdopen(mkstemp(positron_temp_file_filename), "w");
    hadron_temp_file   = fdopen(mkstemp(hadron_temp_file_filename), "w");
    pion_temp_file     = fdopen(mkstemp(pion_temp_file_filename), "w");
    muon_temp_file     = fdopen(mkstemp(muon_temp_file_filename), "w");
    neutrino_temp_file = fdopen(mkstemp(neutrino_temp_file_filename), "w");
    gnuplot_pipe = popen("gnuplot -p", "w");
    fprintf(gnuplot_pipe,"unset object 101\n");
    fprintf(gnuplot_pipe,"set termoption enhanced\n");
    fprintf(gnuplot_pipe,"set xrange [%lg:%lg]\n", st.photons.energy[0] * ELECTRON_ENERGY_EV, st.protons.energy[st.protons.size - 1] * PROTON_ENERGY_EV);
    fprintf(gnuplot_pipe,"set x2range [%lg:%lg]\n", st.photons.energy[0] * ELECTRON_ENERGY / PLANCK_CONSTANT, st.protons.energy[st.protons.size - 1] * PROTON_ENERGY / PLANCK_CONSTANT);
    fprintf(gnuplot_pipe,"set yrange [%lg:%lg]\n", 1e-45, 1e3);

    fprintf(gnuplot_pipe,"set grid\n");
    fprintf(gnuplot_pipe,"set logscale xy\n");
    fprintf(gnuplot_pipe,"set logscale x2\n");
    fprintf(gnuplot_pipe,"set x2tics mirror\n");
#endif

    double tol = st.tol;
    unsigned int i_max = 1e7;

    bool out = false;

    mkdir("./data/", 0744);
    i = 0;
    while(true)
    {
        st.tentative_step_function(&st, i % 2 == 0);

        if(i % 20 == 0)
            out = state_check_steady_state(&st, i, tol);

#if GRAPHICAL_OUTPUT
        if(i % 2 == 0)
        {
            unsigned int j;

            freopen(photon_temp_file_filename, "w", photon_temp_file);
            freopen(lepton_temp_file_filename, "w", lepton_temp_file);
            freopen(positron_temp_file_filename, "w", positron_temp_file);
            freopen(hadron_temp_file_filename, "w", hadron_temp_file);
            freopen(pion_temp_file_filename, "w", pion_temp_file);
            freopen(muon_temp_file_filename, "w", muon_temp_file);
            freopen(neutrino_temp_file_filename, "w", neutrino_temp_file);

            for(j = 0; j < st.photons.size; j++)
            {
                double gains_temp =
                    st.electron_synchrotron.photon_gains[j] +
                    st.proton_synchrotron.photon_gains[j] +
                    st.positive_pion_synchrotron.photon_gains[j] +
                    st.negative_pion_synchrotron.photon_gains[j] +
                    st.positive_left_muon_synchrotron.photon_gains[j] +
                    st.positive_right_muon_synchrotron.photon_gains[j] +
                    st.negative_left_muon_synchrotron.photon_gains[j] +
                    st.negative_right_muon_synchrotron.photon_gains[j] +
                    st.inverse_compton_photon_gains[j];
                double gains_aux = st.photons.population[j] / gains_temp;

                double losses_temp =
                    st.electron_synchrotron.photon_losses[j] +
                    st.proton_synchrotron.photon_losses[j] +
                    st.positive_pion_synchrotron.photon_losses[j] +
                    st.negative_pion_synchrotron.photon_losses[j] +
                    st.positive_left_muon_synchrotron.photon_losses[j] +
                    st.positive_right_muon_synchrotron.photon_losses[j] +
                    st.negative_left_muon_synchrotron.photon_losses[j] +
                    st.negative_right_muon_synchrotron.photon_losses[j] +
                    st.inverse_compton_photon_losses[j] +
                    st.photon_escape.losses[j];
                double losses_aux = st.photons.population[j] / losses_temp;

                fprintf(photon_temp_file,"%lg\t%lg\t%lg\t%lg\t",
                        st.photons.energy[j] * ELECTRON_ENERGY_EV,
                        st.photons.energy[j] * ELECTRON_ENERGY / PLANCK_CONSTANT,
                        st.photons.energy[j] * sqrt(ELECTRON_ENERGY),
                        st.photons.population[j]);

                fprintf(photon_temp_file,
                        "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t",
                        gains_aux * st.electron_synchrotron.photon_gains[j],
                        gains_aux * st.proton_synchrotron.photon_gains[j],
                        gains_aux * st.positive_pion_synchrotron.photon_gains[j],
                        gains_aux * st.negative_pion_synchrotron.photon_gains[j],
                        gains_aux * st.positive_left_muon_synchrotron.photon_gains[j],
                        gains_aux * st.positive_right_muon_synchrotron.photon_gains[j],
                        gains_aux * st.negative_left_muon_synchrotron.photon_gains[j],
                        gains_aux * st.negative_right_muon_synchrotron.photon_gains[j],
                        gains_aux * st.inverse_compton_photon_gains[j]);

                fprintf(photon_temp_file,
                        "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        losses_aux * st.electron_synchrotron.photon_losses[j],
                        losses_aux * st.proton_synchrotron.photon_losses[j],
                        losses_aux * st.positive_pion_synchrotron.photon_losses[j],
                        losses_aux * st.negative_pion_synchrotron.photon_losses[j],
                        losses_aux * st.positive_left_muon_synchrotron.photon_losses[j],
                        losses_aux * st.positive_right_muon_synchrotron.photon_losses[j],
                        losses_aux * st.negative_left_muon_synchrotron.photon_losses[j],
                        losses_aux * st.negative_right_muon_synchrotron.photon_losses[j],
                        losses_aux * st.inverse_compton_photon_losses[j],
                        losses_aux * st.photon_escape.losses[j]);
            }

            for(j = 0; j < st.electrons.size; j++)
                fprintf(lepton_temp_file,"%lg\t%lg\n", st.electrons.energy[j] * ELECTRON_ENERGY_EV, st.electrons.population[j]);

            for(j = 0; j < st.positrons.size; j++)
                fprintf(positron_temp_file,"%lg\t%lg\n", st.positrons.energy[j] * ELECTRON_ENERGY_EV, st.positrons.population[j]);

            for(j = 0; j < st.protons.size; j++)
                fprintf(hadron_temp_file,"%lg\t%lg\t%lg\n",
                        st.protons.energy[j] * PROTON_ENERGY_EV,
                        st.protons.population[j],
                        st.neutrons.population[j]);

            // Remember that we assume that the mass of all the pions is the same
            for(j = 0; j < st.neutral_pions.size; j++)
                fprintf(pion_temp_file,"%lg\t%lg\t%lg\t%lg\n", st.neutral_pions.energy[j] * NEUTRAL_PION_ENERGY_EV,
                        st.neutral_pions.population[j], st.positive_pions.population[j], st.negative_pions.population[j]);

            // Remember that the mass of all the muons is the same
            for(j = 0; j < st.positive_left_muons.size; j++)
                fprintf(muon_temp_file,"%lg\t%lg\t%lg\n", st.positive_left_muons.energy[j] * MUON_ENERGY_EV,
                        st.positive_left_muons.population[j] + st.positive_right_muons.population[j],
                        st.negative_left_muons.population[j] + st.negative_right_muons.population[j]);

            // Remember that we normalize the energy of the neutrinos with the mass of the electron
            for(j = 0; j < st.electron_neutrinos.size; j++)
                fprintf(neutrino_temp_file,"%lg\t%lg\t%lg\t%lg\t%lg\n", st.electron_neutrinos.energy[j] * ELECTRON_ENERGY_EV,
                        st.electron_neutrinos.population[j],
                        st.electron_antineutrinos.population[j],
                        st.muon_neutrinos.population[j],
                        st.muon_antineutrinos.population[j]);

            fflush(photon_temp_file);
            fflush(lepton_temp_file);
            fflush(positron_temp_file);
            fflush(hadron_temp_file);
            fflush(pion_temp_file);
            fflush(muon_temp_file);
            fflush(neutrino_temp_file);

            nanosleep((struct timespec *) &(struct timespec){0,5e7}, NULL);

            fprintf(gnuplot_pipe, "plot \"%s\" u 1:($3*$3*$4) t \"Photon\", ", photon_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 2:($3*$3*$4) axes x2y1 t \"Photon\", ", photon_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"Electron\", ", lepton_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"Positron\", ", positron_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"Proton\", ", hadron_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:3 t \"Neutron\", ", hadron_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"{/Symbol p}^0\", ", pion_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:3 t \"{/Symbol p}^+\", ", pion_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:4 t \"{/Symbol p}^-\",", pion_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"{/Symbol m}^+\", ", muon_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:3 t \"{/Symbol m}^-\", ", muon_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:2 t \"{/Symbol n}_e\", ", neutrino_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:3 t \"{/Symbol n}_e\", ", neutrino_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:4 t \"{/Symbol n}_{/Symbol m}\",", neutrino_temp_file_filename);
            fprintf(gnuplot_pipe, "     \"%s\" u 1:5 t \"{/Symbol n}_{/Symbol m}\"", neutrino_temp_file_filename);
            fprintf(gnuplot_pipe,"\n");
            fflush(gnuplot_pipe);
            nanosleep((struct timespec *) &(struct timespec){0,5e7}, NULL);

        }
#endif

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

        if(st.t > st.t_max)
            out = true;

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
