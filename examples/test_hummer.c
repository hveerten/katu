#include "state.h"
#include "state_step_common.h"
#include "constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <float.h>

#include <assert.h>

#include <unistd.h>

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    unsigned int i;

    state_t st;
    config_t cfg;
    distribution_metadata_t dm;
    double norm;

    char photon_temp_file_filename[128];
    char hadron_temp_file_filename[128];
    char pion_temp_file_filename[128];

    dm.min = 1e-12;
    dm.max = 1;

    config_read_file(&cfg, "default.toml");
    config_read_file(&cfg, argv[1]);

    state_init_from_config(&st, &cfg);

    state_report_general_info(&st);
 
    // GRB case

    for(i = 0; i < st.protons.size; i++)
    {
        double gp = st.protons.energy[i];
        double np = 1e3 * exp(-pow(gp / 7.356e8, 2)) / (gp*gp);

        np = np > 0 ? np : DBL_MIN;

        st.protons.population[i]     =     np;
        st.protons.log_population[i] = log(np);
    }

    dm.dt = broken_power_law;
    dm.gc = 1.9569e-3;
    dm.p1 = 1;
    dm.p2 = 2;
    norm = distribution_norm(&dm);

    generate_distribution(st.photons.population, st.photons.energy, &dm, st.photons.size);

    for(i = 0; i < st.photons.size; i++)
    {
        double e = st.photons.energy[i];
        double n = st.photons.population[i] / norm;

        if(e < 3.9e-7 || e > 0.587)
        {
            st.photons.population[i]     =     DBL_MIN;
            st.photons.log_population[i] = log(DBL_MIN);
        }
        else
        {
            st.photons.population[i]     =     1e3 * n;
            st.photons.log_population[i] = log(1e3 * n);
        }
    }

    step_calculate_processes(&st);

    sprintf(photon_temp_file_filename, "data/photon_hummer_test_GRB.tsv");
    sprintf(hadron_temp_file_filename, "data/hadron_hummer_test_GRB.tsv");
    sprintf(pion_temp_file_filename,   "data/pion_hummer_test_GRB.tsv");

    state_print_data_to_file(&st, photon,        photon_temp_file_filename);
    state_print_data_to_file(&st, proton,        hadron_temp_file_filename);
    state_print_data_to_file(&st, positive_pion, pion_temp_file_filename);


    // AGN case

    for(i = 0; i < st.protons.size; i++)
    {
        double gp = st.protons.energy[i];
        double np = 1e3 * exp(-pow(gp / 2.6e9, 2)) / (gp*gp);

        np = np > 0 ? np : DBL_MIN;

        st.protons.population[i]     =     np;
        st.protons.log_population[i] = log(np);
    }

    dm.dt = broken_power_law;
    dm.gc = 2.74e-4;
    dm.p1 = 1.6;
    dm.p2 = 1.8;
    norm = distribution_norm(&dm);

    generate_distribution(st.photons.population, st.photons.energy, &dm, st.photons.size);

    for(i = 0; i < st.photons.size; i++)
    {
        double e = st.photons.energy[i];
        double n = st.photons.population[i] / norm;

        if(e < 1.9569e-9 || e > 7e-3)
        {
            st.photons.population[i]     =     DBL_MIN;
            st.photons.log_population[i] = log(DBL_MIN);
        }
        else
        {
            st.photons.population[i]     =     1e3 * n;
            st.photons.log_population[i] = log(1e3 * n);
        }
    }

    step_calculate_processes(&st);

    sprintf(photon_temp_file_filename, "data/photon_hummer_test_AGN.tsv");
    sprintf(hadron_temp_file_filename, "data/hadron_hummer_test_AGN.tsv");
    sprintf(pion_temp_file_filename,   "data/pion_hummer_test_AGN.tsv");

    state_print_data_to_file(&st, photon,        photon_temp_file_filename);
    state_print_data_to_file(&st, proton,        hadron_temp_file_filename);
    state_print_data_to_file(&st, positive_pion, pion_temp_file_filename);


    // BB case

    for(i = 0; i < st.protons.size; i++)
    {
        double gp = st.protons.energy[i];
        double np = 1e3 / (gp*gp);

        if(gp < 1.066e6 || gp > 1.066e12)
        {
            st.protons.population[i]     =     DBL_MIN;
            st.protons.log_population[i] = log(DBL_MIN);
        }
        else
        {
            st.protons.population[i]     =     np;
            st.protons.log_population[i] = log(np);
        }
    }

    dm.dt = black_body;
    dm.t = 1.9569e-5;
    norm = distribution_norm(&dm);

    generate_distribution(st.photons.population, st.photons.energy, &dm, st.photons.size);

    for(i = 0; i < st.photons.size; i++)
    {
        double n = st.photons.population[i] / norm;

        st.photons.population[i]     =     1e3 * n;
        st.photons.log_population[i] = log(1e3 * n);
    }

    step_calculate_processes(&st);

    sprintf(photon_temp_file_filename, "data/photon_hummer_test_BB.tsv");
    sprintf(hadron_temp_file_filename, "data/hadron_hummer_test_BB.tsv");
    sprintf(pion_temp_file_filename,   "data/pion_hummer_test_BB.tsv");

    state_print_data_to_file(&st, photon,        photon_temp_file_filename);
    state_print_data_to_file(&st, proton,        hadron_temp_file_filename);
    state_print_data_to_file(&st, positive_pion, pion_temp_file_filename);

#ifdef USE_THREAD_POOL
    thread_pool_clear(&st.thread_pool);
#endif

    return 0;
}
