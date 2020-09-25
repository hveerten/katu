#include "constants.h"
#include "pion_decay.h"
#include "state.h"

#include <stdlib.h>
#include <stdio.h>

static double f_positive_muon(double x, double h)
{
    if(1 < x)
        return 0;

    if(fabs(h) > 1)
        return 0;

    if(h ==  1.0) return 4./3                      - 4 * x * x * x / 3;
    /*if(h == -1.0) return 2             - 6 * x * x + 4 * x * x * x;*/

    if(h == -1.0) return 2 * (1 - x*x * (3 - 2 * x));

    double aux1 = 5./3 - 3 * x * x + 4 * x*x*x / 3;
    double aux2 =-1./3 + 3 * x * x - 8 * x*x*x / 3;

    return aux1 + h * aux2;
}

static double f_positive_electron(double x, double h)
{
    if(1 < x)
        return 0;

    if(fabs(h) > 1)
        return 0;

    /*if(h ==  1.0) return 4 - 12 * x + 12 * x * x -  4 * x * x * x;*/
    if(h ==  1.0) return 4 * (1 - 3 * (x - x * x) - x * x * x);
    if(h == -1.0) return     12 * x - 24 * x * x + 12 * x * x * x;

    double aux1 = 2          -  6 * x * x + 4 * x * x * x;
    double aux2 = 2 - 12 * x + 18 * x * x - 8 * x * x * x;

    return aux1 + h * aux2;
}

void muon_decay(state_t *st)
{
    unsigned int i, j;
    double dlng = st->negative_left_muons.log_energy[1] - st->negative_left_muons.log_energy[0];

    // All muons with the same size and mass
    for(i = 0; i < st->electrons.size; i++)
    {
        double g_electron = st->electrons.energy[i];

        double electron_production = 0;
        double positron_production = 0;

        // Again, positive pions and negative pions have the same size
        for(j = 0; j < st->negative_left_muons.size; j++)
        {
            double g_muon = st->negative_left_muons.energy[j];
            double x = g_electron * ELECTRON_MASS / (g_muon * MUON_MASS);

            double n_positive_left_muons  = st->positive_left_muons.population[j];
            double n_positive_right_muons = st->positive_right_muons.population[j];
            double n_negative_left_muons  = st->negative_left_muons.population[j];
            double n_negative_right_muons = st->negative_right_muons.population[j];

            positron_production += 2 * n_positive_left_muons  * f_positive_muon(x, -1) / g_muon;
            positron_production += 2 * n_positive_right_muons * f_positive_muon(x,  1) / g_muon;
            electron_production += 2 * n_negative_left_muons  * f_positive_muon(x,  1) / g_muon;
            electron_production += 2 * n_negative_right_muons * f_positive_muon(x, -1) / g_muon;
        }

        st->muon_decay_electron_gains[i] = electron_production * dlng / (2 * MUON_LIFETIME);
        st->muon_decay_positron_gains[i] = positron_production * dlng / (2 * MUON_LIFETIME);
    }

    // All neutrinos with the same mass
    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        unsigned int index_base = i * st->negative_left_muons.size;

        double g_neutrino = st->muon_neutrinos.energy[i];

        double electron_neutrino_production     = 0;
        double electron_antineutrino_production = 0;
        double muon_neutrino_production         = 0;
        double muon_antineutrino_production     = 0;

        // Again, positive pions and negative pions have the same size
        for(j = 0; j < st->negative_left_muons.size; j++)
        {
            /*double g_muon = st->negative_left_muons.energy[j];*/
            /*double x = g_neutrino * ELECTRON_MASS / (g_muon * MUON_MASS);*/

            double n_positive_left_muons  = st->positive_left_muons.population[j];
            double n_positive_right_muons = st->positive_right_muons.population[j];
            double n_negative_left_muons  = st->negative_left_muons.population[j];
            double n_negative_right_muons = st->negative_right_muons.population[j];

            /*
             *electron_neutrino_production     += 2 * n_positive_left_muons  * f_positive_electron(x, -1) / g_muon;
             *electron_neutrino_production     += 2 * n_positive_right_muons * f_positive_electron(x,  1) / g_muon;
             *electron_antineutrino_production += 2 * n_negative_left_muons  * f_positive_electron(x,  1) / g_muon;
             *electron_antineutrino_production += 2 * n_negative_right_muons * f_positive_electron(x, -1) / g_muon;
             *muon_neutrino_production         += 2 * n_negative_left_muons  * f_positive_muon(x,  1)     / g_muon;
             *muon_neutrino_production         += 2 * n_negative_right_muons * f_positive_muon(x, -1)     / g_muon;
             *muon_antineutrino_production     += 2 * n_positive_left_muons  * f_positive_muon(x, -1)     / g_muon;
             *muon_antineutrino_production     += 2 * n_positive_right_muons * f_positive_muon(x,  1)     / g_muon;
             */

            electron_neutrino_production     += 2 * n_positive_left_muons  * st->muon_decay_LUT_positive_electron_minus_1[index_base + j];
            electron_neutrino_production     += 2 * n_positive_right_muons * st->muon_decay_LUT_positive_electron_plus_1 [index_base + j];
            electron_antineutrino_production += 2 * n_negative_left_muons  * st->muon_decay_LUT_positive_electron_plus_1 [index_base + j];
            electron_antineutrino_production += 2 * n_negative_right_muons * st->muon_decay_LUT_positive_electron_minus_1[index_base + j];
            muon_neutrino_production         += 2 * n_negative_left_muons  * st->muon_decay_LUT_positive_muon_plus_1 [index_base + j];
            muon_neutrino_production         += 2 * n_negative_right_muons * st->muon_decay_LUT_positive_muon_minus_1[index_base + j];
            muon_antineutrino_production     += 2 * n_positive_left_muons  * st->muon_decay_LUT_positive_muon_minus_1[index_base + j];
            muon_antineutrino_production     += 2 * n_positive_right_muons * st->muon_decay_LUT_positive_muon_plus_1 [index_base + j];
        }

        st->muon_decay_electron_neutrino_gains[i]     = electron_neutrino_production     * dlng / (2 * MUON_LIFETIME);
        st->muon_decay_electron_antineutrino_gains[i] = electron_antineutrino_production * dlng / (2 * MUON_LIFETIME);
        st->muon_decay_muon_neutrino_gains[i]         = muon_neutrino_production         * dlng / (2 * MUON_LIFETIME);
        st->muon_decay_muon_antineutrino_gains[i]     = muon_antineutrino_production     * dlng / (2 * MUON_LIFETIME);
    }
}

void init_muon_decay_LUT_neutrino_functions(state_t *st)
{
    size_t size = st->muon_neutrinos.size * st->positive_left_muons.size;

    posix_memalign((void **) &st->muon_decay_LUT_positive_electron_minus_1, 32, sizeof(double)* size);
    posix_memalign((void **) &st->muon_decay_LUT_positive_electron_plus_1,  32, sizeof(double)* size);
    posix_memalign((void **) &st->muon_decay_LUT_positive_muon_minus_1,     32, sizeof(double)* size);
    posix_memalign((void **) &st->muon_decay_LUT_positive_muon_plus_1,      32, sizeof(double)* size);
}

void calculate_muon_decay_LUT_neutrino_functions(state_t *st)
{
    unsigned int i, j;

    // As usual, we assume same limits and size for all neutrinos and muons
    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        double g_neutrino = st->muon_neutrinos.energy[i];

        unsigned int index_base = i * st->negative_left_muons.size;

        for(j = 0; j < st->negative_left_muons.size; j++)
        {
            double g_muon = st->negative_left_muons.energy[j];
            double x = g_neutrino * ELECTRON_MASS / (g_muon * MUON_MASS);

            st->muon_decay_LUT_positive_electron_minus_1[index_base + j] = f_positive_electron(x, -1) / g_muon * (ELECTRON_MASS / MUON_MASS);
            st->muon_decay_LUT_positive_electron_plus_1 [index_base + j] = f_positive_electron(x,  1) / g_muon * (ELECTRON_MASS / MUON_MASS);
            st->muon_decay_LUT_positive_muon_plus_1     [index_base + j] = f_positive_muon(x,  1)     / g_muon * (ELECTRON_MASS / MUON_MASS);
            st->muon_decay_LUT_positive_muon_minus_1    [index_base + j] = f_positive_muon(x, -1)     / g_muon * (ELECTRON_MASS / MUON_MASS);

            if(f_positive_electron(x, -1) < 0) fprintf(stderr,"ERROR generating tables for muon decay +e -1\n");
            if(f_positive_electron(x,  1) < 0) fprintf(stderr,"ERROR generating tables for muon decay +e  1\n");
            if(f_positive_muon(x, -1) < 0)     fprintf(stderr,"ERROR generating tables for muon decay +m -1\n");
            if(f_positive_muon(x,  1) < 0)     fprintf(stderr,"ERROR generating tables for muon decay +m  1\n");
        }
    }
}
