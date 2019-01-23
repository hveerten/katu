#include "constants.h"
#include "pion_decay.h"
#include "state.h"

#include <stdio.h>
#include <stdlib.h>

#define MASSES_RATIO ((MUON_MASS * MUON_MASS) / (PION_MASS * PION_MASS))

static double f_positive_right(double x)
{
    if(x < MASSES_RATIO || 1 < x)
        return 0;

    double aux = MASSES_RATIO / ((1 - MASSES_RATIO) * (1 - MASSES_RATIO));

    return aux * (1 - x) / x;
}

static double f_positive_left(double x)
{
    if(x < MASSES_RATIO || 1 < x)
        return 0;

    double aux = 1 / ((1 - MASSES_RATIO) * (1 - MASSES_RATIO));

    return aux * (x - MASSES_RATIO) / x;
}

void charged_pion_decay(state_t *st)
{
    unsigned int i, j;
    double dlng = st->positive_pions.log_energy[1] - st->positive_pions.log_energy[0];

    // All muons with the same size and mass
    for(i = 0; i < st->negative_left_muons.size; i++)
    {
        unsigned int index_base = i * st->positive_pions.size;

        /*double g_muon = st->negative_left_muons.energy[i];*/
        
        double positive_left_muon_production  = 0;
        double positive_right_muon_production = 0;
        double negative_left_muon_production  = 0;
        double negative_right_muon_production = 0;

        // Again, positive pions and negative pions have the same size
        for(j = 0; j < st->positive_pions.size; j++)
        {
            /*double g_pion = st->positive_pions.energy[j];*/
            /*double x = g_muon * MUON_MASS / (g_pion * PION_MASS);*/

            double n_positive_pion = st->positive_pions.population[j];
            double n_negative_pion = st->negative_pions.population[j];

            /*
             *positive_right_muon_production += 2 * n_positive_pion * f_positive_right(x) / g_pion;
             *positive_left_muon_production  += 2 * n_positive_pion * f_positive_left(x)  / g_pion;
             *negative_right_muon_production += 2 * n_negative_pion * f_positive_left(x)  / g_pion;
             *negative_left_muon_production  += 2 * n_negative_pion * f_positive_right(x) / g_pion;
             */

            positive_right_muon_production += 2 * n_positive_pion * st->pion_decay_LUT_positive_right[index_base + j];
            positive_left_muon_production  += 2 * n_positive_pion * st->pion_decay_LUT_positive_left [index_base + j];
            negative_right_muon_production += 2 * n_negative_pion * st->pion_decay_LUT_positive_left [index_base + j];
            negative_left_muon_production  += 2 * n_negative_pion * st->pion_decay_LUT_positive_right[index_base + j];
        }

        st->pion_decay_positive_left_muon_gains[i]  = positive_right_muon_production * dlng / (2 * CHARGED_PION_LIFETIME);
        st->pion_decay_positive_right_muon_gains[i] = positive_left_muon_production  * dlng / (2 * CHARGED_PION_LIFETIME);
        st->pion_decay_negative_left_muon_gains[i]  = negative_right_muon_production * dlng / (2 * CHARGED_PION_LIFETIME);
        st->pion_decay_negative_right_muon_gains[i] = negative_left_muon_production  * dlng / (2 * CHARGED_PION_LIFETIME);
    }

    // All neutrinos with the same mass
    for(i = 0; i < st->muon_neutrinos.size; i++)
    {
        double g_neutrino = st->muon_neutrinos.energy[i];

        double muon_neutrino_production     = 0;
        double muon_antineutrino_production = 0;

        // Again, positive pions and negative pions have the same size
        for(j = 0; j < st->positive_pions.size; j++)
        {
            double g_pion = st->positive_pions.energy[j];
            double x = g_neutrino * ELECTRON_MASS / (g_pion * PION_MASS);

            double n_positive_pion = st->positive_pions.population[j];
            double n_negative_pion = st->negative_pions.population[j];

            if(x < 1 - MASSES_RATIO && x < 1)
            {
                muon_neutrino_production     += 2 * n_positive_pion / (g_pion * g_pion * (1 - MASSES_RATIO)) * (ELECTRON_MASS / PION_MASS);
                muon_antineutrino_production += 2 * n_negative_pion / (g_pion * g_pion * (1 - MASSES_RATIO)) * (ELECTRON_MASS / PION_MASS);
            }
        }

        st->pion_decay_muon_neutrino_gains[i]     = muon_neutrino_production     * dlng / (2 * CHARGED_PION_LIFETIME);
        st->pion_decay_muon_antineutrino_gains[i] = muon_antineutrino_production * dlng / (2 * CHARGED_PION_LIFETIME);
    }
}

/* NOTE:
 * Highly experimental and only works as an integral of the neutral pions!
 */
void neutral_pion_decay(state_t *st)
{
    unsigned int i, j;

    double factor = ELECTRON_MASS / PION_MASS;
    double dlng = st->neutral_pions.log_energy[1] - st->neutral_pions.log_energy[0];

    for(i = 0; i < st->photons.size; i++)
    {
        double e = st->photons.energy[i];

        double aux = 2 * e * factor;

        double g_min = 0;
        if(factor < 1)
            g_min = 1 / (2 * aux);
        else
            g_min = aux / 2;

        g_min = fmin(g_min, 1);

        for(j = 0; j < st->neutral_pions.size && st->neutral_pions.energy[j] < g_min; j++)
        {
        }

        double gains = 0;
        for(     ; j < st->neutral_pions.size; j++)
        {
            gains += st->dt *
                (st->multi_resonances_neutral_pion_gains[i] +
                 st->direct_neutral_pion_gains[j]);
        }
        gains *= dlng;

        st->pion_decay_photon_gains[i] = gains * 4 * M_PI / PION_MASS;
    }
}

void init_pion_decay_LUT_muon_functions(state_t *st)
{
    size_t size = st->negative_left_muons.size * st->positive_pions.size;

    posix_memalign((void **) &st->pion_decay_LUT_positive_right, 32, sizeof(double) * size);
    posix_memalign((void **) &st->pion_decay_LUT_positive_left,  32, sizeof(double) * size);
}

void calculate_pion_decay_LUT_muon_functions(state_t *st)
{
    unsigned int i, j;

    // As usual, we assume same limits and size for all muons and pions
    for(i = 0; i < st->negative_left_muons.size; i++)
    {
        double g_muon = st->negative_left_muons.energy[i];

        unsigned int index_base = i * st->positive_pions.size;

        for(j = 0; j < st->positive_pions.size; j++)
        {
            double g_pion = st->positive_pions.energy[j];
            double x = g_muon * MUON_MASS / (g_pion * PION_MASS);

            st->pion_decay_LUT_positive_right[index_base + j] = f_positive_right(x) / (g_pion * g_pion) * (MUON_MASS / PION_MASS);
            st->pion_decay_LUT_positive_left [index_base + j] = f_positive_left(x)  / (g_pion * g_pion) * (MUON_MASS / PION_MASS);
        }
    }
}
