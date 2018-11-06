#include "pair_production.h"
#include "constants.h"

#include <gsl/gsl_sf_dilog.h>

#include <math.h>
#include <stdlib.h>

double R(double x)
{
    double a = sqrt(1 - 1/x);

    double aux1 = a - 2 * x * a;
    double aux2 = 2 - 2 * x - 1 / x - log(4 * x);
    double aux3 = atanh(a);

    // Happens when 1 << x
    if(a==1) aux3 = log(2) + log(x)/2;

    double aux4 = gsl_sf_dilog((1 - a) / 2) -
                  gsl_sf_dilog((1 + a) / 2);

    return (aux1 - aux2 * aux3 + aux4) / (x * x);
}

void pair_production_process_photon_losses(state_t *st)
{
    unsigned int i, j;

    double factor = 3 * THOMSON_CROSS_SECTION / 4;
    double dlnx = st->photons.log_energy[1] - st->photons.log_energy[0];
    
    for(i = 0; i < st->photons.size; i++)
    {
        double n = st->photons.population[i];
        
        unsigned int index_e_min = st->pair_production_LUT_index_e_min[i];

        double losses = 0;
        for(j = index_e_min; j < st->photons.size; j++)
        {
            unsigned int index_base = i * st->photons.size;
            double nn = st->photons.population[j];

            losses += 2 * nn * st->pair_production_LUT_R[index_base + j];
        }

        losses *= dlnx;

        st->pair_production_losses[i] = -factor * n * losses;
    }
}

void init_pair_production_LUT(state_t *st)
{
    size_t size = st->photons.size * st->photons.size;

    posix_memalign((void **) &st->pair_production_LUT_R, 32, sizeof(double) * size);
    posix_memalign((void **) &st->pair_production_LUT_index_e_min, 32, sizeof(double) * st->photons.size);
}

void calculate_pair_production_LUT(state_t *st)
{
    unsigned int i, j;

    for(i = 0; i < st->photons.size; i++)
    {
        double e = st->photons.energy[i];

        unsigned int index_e_min = 0;

        for(; index_e_min < st->photons.size && st->photons.energy[index_e_min] < 1 / e; index_e_min++)
        {
        }

        st->pair_production_LUT_index_e_min[i] = index_e_min;

        for(j = 0; j < st->photons.size; j++)
        {
            unsigned int index_base = i * st->photons.size;

            double ee = st->photons.energy[j];

            st->pair_production_LUT_R[index_base + j] = R(e * ee) * ee;
        }
    }
}
