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

static double H_0(double e, double ee, double d)
{
    double a = e / (ee * sqrt(ee)) * (e*e/12 - d/8);
    double b = 1 / sqrt(ee) * (e*e*e/6 + e/2 + 1/(4*e));

    return a + b;
}

static double Hpm(double e, double ee, double c, double d)
{
    double aux = sqrt(ee + c * e * e);

    double Ipm = c > 0 ?
                    1 / sqrt(c) * log(e * sqrt(c) + aux) :
                    1 / sqrt(-c) * asin(e * sqrt(-c/ee));

    double aux1 = - e / (8 * aux) * (d / ee + 2 / c);
    double aux2 = (2 - (ee - 1) / c) * Ipm / 4;
    double aux3 = aux / 4 * (e / c + 1 / (e * ee));

    return aux1 + aux2 + aux3;
}

static double rate_lepton_gains(double e1, double e2, double g)
{
    double ee = e1 * e2;
    double E  = e1 + e2;

    double cp = pow(e1 - g, 2) - 1;
    double cm = pow(e2 - g, 2) - 1;
    double dp = e1 * e1 + ee + g * (e2 - e1);
    double dm = e2 * e2 + ee - g * (e2 - e1);

    double aux = g * (E - g) + 1;
    double e_min_aux = sqrt((aux - sqrt(aux * aux - E * E)) / 2);
    double e_max_aux = sqrt((aux + sqrt(aux * aux - E * E)) / 2);

    double e_min = fmax(1            , e_min_aux);
    double e_max = fmin(sqrt(e1 * e2), e_max_aux);


    double aux1 = sqrt(E*E - 4 * e_max * e_max) / 4;
    double aux2 = cp != 0 ? Hpm(e_max, ee, cp, dp) : H_0(e_max, ee, dp);
    double aux3 = cm != 0 ? Hpm(e_max, ee, cm, dm) : H_0(e_max, ee, dm);

    double aux4 = sqrt(E*E - 4 * e_min * e_min) / 4;
    double aux5 = cp != 0 ? Hpm(e_min, ee, cp, dp) : H_0(e_min, ee, dp);
    double aux6 = cm != 0 ? Hpm(e_min, ee, cm, dm) : H_0(e_min, ee, dm);

    return (aux1 + aux2 + aux3 - aux4 - aux5 - aux6) / ee;
}

void pair_production_process_photon_losses(state_t *st)
{
    unsigned int i, j;

    double factor = 3 * THOMSON_CROSS_SECTION * M_PI * LIGHT_SPEED / 8;
    double dlnx = st->photons.log_energy[1] - st->photons.log_energy[0];

    for(i = 0; i < st->photons.size; i++)
    {
        double n = st->photons.population[i];

        unsigned int index_base = i * st->photons.size;
        unsigned int index_e_min = st->pair_production_LUT_photon_losses_index_e_min[i];
        unsigned int index_e_max = st->photons.size - 1;

        if(index_e_min >= index_e_max) continue;

        double losses = 0;

        losses += st->photons.population[index_e_min] * st->pair_production_LUT_photon_losses_R[index_base + index_e_min];
        losses += st->photons.population[index_e_max] * st->pair_production_LUT_photon_losses_R[index_base + index_e_max];

        for(j = index_e_min + 1; j < index_e_max; j++)
        {
            double nn = st->photons.population[j];

            losses += 2 * nn * st->pair_production_LUT_photon_losses_R[index_base + j];
        }

        losses *= dlnx / 2;

        st->pair_production_photon_losses[i] = -factor * n * losses;
    }
}

void pair_production_process_lepton_gains(state_t *st)
{
    unsigned int i, j, k;

    double factor = 3 * THOMSON_CROSS_SECTION * LIGHT_SPEED / 4;
    double dlnx = st->photons.log_energy[1] - st->photons.log_energy[0];

    double gains_inner[st->photons.size];

    for(i = 0; i < st->electrons.size; i++)
    {
        unsigned int index_base1 = i * st->photons.size;

        for(j = 0; j < st->photons.size; j++)
        {
            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            unsigned int index_e_min = st->pair_production_LUT_lepton_gains_index_e_min[index_base1 + j];
            unsigned int index_e_max = st->photons.size - 1;

            if(index_e_min >= index_e_max) continue;

            gains_inner[j] = 0;
            gains_inner[j] += st->photons.population[index_e_min] * st->pair_production_LUT_lepton_gains_R[index_base2 + index_e_min];
            gains_inner[j] += st->photons.population[index_e_max] * st->pair_production_LUT_lepton_gains_R[index_base2 + index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double n = st->photons.population[k];

                gains_inner[j] += 2 * n * st->pair_production_LUT_lepton_gains_R[index_base2 + k];
            }

            gains_inner[j] *= dlnx / 2;
        }

        double gains = 0;
        gains += st->photons.population[0]                    * gains_inner[0];
        gains += st->photons.population[st->photons.size - 1] * gains_inner[st->photons.size - 1];

        for(j = 1; j < st->photons.size - 1; j++)
        {
            gains += 2 * st->photons.population[j] * gains_inner[j];
        }
        gains *= dlnx / 2;

        st->pair_production_electron_gains[i] = factor * gains;
        st->pair_production_positron_gains[i] = factor * gains;
    }
}

void init_pair_production_LUT_photon_losses(state_t *st)
{
    size_t size = st->photons.size * st->photons.size;

    posix_memalign((void **) &st->pair_production_LUT_photon_losses_R, 32, sizeof(double) * size);
    posix_memalign((void **) &st->pair_production_LUT_photon_losses_index_e_min, 32, sizeof(double) * st->photons.size);
}

void init_pair_production_LUT_lepton_gains(state_t *st)
{
    size_t size1 = st->electrons.size * st->photons.size;
    size_t size2 = st->electrons.size * st->photons.size * st->photons.size;

    posix_memalign((void **) &st->pair_production_LUT_lepton_gains_R, 32, sizeof(double) * size2);
    posix_memalign((void **) &st->pair_production_LUT_lepton_gains_index_e_min, 32, sizeof(double) * size1);
}

void calculate_pair_production_LUT_photon_losses(state_t *st)
{
    unsigned int i, j;

    for(i = 0; i < st->photons.size; i++)
    {
        double e = st->photons.energy[i];

        unsigned int index_e_min = 0;

        for(; index_e_min < st->photons.size && st->photons.energy[index_e_min] < 1 / e; index_e_min++)
        {
        }

        st->pair_production_LUT_photon_losses_index_e_min[i] = index_e_min;

        for(j = 0; j < st->photons.size; j++)
        {
            unsigned int index_base = i * st->photons.size;

            double ee = st->photons.energy[j];

            st->pair_production_LUT_photon_losses_R[index_base + j] = R(e * ee) * ee;
        }
    }
}

void calculate_pair_production_LUT_lepton_gains(state_t *st)
{
    unsigned int i, j, k;

    for(i = 0; i < st->electrons.size; i++)
    {
        double g = st->electrons.energy[i];

        unsigned int index_base1 = i * st->photons.size;

        for(j = 0; j < st->photons.size; j++)
        {
            double e1 = st->photons.energy[j];

            unsigned int index_e_min = 0;
            double e_min = fmax(1 / e1, g + 1 - e1);

            for(; index_e_min < st->photons.size && st->photons.energy[index_e_min] < e_min; index_e_min++)
            {
            }

            st->pair_production_LUT_lepton_gains_index_e_min[index_base1 + j] = index_e_min;

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            for(k = 0; k < st->photons.size; k++)
            {
                unsigned int index = index_base2 + k;

                double e2 = st->photons.energy[k];

                st->pair_production_LUT_lepton_gains_R[index] = rate_lepton_gains(e1, e2, g);
            }
        }
    }
}
