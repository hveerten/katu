#include "pair_annihilation.h"
#include "constants.h"

#include "utils.h"

#include <gsl/gsl_sf_dilog.h>

#include <math.h>
#include <stdlib.h>

#include <stdio.h>
#include <assert.h>

static double R(double x)
{
    (void) x;
    return 0;
}

static double H_0(double g, double g2, double g_aux, double g_aux2, double d)
{
    double a = 1/ g_aux           * (2./3. * g * g2 + 2 * g + 1 / g);
    double b = 1/(g_aux * g_aux2) * (2./3. * g * g2 - d * g);

    return a + b;
}

static double H_pm(double g, double g2, double g_aux, double g_aux2, double ge, double d)
{
    double c = ge*ge - 1;
    /*double u = sqrt(c*g2 + g_aux2);*/
    double u = sqrt(ge*ge*g2 + g_aux2 - g2);

    double I = c > 0 ? 1 / sqrt( c) *  log(g * sqrt( c) + u) :
                       1 / sqrt(-c) * asin(g * sqrt(-c) / g_aux);

    double aux1 = I * (2 + (1 - g_aux2) / c);
    double aux2 = (1 / g - g / c + g * (2 * c - d) / (2 * g_aux2)) / u;
    double aux3 = g * u / c;

    return aux1 + aux2 + aux3;
}

static double rate_photon_gains(double g1, double g2, double e)
{
    double b1 = sqrt(1 - 1/(g1*g1));
    double b2 = sqrt(1 - 1/(g2*g2));

    double g_max2_aux = (1 + g1*g2 * (1 + b1*b2)) / 2;
    double g_min2_aux = (1 + g1*g2 * (1 - b1*b2)) / 2;

    /*double cp = pow(g2 - e, 2) - 1;*/
    /*double cm = pow(g1 - e, 2) - 1;*/
    double dp = g2 * (g1 + g2) + e * (g1 - g2);
    double dm = g1 * (g1 + g2) - e * (g1 - g2);

    double g_aux2 = e * (g1 + g2 - e);

    double g_max2 = fmin(g_max2_aux, g_aux2);
    double g_min2 =      g_min2_aux;            // is this always > 1?

    // protection from g_aux2 + c * g_max2 = 0
    if(g_max2 == g_aux2 && ((g1 == e) || (g2 == e)))
        g_max2 *= (1 - 1e-15);

    double g_max = sqrt(g_max2);
    double g_min = sqrt(g_min2);
    double g_aux = sqrt(g_aux2);

    /*fprintf(stderr,"g_max %lg\n", g_max);*/
    /*fprintf(stderr,"g_min %lg\n", g_min);*/

    double aux1_max = sqrt(pow(g1 + g2, 2) - 4 * g_max2);
    /*double aux1_min = sqrt(pow(g1 + g2, 2) - 4 * g_min2);*/
    double aux1_min = fabs(b1*g1 - b2*g2);

    double aux2_max = fabs(g2 - e) != 1 ? H_pm(g_max, g_max2, g_aux, g_aux2, g2 - e, dp) : H_0(g_max, g_max2, g_aux, g_aux2, dp);
    double aux2_min = fabs(g2 - e) != 1 ? H_pm(g_min, g_min2, g_aux, g_aux2, g2 - e, dp) : H_0(g_min, g_min2, g_aux, g_aux2, dp);

    double aux3_max = fabs(g1 - e) != 1 ? H_pm(g_max, g_max2, g_aux, g_aux2, g1 - e, dm) : H_0(g_max, g_max2, g_aux, g_aux2, dm);
    double aux3_min = fabs(g1 - e) != 1 ? H_pm(g_min, g_min2, g_aux, g_aux2, g1 - e, dm) : H_0(g_min, g_min2, g_aux, g_aux2, dm);

    if(g_max2 == g_aux2)
        aux1_max = fabs(b1*g1 + b2*g2);

    /*fprintf(stderr,"1 %lg - %lg %lg\n", aux1_max, aux1_min, aux1_max - aux1_min);*/
    /*fprintf(stderr,"2 %lg - %lg %lg\n", aux2_max, aux2_min, aux2_max - aux2_min);*/
    /*fprintf(stderr,"3 %lg - %lg %lg\n", aux3_max, aux3_min, aux3_max - aux3_min);*/

    return (aux1_max - aux1_min + aux2_max - aux2_min + aux3_max - aux3_min) / (b1*b2*g1*g2);
}

void pair_annihilation_process_photon_gains(state_t *st)
{
    unsigned int i, j, k;

    double factor = 3 * THOMSON_CROSS_SECTION * LIGHT_SPEED / 8;
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    double gains_inner[st->electrons.size];

    for(i = 0; i < st->photons.size; i++)
    {
        unsigned int index_base1 = i * st->photons.size;

        unsigned int index_g1_min = st->pair_annihilation_LUT_photon_gains_index_g1_min[index_base1 + j];
        unsigned int index_g1_max = st->electrons.size - 1;

        for(j = index_g1_min; j < st->electrons.size; j++)
        {
            unsigned int index_base2 = (index_base1 + j) * st->electrons.size;

            unsigned int index_g2_min = st->pair_annihilation_LUT_photon_gains_index_g1_min[index_base1 + j];
            unsigned int index_g2_max = st->positrons.size - 1;

            if(index_g2_min >= index_g2_max) continue;

            gains_inner[j] = 0;
            gains_inner[j] += st->positrons.population[index_g2_min] * st->pair_annihilation_LUT_photon_gains_R[index_base2 + index_g2_min];
            gains_inner[j] += st->positrons.population[index_g2_max] * st->pair_annihilation_LUT_photon_gains_R[index_base2 + index_g2_max];

            for(k = index_g2_min + 1; k < index_g2_max; k++)
            {
                double n = st->positrons.population[k];

                gains_inner[j] += 2 * n * st->pair_annihilation_LUT_photon_gains_R[index_base2 + k];
            }

            gains_inner[j] *= dlng / 2;
        }

        double gains = 0;
        gains += st->electrons.population[index_g1_min] * gains_inner[index_g1_min];
        gains += st->electrons.population[index_g1_max] * gains_inner[index_g1_max];

        for(j = index_g1_min + 1; j < st->electrons.size - 1; j++)
            gains += 2 * st->electrons.population[j] * gains_inner[j];

        gains *= dlng / 2;

        st->pair_annihilation_photon_gains[i] = factor * gains;
    }
}

void pair_annihilation_process_lepton_losses(state_t *st)
{
    (void) st;
}

void init_pair_annihilation_LUT_photon_gains(state_t *st)
{
    size_t size0 = st->photons.size;
    size_t size1 = st->photons.size * st->electrons.size;
    size_t size2 = st->photons.size * st->electrons.size * st->positrons.size;

    MEM_PREPARE(st->pair_annihilation_LUT_photon_gains_R, size2, double);
    MEM_PREPARE(st->pair_annihilation_LUT_photon_gains_index_g1_min, size0, unsigned int);
    MEM_PREPARE(st->pair_annihilation_LUT_photon_gains_index_g2_min, size1, unsigned int);
}

void init_pair_annihilation_LUT_lepton_losses(state_t *st)
{
    size_t size = st->electrons.size * st->positrons.size;

    MEM_PREPARE(st->pair_annihilation_LUT_lepton_losses_R, size, double);
    MEM_PREPARE(st->pair_annihilation_LUT_lepton_losses_index_g_min, size, unsigned int);
}

void calculate_pair_annihilation_LUT_photon_gains(state_t *st)
{
    unsigned int i, j, k;

    for(i = 0; i < st->photons.size; i++)
    {
        unsigned int index_base1 = i * st->photons.size;

        double e = st->photons.energy[i];

        double ga = (4*e*e+1) / (4*e);
        double gb = (2*e*e - 2*e + 1) / (2*e - 1);

        double g1_min = e <= 0.5 ? ga : 1;
        unsigned int index_g1_min = 0;

        for(index_g1_min = 0;
                index_g1_min < st->electrons.size &&
                st->electrons.energy[index_g1_min] < g1_min;
            index_g1_min++)
        { }

        st->pair_annihilation_LUT_photon_gains_index_g1_min[i] = index_g1_min;

        for(j = index_g1_min; j < st->electrons.size; j++)
        {
            double g1 = st->electrons.energy[j];
            double b1 = sqrt(1 - 1/(g1*g1));

            double F_plus  = 2 * e - g1 * (1 + b1);
            double F_minus = 2 * e - g1 * (1 - b1);

            double g2_min = 1;
            unsigned int index_g2_min = 0;

            if     (g1 >= gb && e >= 0.5) g2_min = 1;
            else if(g1 <  gb && e >= 1.0) g2_min = 0.5 * (F_plus  + 1/F_plus);
            else                          g2_min = 0.5 * (F_minus + 1/F_minus);

            for(index_g2_min = 0;
                    index_g2_min < st->electrons.size &&
                    st->electrons.energy[index_g2_min] < g2_min;
                index_g2_min++)
            { }

            st->pair_annihilation_LUT_photon_gains_index_g1_min[index_base1 + j] = index_g2_min;

            unsigned int index_base2 = (index_base1 + j) * st->electrons.size;

            for(k = index_g2_min; k < st->electrons.size; k++)
            {
                unsigned int index = index_base2 + k;

                double g2 = st->electrons.energy[k];

                double r = rate_photon_gains(g1, g2, e);

                st->pair_annihilation_LUT_photon_gains_R[index] = r;

                if(!isnormal(r))
                {
                    fprintf(stderr,"%u %u %u %lg\n", i, j, k, r);
                    break;
                }
            }
        }
    }
}

void calculate_pair_annihilation_LUT_lepton_losses(state_t *st)
{
}
