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

    if(g2 == g_aux2)
        I = c > 0 ? 1 / sqrt( c) * (0.5 * log(g_aux2) + asinh(sqrt( c))) :
                    1 / sqrt(-c) *                       asin(sqrt(-c));

    double aux1 = I * (2 + (1 - g_aux2) / c);
    double aux2 = (1 / g - g / c + g * (2 * c - d) / (2 * g_aux2)) / u;
    double aux3 = g * u / c;

    aux2 = g/u * (g2 + 1/g2 + (g_aux2-1)/c + (2 * c - d) / (2 * g_aux2));
    aux3 = 0;

    return aux1 + aux2 + aux3;
}

static double J_pm(double g_min2, double g_max2, double g_aux2, double ge, double d, double gg)
{
    double c = ge*ge - 1;

    double u_min = sqrt(ge*ge*g_min2 + g_aux2 - g_min2);
    double u_max = sqrt(ge*ge*g_max2 + g_aux2 - g_max2);

    double g_min = sqrt(g_min2);
    double g_max = sqrt(g_max2);
    double g_aux = sqrt(g_aux2);

    double gu_min = g_min / u_min;
    double gu_max = g_max / u_max;

    double I_min = c > 0 ? 1 / sqrt( c) *  log(g_min * sqrt( c) + u_min) :
                           1 / sqrt(-c) * asin(g_min * sqrt(-c) / g_aux);
    double I_max = c > 0 ? 1 / sqrt( c) *  log(g_max * sqrt( c) + u_max) :
                           1 / sqrt(-c) * asin(g_max * sqrt(-c) / g_aux);

    double aux_max = g_aux2 / (c * g_max2);
    double aux_min = g_aux2 / (c * g_min2);

    double I_max_min;

    if(g_max2 == g_aux2)
        gu_max = 1 / fabs(ge);
    if(g_max2 == g_aux2 && ge == 0)
        gu_max = 1e200;

    if(g_max2 == g_aux2)
        aux_max = 1 / c;

    if(c > 0)
        I_max_min = (0.5 * log(g_max2 / g_min2) + log((1 + sqrt(1 + aux_max)) / (1 + sqrt(1 + aux_min)))) / sqrt(c);
    else
    {
        if(g_max2 == g_aux2)
            I_max_min = (asin(        sqrt(-c)        ) - asin(g_min * sqrt(-c) / g_aux)) / sqrt(-c);
        else
            I_max_min = (asin(g_max * sqrt(-c) / g_aux) - asin(g_min * sqrt(-c) / g_aux)) / sqrt(-c);
    }

    if(fabs(aux_min) < 1e-5)
        gu_min = 1 / (sqrt(c) * (1 + aux_min / 2));


    double aux01 = (2 + (1 - g_aux2) / c);
    double aux02 = I_max_min;

    double aux11 = gu_max * (g_max2 + 1/g_max2 - 0.5);
    double aux12 = gu_min * (g_min2 + 1/g_min2 - 0.5);

    double aux21 = gu_max - gu_min;
    double aux22 = (g_aux2 - 1) / c - 1 / (2 * g_aux2);

    double aux3 = aux21 * (c - gg)/(2 * g_aux2);

    if(fabs(aux_min) < 1e-5 && fabs(aux_max) < 1e5)
        aux3 = 1/(4*sqrt(c)) * (g_min2 - g_max2) / (g_min2*g_max2) / (sqrt(1 + aux_max) * sqrt(1 + aux_min)) * (1 - gg/c);

    if(gu_max - gu_min < 1e-5)
        aux21 = (g_max*u_min - g_min*u_max) / (u_min * u_max);

    return aux01 * aux02 + aux11 - aux12 + aux21 * aux22 + aux3;
}

static double rate_photon_gains(double g1, double g2, double e)
{
    double b1 = sqrt(1 - 1/(g1*g1));
    double b2 = sqrt(1 - 1/(g2*g2));

    double g_max2_aux = (1 + g1*g2 * (1 + b1*b2)) / 2;
    double g_min2_aux = (1 + g1*g2 * (1 - b1*b2)) / 2;

    if(b1*b2 > 1-1e-10)
        g_min2_aux = (1 + 0.5*(g1/g2 + g2/g1))/2;

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

    if(g_min > g_max) return 0;

    double aux1_max = sqrt(pow(g1 + g2, 2) - 4 * g_max2);
    /*double aux1_min = sqrt(pow(g1 + g2, 2) - 4 * g_min2);*/
    double aux1_min = b1*g1 + b2*g2;

    double aux2_max = fabs(g2 - e) != 1 ? H_pm(g_max, g_max2, g_aux, g_aux2, g2 - e, dp) : H_0(g_max, g_max2, g_aux, g_aux2, dp);
    double aux2_min = fabs(g2 - e) != 1 ? H_pm(g_min, g_min2, g_aux, g_aux2, g2 - e, dp) : H_0(g_min, g_min2, g_aux, g_aux2, dp);

    double aux3_max = fabs(g1 - e) != 1 ? H_pm(g_max, g_max2, g_aux, g_aux2, g1 - e, dm) : H_0(g_max, g_max2, g_aux, g_aux2, dm);
    double aux3_min = fabs(g1 - e) != 1 ? H_pm(g_min, g_min2, g_aux, g_aux2, g1 - e, dm) : H_0(g_min, g_min2, g_aux, g_aux2, dm);

    if(fabs(g2 - e) != 1)
    {
        aux2_max = J_pm(g_min2, g_max2, g_aux2, g2 - e, dp, g1*g2);
        aux2_min = 0;
    }

    if(fabs(g1 - e) != 1)
    {
        aux3_max = J_pm(g_min2, g_max2, g_aux2, g1 - e, dm, g1*g2);
        aux3_min = 0;
    }

    if(g_max2 == g_max2_aux)
    {
        /*aux1_max = fabs(b1*g1 - b2*g2);*/

        if(b1*g1 > b2*g2)
        {
            aux1_max = -2 * b2*g2;
            aux1_min = 0;
        }
        else
        {
            aux1_max = -2 * b1*g1;
            aux1_min = 0;
        }
    }
    else if(4*g_aux2 / ((g1+g2)*(g1+g2)) < 1e-5)
    {
        double aa = g1 > 1e5 ? 1 / (2*g1) : g1 * (1 - b1);
        double bb = g2 > 1e5 ? 1 / (2*g2) : g2 * (1 - b2);

        aux1_max = aa + bb;
        aux1_min = 2 * g_aux2 / (g1 + g2);
    }

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
        unsigned int index_base1 = i * st->electrons.size;

        unsigned int index_g1_min = st->pair_annihilation_LUT_photon_gains_index_g1_min[i];
        unsigned int index_g1_max = st->electrons.size - 1;

        if(index_g1_min >= index_g1_max) continue;

        for(j = index_g1_min; j < st->electrons.size; j++)
        {
            unsigned int index_base2 = (index_base1 + j) * st->electrons.size;

            unsigned int index_g2_min = st->pair_annihilation_LUT_photon_gains_index_g2_min[index_base1 + j];
            unsigned int index_g2_max = st->positrons.size - 1;

            gains_inner[j] = 0;

            if(index_g2_min >= index_g2_max) continue;

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
        unsigned int index_base1 = i * st->electrons.size;

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

            if(e == g1)
                F_plus = e * (1 - b1);
            if(e == g1 && g1 > 1e5)
                F_plus = 1 / (2 * e);

            if(g1 > 1e5)
                F_minus = 2*e - 0.5 * (1/g1 + 0.25*1/(g1*g1*g1));

            double g2_min = 1;
            unsigned int index_g2_min = 0;

            if     (g1 >= gb && e >= 0.5) g2_min = 1;
            else if(g1 <  gb && e >= 1.0) g2_min = 0.5 * (F_plus  + 1/F_plus);
            else                          g2_min = 0.5 * (F_minus + 1/F_minus);

            /* NOTE: Be careful with the case that g2_min ~ a value of the grid */
            g2_min *= 1 + 1e-7;

            for(index_g2_min = 0;
                    index_g2_min < st->electrons.size &&
                    st->electrons.energy[index_g2_min] < g2_min;
                index_g2_min++)
            { }

            st->pair_annihilation_LUT_photon_gains_index_g2_min[index_base1 + j] = index_g2_min;

            unsigned int index_base2 = (index_base1 + j) * st->electrons.size;

            for(k = index_g2_min; k < st->electrons.size; k++)
            {
                unsigned int index = index_base2 + k;

                double g2 = st->electrons.energy[k];

                double r = rate_photon_gains(g1, g2, e);

                st->pair_annihilation_LUT_photon_gains_R[index] = r;

                if(!isfinite(r) || r < 0)
                {
                    fprintf(stderr,"Error in Pair Annihilation LUT: %u %u %u %lg\n", i, j, k, r);
                }
            }
        }
    }
}

void calculate_pair_annihilation_LUT_lepton_losses(state_t *st)
{
}
