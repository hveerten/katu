#include "pair_production.h"
#include "constants.h"

#include <gsl/gsl_sf_dilog.h>

#include <math.h>
#include <stdlib.h>

static double R(double x)
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

static double Jpm(double e_min2, double e_max2, double ee, double eg, double d)
{
    double c = eg * eg - 1;
    double e_min = sqrt(e_min2);
    double e_max = sqrt(e_max2);

    double A_min = ee - e_min2 + eg * eg * e_min2;
    double A_max = ee - e_max2 + eg * eg * e_max2;

    double B_min = sqrt(c) * sqrt(1 + ee/e_min2/c);
    double B_max = sqrt(c) * sqrt(1 + ee/e_max2/c);

    double I_min = c > 0 ?
                    1 / sqrt(c) * log(e_min * sqrt(c) + sqrt(A_min)) :
                    1 / sqrt(-c) * asin(e_min * sqrt(-c/ee));
    double I_max = c > 0 ?
                    1 / sqrt(c) * log(e_max * sqrt(c) + sqrt(A_max)) :
                    1 / sqrt(-c) * asin(e_max * sqrt(-c/ee));

    if(c < 0)
    {
        B_min = sqrt(A_min) / e_min;
        B_max = sqrt(A_max) / e_max;
    }
    if(c == -1)
        B_max = sqrt(ee - e_max2) / e_max;

    if(ee == e_max2)
    {
        I_max = c > 0 ?
                    1 / sqrt(c) * (log(ee)/2 + asinh(sqrt(c))) :
                    1 / sqrt(-c) * asin(sqrt(-c));

        B_max = sqrt(1 + c);
    }


    double aux1 = (e_max2 + 1 / e_max2) / B_max -
                  (e_min2 + 1 / e_min2) / B_min;

    double aux2_factor = ((ee - 1)/c + (2*c-d)/(2*ee));

    double aux2 = aux2_factor * (B_min - B_max) / (B_min * B_max);

    if(c > 0)
    {
        if(ee/(e_min2*c) < 1e-5 && ee/(e_max2*c) < 1e-5)
            aux2 = aux2_factor * (ee/(2*c) * (1/e_min2 - 1/e_max2)) / (B_min * sqrt(1 + ee/e_max2/c));
    }

    double aux3 = (I_max - I_min) * (2 - (ee - 1) / c);

    /*fprintf(stderr,"Jpm\n");*/
    /*fprintf(stderr,"A  %lg %lg\n", A_min, A_max);*/
    /*fprintf(stderr,"B  %lg %lg\n", B_min, B_max);*/
    /*fprintf(stderr,"c  %lg\n", c);*/
    /*fprintf(stderr,"1  %lg\n", aux1);*/
    /*fprintf(stderr,"2  %lg\n", aux2);*/
    /*fprintf(stderr,"3  %lg\n", aux3);*/

    return (aux1 + aux2 + aux3) / 4;
}

static double rate_lepton_gains(double e1, double e2, double g)
{
    double ee = e1 * e2;
    double E  = e1 + e2;

    /*double cp = pow(e1 - g, 2) - 1;*/
    /*double cm = pow(e2 - g, 2) - 1;*/
    double dp = e1 * e1 + ee + g * (e2 - e1);
    double dm = e2 * e2 + ee - g * (e2 - e1);

    double aux = g * (E - g) + 1;
    /*double e_min_aux2 = (aux - sqrt(aux * aux - E * E)) / 2;*/
    /*double e_max_aux2 = (aux + sqrt(aux * aux - E * E)) / 2;*/
    double e_min_aux2 = aux * (1 - sqrt(1 - E * E / aux / aux)) / 2;
    double e_max_aux2 = aux * (1 + sqrt(1 - E * E / aux / aux)) / 2;

    if(E / aux < 2e-5)
    {
        double EE   = E/aux;
        double eeee = EE*EE;

        double t = E*EE/4;
        double tt = t * (1 + eeee/4);

        e_min_aux2 = tt;
        e_max_aux2 = aux * (2 - E * E / aux / aux/2) / 2;

        if(e1 == g || e2 == g)
            e_max_aux2 = ee + (4 - (e1 - e2)*(e1-e2)) / (4*(ee+1));
    }

    if(ee == e_max_aux2 && ((e1 == g) || (e2 == g)))
        e_max_aux2 = e_max_aux2 * (1 - 1e-15);

    double e_min2 = fmax(1 , e_min_aux2);
    double e_max2 = fmin(ee, e_max_aux2);
    double e_min  = sqrt(e_min2);
    double e_max  = sqrt(e_max2);

    /*fprintf(stderr,"%lg %lg\n", ee, e_max_aux2);*/
    /*fprintf(stderr,"%lg %lg\n", e_min, e_max);*/

    if(e_min > e_max) return 0;

    double aux1 = sqrt(E*E - 4 * e_max2) / 4;
    double aux4 = sqrt(E*E - 4 * e_min2) / 4;

    double aux2 = fabs(e1 - g) != 1 ? Jpm(e_min2, e_max2, ee, e1 - g, dp) : H_0(e_max, ee, dp) - H_0(e_min, ee, dp);
    double aux3 = fabs(e2 - g) != 1 ? Jpm(e_min2, e_max2, ee, e2 - g, dm) : H_0(e_max, ee, dm) - H_0(e_min, ee, dm);

    /* This is a math result, not an approximation */
    if(ee < e_max_aux2)
        aux1 = fabs(e1 - e2) / 4;

    if(E / aux < 5e-6)
    {
        aux4 = E * sqrt(1 - 1 / aux) / 4;
        /*aux4 = E * sqrt(1 - (1+E*E/aux/aux/4) / aux) / 4;*/

        if(e_max2 == ee)
        {
            aux1 = (-2*fmin(e1,e2) + E/(2*aux)) / 4;
            aux4 = 0;
        }
    }

    /*fprintf(stderr,"cp cm %lg %lg\n", cp, cm);*/
    /*fprintf(stderr,"dp dm %lg %lg\n", dp, dm);*/

    /*fprintf(stderr,"%lg\n", aux1);*/
    /*fprintf(stderr,"%lg\n", aux2);*/
    /*fprintf(stderr,"%lg\n", aux3);*/
    /*fprintf(stderr,"%lg\n", aux4);*/

    /*fprintf(stderr,"1 - 4 %lg\n", aux1 - aux4);*/
    /*fprintf(stderr,"2 %lg\n", aux2);*/
    /*fprintf(stderr,"3 %lg\n", aux3);*/

    return (aux1 - aux4 + aux2 + aux3) / ee;
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

            gains_inner[j] = 0;

            if(index_e_min >= index_e_max) continue;

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
            gains += 2 * st->photons.population[j] * gains_inner[j];

        gains *= dlnx / 2;

        st->pair_production_lepton_gains[i] = factor * gains;
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

            for(k = index_e_min; k < st->photons.size; k++)
            {
                unsigned int index = index_base2 + k;

                double e2 = st->photons.energy[k];

                st->pair_production_LUT_lepton_gains_R[index] = rate_lepton_gains(e1, e2, g);
            }
        }
    }
}
