#include "inverse_compton.h"
#include "constants.h"

#include <gsl/gsl_sf_dilog.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <float.h>

#define MAX(x,y) ((x) < (y) ? (y) : (x))

double g1(double g, double e, double es)
{
    double aux1 = 2 * g * g;
    double q = 2 * aux1 * es / e;

    return ((q - 1) * (1 + 2 / q) - 2 * log(q)) / aux1;
}

double g2(double g, double e, double es)
{
    double aux1 = es / (g - es);
    double q = aux1 / (4 * e * g);

    return 2 * (2 * q * log(q) + (1 + 2 * q) * (1 - q) + aux1 * es / g * (1 - q) / 2);
}

void inverse_compton_process_head_on(state_t *st)
{
    double factor = LIGHT_SPEED * ELECTRON_RADIUS * ELECTRON_RADIUS * M_PI;
    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    unsigned int i, j, k;

    for(i = 0; i < st->photons.size; i++)
    {
        unsigned int index_base1 = i * st->electrons.size;

        // double es = st->photons.energy[i];

        // Photon downscattering
        double downscattering = 0;

        for(j = 0; j < st->electrons.size; j++)
        {
            // double g = st->electrons.energy[j];
            // double e_min = es;
            // double e_max = 4 * g * g * es;

            unsigned int index_e_min = st->inverse_compton_LUT_downscattering_e_min_index[index_base1 + j];
            unsigned int index_e_max = st->inverse_compton_LUT_downscattering_e_max_index[index_base1 + j];

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            st->inverse_compton_inner_downscattering[j] = 0;

            st->inverse_compton_inner_downscattering[j] += st->photons.population[index_e_min] * st->inverse_compton_LUT_g1[index_base2 + index_e_min];
            st->inverse_compton_inner_downscattering[j] += st->photons.population[index_e_max] * st->inverse_compton_LUT_g1[index_base2 + index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double n = st->photons.population[k];

                unsigned int index = index_base2 + k;

                /*st->inverse_compton_inner_downscattering[j] += 2 * n * g1(g, e, es);*/
                st->inverse_compton_inner_downscattering[j] += 2 * n * st->inverse_compton_LUT_g1[index];
            }

            st->inverse_compton_inner_downscattering[j] *= dlne / 2;
        }

        downscattering += st->electrons.population[0] * st->inverse_compton_inner_downscattering[0];
        downscattering += st->electrons.population[st->electrons.size - 1] * st->inverse_compton_inner_downscattering[st->electrons.size - 1];

        for(j = 1; j < st->electrons.size - 1; j++)
        {
            downscattering += 2 * st->electrons.population[j] * st->inverse_compton_inner_downscattering[j];
        }

        st->inverse_compton_photon_gains[i] = factor * downscattering * dlng / 2;


        // Upscattering
        double upscattering = 0;

        /*double g_min = MAX(1, es);*/
        // double g_max = INFINITY;

        unsigned int index_g_min = st->inverse_compton_LUT_upscattering_g_min_index[i];
        unsigned int index_g_max = st->electrons.size - 1;

        for(j = index_g_min; j < index_g_max + 1; j++)
        {
            // double g = st->electrons.energy[j];
            // double e_min = es / (4 * g * (g - es));
            // double e_max = es;

            unsigned int index_e_min = st->inverse_compton_LUT_upscattering_e_min_index[index_base1 + j];
            unsigned int index_e_max = st->inverse_compton_LUT_upscattering_e_max_index[index_base1 + j]; // TODO

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            if(index_e_min >= index_e_max) continue;

            st->inverse_compton_inner_upscattering[j] = 0;

            st->inverse_compton_inner_upscattering[j] += st->photons.population[index_e_min] * st->inverse_compton_LUT_g2[index_base2 + index_e_min];
            st->inverse_compton_inner_upscattering[j] += st->photons.population[index_e_max] * st->inverse_compton_LUT_g2[index_base2 + index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double n = st->photons.population[k];

                unsigned int index = index_base2 + k;

                /*st->inverse_compton_inner_upscattering[j] += 2 * n * g2(g, e, es);*/
                st->inverse_compton_inner_upscattering[j] += 2 * n * st->inverse_compton_LUT_g2[index];
            }

            st->inverse_compton_inner_upscattering[j] *= dlne / 2;

            /*fprintf(stderr,"AAA%u:\t%lg\n", j, st->inverse_compton_inner_upscattering[j]);*/
        }

        upscattering += st->electrons.population[index_g_min] * st->inverse_compton_inner_upscattering[index_g_min];
        upscattering += st->electrons.population[index_g_max] * st->inverse_compton_inner_upscattering[index_g_max];

        for(j = index_g_min + 1; j < index_g_max; j++)
        {
            upscattering += 2 * st->electrons.population[j] * st->inverse_compton_inner_upscattering[j];
        }

        st->inverse_compton_photon_gains[i] += factor * upscattering * dlng / 2;

        /*fprintf(stderr,"%u:\t%lg\t%lg\n", i, downscattering, upscattering);*/
    }
}

void inverse_compton_process_head_on_downscattering(state_t *st)
{
    double factor = LIGHT_SPEED * ELECTRON_RADIUS * ELECTRON_RADIUS * M_PI;
    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    unsigned int i, j, k;

    for(i = 0; i < st->photons.size; i++)
    {
        unsigned int index_base1 = i * st->electrons.size;

        // double es = st->photons.energy[i];

        for(j = 0; j < st->electrons.size; j++)
        {
            // double g = st->electrons.energy[j];
            // double e_min = es;
            // double e_max = 4 * g * g * es;

            unsigned int index_e_min = st->inverse_compton_LUT_downscattering_e_min_index[index_base1 + j];
            unsigned int index_e_max = st->inverse_compton_LUT_downscattering_e_max_index[index_base1 + j]; // TODO

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            st->inverse_compton_inner_downscattering[j] = 0;
            st->inverse_compton_inner_downscattering[j] += st->photons.population[index_e_min] * st->inverse_compton_LUT_g1[index_base2 + index_e_min];
            st->inverse_compton_inner_downscattering[j] += st->photons.population[index_e_max] * st->inverse_compton_LUT_g1[index_base2 + index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double n = st->photons.population[k];
                // double e = st->photons.energy[k];

                unsigned int index = index_base2 + k;

                /*st->inverse_compton_inner_downscattering[j] += 2 * n * g1(g, e, es);*/
                st->inverse_compton_inner_downscattering[j] += 2 * n * st->inverse_compton_LUT_g1[index];
            }

            st->inverse_compton_inner_downscattering[j] *= dlne / 2;
        }

        double downscattering = 0;
        downscattering += st->electrons.population[0] * st->inverse_compton_inner_downscattering[0];
        downscattering += st->electrons.population[st->electrons.size - 1] * st->inverse_compton_inner_downscattering[st->electrons.size - 1];

        for(j = 1; j < st->electrons.size - 1; j++)
        {
            downscattering += 2 * st->electrons.population[j] * st->inverse_compton_inner_downscattering[j];
        }

        st->inverse_compton_photon_gains_downscattering[i] = factor * downscattering * dlng / 2;
    }
}

void inverse_compton_process_head_on_upscattering(state_t *st)
{
    double factor = LIGHT_SPEED * ELECTRON_RADIUS * ELECTRON_RADIUS * M_PI;
    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    unsigned int i, j, k;

    for(i = 0; i < st->photons.size; i++)
    {
        unsigned int index_base1 = i * st->electrons.size;

        // double g_min = MAX(1, es);
        // double g_max = INFINITY;

        unsigned int index_g_min = st->inverse_compton_LUT_upscattering_g_min_index[i];
        unsigned int index_g_max = st->electrons.size - 1;

        for(j = index_g_min; j < index_g_max + 1; j++)
        {
            // double g = st->electrons.energy[j];
            // double e_min = es / (4 * g * (g - es));
            // double e_max = es;

            unsigned int index_e_min = st->inverse_compton_LUT_upscattering_e_min_index[index_base1 + j];
            unsigned int index_e_max = st->inverse_compton_LUT_upscattering_e_max_index[index_base1 + j];

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            if(index_e_min >= index_e_max) continue;

            st->inverse_compton_inner_upscattering[j] = 0;
            st->inverse_compton_inner_upscattering[j] += st->photons.population[index_e_min] * st->inverse_compton_LUT_g2[index_base2 + index_e_min];
            st->inverse_compton_inner_upscattering[j] += st->photons.population[index_e_max] * st->inverse_compton_LUT_g2[index_base2 + index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double n = st->photons.population[k];
                /*double e = st->photons.energy[k];*/

                unsigned int index = index_base2 + k;

                /*st->inverse_compton_inner_upscattering[j] += 2 * n * g2(g, e, es);*/
                st->inverse_compton_inner_upscattering[j] += 2 * n * st->inverse_compton_LUT_g2[index];
            }

            st->inverse_compton_inner_upscattering[j] *= dlne / 2;

            /*fprintf(stderr,"AAA%u:\t%lg\n", j, st->inverse_compton_inner_upscattering[j]);*/
        }

        double upscattering = 0;
        upscattering += st->electrons.population[index_g_min] * st->inverse_compton_inner_upscattering[index_g_min];
        upscattering += st->electrons.population[index_g_max] * st->inverse_compton_inner_upscattering[index_g_max];

        for(j = index_g_min + 1; j < index_g_max; j++)
        {
            upscattering += 2 * st->electrons.population[j] * st->inverse_compton_inner_upscattering[j];
        }

        st->inverse_compton_photon_gains_upscattering[i] = factor * upscattering * dlng / 2;

        /*fprintf(stderr,"%u:\t%lg\t%lg\n", i, downscattering, upscattering);*/
    }
}

void inverse_compton_process_photon_losses(state_t *st)
{
    unsigned int i,j;

    double factor = 3 * THOMSON_CROSS_SECTION / 16 * LIGHT_SPEED;
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    for(i = 0; i < st->photons.size; i++)
    {
        double n = st->photons.population[i];

        unsigned int index_base = i * st->electrons.size;

        double losses = 0;
        for(j = 0; j < st->electrons.size; j++)
        {
            double reaction_rate = st->inverse_compton_LUT_losses_reaction_rate[index_base + j];
            double ne = st->electrons.population[j];

            losses += 2 * ne * reaction_rate;
        }

        losses *= dlng;

        st->inverse_compton_photon_losses[i] = -factor * n * losses;
    }
}

void inverse_compton_process_electron_losses(state_t *st)
{
    unsigned int i;
    double factor;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];
    double dlng = st->electrons.log_energy[1] - st->electrons.log_energy[0];

    double photon_energy = 0;
    for(i = 0; i < st->photons.size; i++)
    {
        double e = st->photons.energy[i];
        double n = st->photons.population[i];

        photon_energy += 2 * e * e * n;
    }
    photon_energy *= dlne;

    factor = (4 * LIGHT_SPEED * THOMSON_CROSS_SECTION) / 3 * photon_energy;
    st->photon_energy = photon_energy;
    st->inverse_compton_electron_losses_factor = factor;

    double log_n_new = st->electrons.log_population[st->electrons.size - 1];
    double n_new = exp(log_n_new);
    for(i = st->electrons.size - 2; i < st->electrons.size; i--)
    {
        double aux1 = st->dt * st->electrons.energy[i] * factor / dlng;

        log_n_new = (st->electrons.log_population[i] + aux1 * (log_n_new + 2 * dlng)) / (1 + aux1);

        if(log_n_new < log(DBL_MIN))
        {
            log_n_new = log(DBL_MIN);
            n_new = DBL_MIN;
        }
        else
            n_new = exp(log_n_new);

        st->inverse_compton_electron_losses[i] = (n_new - st->electrons.population[i]) / st->dt;
    }
    st->inverse_compton_electron_losses[st->electrons.size - 1] = 0;
}

void init_inverse_compton_LUT_g1_g2(state_t *st)
{
    size_t size = st->photons.size * st->electrons.size * st->photons.size;

    posix_memalign((void **) &st->inverse_compton_LUT_g1, 32, sizeof(double) * size);
    posix_memalign((void **) &st->inverse_compton_LUT_g2, 32, sizeof(double) * size);

    posix_memalign((void **) &st->inverse_compton_LUT_downscattering_e_min_index, 32, sizeof(double) * st->photons.size * st->electrons.size);
    posix_memalign((void **) &st->inverse_compton_LUT_downscattering_e_max_index, 32, sizeof(double) * st->photons.size * st->electrons.size);
    posix_memalign((void **) &st->inverse_compton_LUT_upscattering_e_min_index,   32, sizeof(double) * st->photons.size * st->electrons.size);
    posix_memalign((void **) &st->inverse_compton_LUT_upscattering_e_max_index,   32, sizeof(double) * st->photons.size * st->electrons.size);
    posix_memalign((void **) &st->inverse_compton_LUT_upscattering_g_min_index,   32, sizeof(double) * st->photons.size);
}

void init_inverse_compton_LUT_losses_reaction_rate(state_t *st)
{
    size_t size = st->photons.size * st->electrons.size;

    posix_memalign((void **) &st->inverse_compton_LUT_losses_reaction_rate, 32, sizeof(double) * size);
}

void calculate_inverse_compton_LUT_g1_g2(state_t *st)
{
    unsigned int i, j, k;

    for(i = 0; i < st->photons.size; i++)
    {
        double es = st->photons.energy[i];

        unsigned int index_base1 = i * st->electrons.size;

        double upscattering_g_min = MAX(1, es);
        unsigned int index_g_min = 0;
        for(/* EMPTY */; index_g_min < st->electrons.size && st->electrons.energy[index_g_min] < upscattering_g_min; index_g_min++);

        st->inverse_compton_LUT_upscattering_g_min_index[i] = index_g_min;

        for(j = 0; j < st->electrons.size; j++)
        {
            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            double g = st->electrons.energy[j];

            double downscattering_e_max = 4 * g * g * es;

            unsigned int index_e_min = i;
            unsigned int index_e_max = i;
            for(/* START */; index_e_max < st->photons.size && st->photons.energy[index_e_max] < downscattering_e_max; index_e_max++);
            index_e_max--;

            st->inverse_compton_LUT_downscattering_e_min_index[index_base1 + j] = index_e_min;
            st->inverse_compton_LUT_downscattering_e_max_index[index_base1 + j] = index_e_max;

            double upscattering_e_min = es / (4 * g * (g - es));
            index_e_min = 0;
            index_e_max = i;

            for(/* START */; index_e_min < st->photons.size && st->photons.energy[index_e_min] < upscattering_e_min; index_e_min++);

            st->inverse_compton_LUT_upscattering_e_min_index[index_base1 + j] = index_e_min;
            st->inverse_compton_LUT_upscattering_e_max_index[index_base1 + j] = index_e_max;

            for(k = 0; k < st->photons.size; k++)
            {
                unsigned int index = index_base2 + k;

                double e = st->photons.energy[k];

                st->inverse_compton_LUT_g1[index] = g1(g, e, es) / g;
                st->inverse_compton_LUT_g2[index] = g2(g, e, es) / g;
            }
        }
    }
}

void calculate_inverse_compton_LUT_losses_reaction_rate(state_t *st)
{
    unsigned int i, j;

    for(i = 0; i < st->photons.size; i++)
    {
        double e = st->photons.energy[i];

        for(j = 0; j < st->electrons.size; j++)
        {
            unsigned int index_base = i * st->electrons.size;

            double g = st->electrons.energy[j];
            double b = sqrt(1 - 1/(g*g));

            double x_min = 2 * e * g * (1 + b);
            double x_max = 2 * e * g * (1 - b);

            double a2 = 1 / (4 * (x_max + 1));
            double a3 = (x_max + 9  + 8 / x_max) / 2;
            double a4 = log1p(x_max);
            double a5 = 2 * gsl_sf_dilog(-x_max);

            double b2 = 1 / (4 * (x_min + 1));
            double b3 = (x_min + 9 + 8 / x_min) / 2;
            double b4 = log1p(x_min);
            double b5 = 2 * gsl_sf_dilog(-x_min);

            // Happens when 1 << g
            if(x_max == 0)
            {
                a3 = 1;
                a4 = 4;
            }

            double res;
            res = -b / (g * e);
            res += (b2 - a2) / (e * e * g * g);
            res += (b3 * b4 - a3 * a4) / (e * e * g * g);
            res += (b5 - a5) / (e * e * g * g);

            if(e * g < 1e-6) res = 16./3 * b;

            res *= g/b;

            st->inverse_compton_LUT_losses_reaction_rate[index_base + j] = res;
        }
    }
}
