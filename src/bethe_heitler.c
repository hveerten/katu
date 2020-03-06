#include "bethe_heitler.h"
#include "constants.h"

#include "utils.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <stdio.h>

static double rate_bethe_heitler_low(double e)
{
    double eta = (e - 2) / (e + 2);
    double rho = 2 * eta / (1 + sqrt(1 - eta * eta));

    double aux0 = pow((e - 2) / e, 3);

    double aux1 = 0.5 * rho;
    double aux2 = 23/ 40. * pow(rho, 2);
    double aux3 = 11/ 60. * pow(rho, 3);
    double aux4 = 29/960. * pow(rho, 4);

    return aux0 * (1 + aux1 + aux2 + aux3 + aux4);
}

static double rate_bethe_heitler_high(double e)
{
    double log2e = log(2 * e);
    double mu2 = 4 / (e*e);

    double a = 28/9. * log2e - 218/27.;

    double b_aux0 = M_PI * M_PI / 6 - 3.5 + 2 * 1.2020569;
    double b_aux1 = (6 - M_PI * M_PI / 3) * log2e;
    double b_aux2 = pow(log2e, 2);
    double b_aux3 = 2 / 3. * pow(log2e, 3);

    double c = 3 / 16. * log2e + 1 / 8.;

    double d = 29 / 2304. * log2e - 77 / 13824.;

    return a +     mu2     * (b_aux0 + b_aux1 - b_aux2 + b_aux3) - \
               pow(mu2, 2) * c - \
               pow(mu2, 3) * d;
}

// NOTE: Check that the global factor is correct
static double rate_bethe_heitler(double e)
{
    double factor = FINE_STRUCTURE_CONSTANT * (ELECTRON_RADIUS * ELECTRON_RADIUS);
    double factor_low  = factor * (2 * M_PI) / 3;
    double factor_high = factor * 1;

    if(e < 3.4)
        return factor_low  * rate_bethe_heitler_low(e);
    else
        return factor_high * rate_bethe_heitler_high(e);
}

static double inelasticity_bethe_heitler_low(double e)
{
    double aux0 = pow((e - 2) / e, 2);
    double eta  = (e - 2) / (e + 2);

    double aux1 =  0.95 * eta;
    double aux2 = -0.88 * eta * eta;

    return M_PI * aux0 * (1 + aux1 + aux2);
}

static double inelasticity_bethe_heitler_mid(double e)
{
    double e_max = 21.21;
    double inelasticity_max = 2.818;

    double eta  = (e - e_max) / (e + e_max);

    double aux1 = -0.75 * eta * eta;
    double aux2 =  0.46 * eta * eta * eta;
    double aux3 = -0.40 * eta * eta * eta * eta;

    return inelasticity_max * (1 + aux1 + aux2 + aux3);
}

static double inelasticity_bethe_heitler_high(double e)
{
    double loge = log(e);

    double aux0 = pow(loge, 3) / e;

    double aux1 = -1.35 / loge;
    double aux2 =  2.5  / (loge * loge);

    return 2.542 * aux0 * (1 + aux1 + aux2);
    // 2.557    -1.57   3.5
}

// TODO: Find the correct global factor
static double inelasticity_bethe_heitler(double e)
{
    double factor = FINE_STRUCTURE_CONSTANT * ELECTRON_RADIUS * ELECTRON_RADIUS / 2;
    factor *= 1 / PROTON_MASS;
    /*factor *= 1 / ELECTRON_MASS;*/

    if(e < 8.2)
        return factor * inelasticity_bethe_heitler_low(e);
    else if(e < 121)
        return factor * inelasticity_bethe_heitler_mid(e);
    else
        return factor * inelasticity_bethe_heitler_high(e);
}


// TODO: Make sure that the global factor is correct
void bethe_heitler_process_lepton_gains(state_t *st)
{
    unsigned int i, j;

    /*double factor = ELECTRON_MASS * LIGHT_SPEED * LIGHT_SPEED * LIGHT_SPEED;*/
    /*double factor = pow(ELECTRON_MASS, 2) * pow(LIGHT_SPEED, 5);*/
    double factor = LIGHT_SPEED;
    double dlnx = st->photons.log_energy[1] - st->photons.log_energy[0];

    double *proton_interpolated;
    gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_steffen, st->protons.size);
    gsl_interp_accel *proton_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, st->protons.log_energy, st->protons.log_population, st->protons.size);

    posix_memalign((void **) &proton_interpolated, 32, st->photons.size * sizeof(double));

    for(i = 0; i < st->electrons.size; i++)
    {
        unsigned int index_base = i * st->photons.size;
        double ge = st->electrons.energy[i];

        for(j = 0; j < st->photons.size; j++)
        {
            double g_proton = st->bethe_heitler_LUT_proton_gamma[index_base + j];

            if(g_proton < st->protons.energy[0] ||
               g_proton > st->protons.energy[st->protons.size - 1])
                proton_interpolated[j] = DBL_MIN;
            else
                proton_interpolated[j] = exp(gsl_spline_eval(proton_spline, log(g_proton), proton_accelerator));
        }

        unsigned int index_photon_min = 0;
        unsigned int index_photon_max = st->photons.size - 1;

        double gains = 0;
        /*gains += proton_interpolated[index_photon_min] * st->photons.population[index_photon_min] * st->photons.energy[index_photon_min] **/
                    /*st->bethe_heitler_LUT_reaction_rate[index_base + index_photon_min] / st->bethe_heitler_LUT_inelasticity[index_base + index_photon_min];*/
        /*gains += proton_interpolated[index_photon_max] * st->photons.population[index_photon_max] * st->photons.energy[index_photon_max] **/
                    /*st->bethe_heitler_LUT_reaction_rate[index_base + index_photon_max] / st->bethe_heitler_LUT_inelasticity[index_base + index_photon_max];*/

        gains += proton_interpolated[index_photon_min] * st->photons.population[index_photon_min] * st->photons.energy[index_photon_min] *
                    st->bethe_heitler_LUT_reaction_rate[index_base + index_photon_min] * st->bethe_heitler_LUT_proton_gamma[index_base + index_photon_min];
        gains += proton_interpolated[index_photon_max] * st->photons.population[index_photon_max] * st->photons.energy[index_photon_max] *
                    st->bethe_heitler_LUT_reaction_rate[index_base + index_photon_max] * st->bethe_heitler_LUT_proton_gamma[index_base + index_photon_max];

        for(j = 1; j < st->photons.size - 1; j++)
        {
            double e        = st->photons.energy[j];
            double n        = st->photons.population[j];
            double np       = proton_interpolated[j];
            double g_proton = st->bethe_heitler_LUT_proton_gamma[index_base + j];
            /*double xi       = st->bethe_heitler_LUT_inelasticity[index_base + j];*/
            double sigma    = st->bethe_heitler_LUT_reaction_rate[index_base + j];

            /*gains += 2 * np * n * e * sigma / xi;*/
            gains += 2 * np * n * e * sigma * g_proton;
        }

        gains *= dlnx / 2;

        st->bethe_heitler_lepton_gains[i] = factor * gains / ge;
    }

    free(proton_interpolated);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(proton_accelerator);
}


void init_bethe_heitler_LUT_lepton_gains(state_t *st)
{
    size_t size1 = st->electrons.size * st->photons.size;

    MEM_PREPARE(st->bethe_heitler_LUT_proton_gamma,  size1, double);
    MEM_PREPARE(st->bethe_heitler_LUT_inelasticity,  size1, double);
    MEM_PREPARE(st->bethe_heitler_LUT_reaction_rate, size1, double);
}

static double f(double loggp, double ge, double e)
{
    double gp = exp(loggp);
    /*return gp - ge * ELECTRON_MASS / PROTON_MASS / inelasticity_bethe_heitler(gp * e);*/
    return loggp - log(ge) - log(ELECTRON_MASS / PROTON_MASS) + log(inelasticity_bethe_heitler(gp * e));
}
void calculate_bethe_heitler_LUT_lepton_gains(state_t *st)
{
    unsigned int i, j;

    for(i = 0; i < st->electrons.size; i++)
    {
        double g = st->electrons.energy[i];

        unsigned int index_base = i * st->photons.size;

        for(j = 0; j < st->photons.size; j++)
        {
            double e = st->photons.energy[j];

            double logg_proton_min = fmax(st->protons.log_energy[0], log(2.001 / e));
            double logg_proton_mid = log(21.21 / e);
            double logg_proton_max = st->protons.log_energy[st->protons.size - 1];

            double value_min = f(logg_proton_min, g, e);
            double value_mid = f(logg_proton_mid, g, e);
            double value_max = f(logg_proton_max, g, e);

            st->bethe_heitler_LUT_proton_gamma[index_base + j]  = st->protons.energy[0];
            st->bethe_heitler_LUT_inelasticity[index_base + j]  = INFINITY;
            st->bethe_heitler_LUT_reaction_rate[index_base + j] = 0;

            /* Here, mid is the point where the function xi has a maximum,
             * so it is possible that f also has a maximum or minimum.
             * This gives the possibility that there are two (?) solutions
             * but as min and max would have the same sign we would find
             * neither.
             * Get at least one of them
             */
            if(logg_proton_min < logg_proton_mid &&
               logg_proton_mid < logg_proton_max)
            {
                if(value_min * value_mid > 0)
                {
                    logg_proton_min = logg_proton_mid;
                    value_min = value_mid;
                }
                else if(value_mid * value_max > 0)
                {
                    logg_proton_max = logg_proton_mid;
                    value_max = value_mid;
                }
                else
                    continue;
            }

            if(logg_proton_min < logg_proton_max && value_min * value_max < 0)
            {
                while(true)
                {
                    logg_proton_mid = (logg_proton_min + logg_proton_max) / 2;
                    value_mid = f(logg_proton_mid, g, e);

                    if(fabs(value_mid) < 1e-4)
                        break;

                    if(value_min * value_mid > 0)
                    {
                        logg_proton_min = logg_proton_mid;
                        value_min = value_mid;
                    }
                    else
                    {
                        logg_proton_max = logg_proton_mid;
                        value_max = value_mid;
                    }
                }

                double g_proton_mid = exp(logg_proton_mid);

                st->bethe_heitler_LUT_proton_gamma[index_base + j]  = g_proton_mid;
                st->bethe_heitler_LUT_inelasticity[index_base + j]  = inelasticity_bethe_heitler(g_proton_mid * e);
                st->bethe_heitler_LUT_reaction_rate[index_base + j] = rate_bethe_heitler(g_proton_mid * e);

                double r  = st->bethe_heitler_LUT_reaction_rate[index_base + j];
                double xi = inelasticity_bethe_heitler(g_proton_mid * e);

                if(!isnormal(r) || r < 0 || !isnormal(xi) || xi <= 0 || xi > 1)
                {
                    fprintf(stderr,"%u %u %lg %lg\n", i, j, r ,xi);
                    break;
                }
            }
        }
    }
}
