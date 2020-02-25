#include "distribution.h"

#include <assert.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

void generate_distribution(double *population, double *energy,
        distribution_metadata_t *dm, unsigned int size)
{
    switch(dm->dt)
    {
        case maxwell_juttner:
            generate_maxwell_juttner(population, energy, dm->t, size);
            break;

        case power_law:
            generate_power_law(population, energy, dm->p, size);
            break;

        case broken_power_law:
            generate_broken_power_law(population, energy, dm->gc, dm->p1, dm->p2, size);
            break;

        case power_law_with_exponential_cutoff:
            generate_power_law_with_exponential_cutoff(population, energy, dm->p, dm->e, size);
            break;

        case hybrid:
            generate_hybrid(population, energy, dm->t, size);
            break;

        case connected_power_law:
            generate_connected_power_law(population, energy, dm->gc, dm->p1, dm->p2, size);
            break;

        default:
            break;
    }
}

void generate_maxwell_juttner(double *population, double *energy,
        double theta, unsigned int size)
{
    unsigned int i;
    double norm = 1 / (theta * gsl_sf_bessel_Kn(2, 1 / theta));

    for(i = 0; i < size; i++)
    {
        double aux0 = energy[i] * sqrt(energy[i] * energy[i] - 1);

        population[i] = aux0 * norm * exp(-energy[i] / theta);
    }
}

void generate_power_law(double *population, double *energy,
        double p, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    if(p == 1)
    {
        double norm = log(energy_min / energy_max);

        for(i = 0; i < size; i++)
            population[i] = norm / energy[i];
    }
    else
    {
        double norm = (1 - p) / (pow(energy_max, 1-p) - pow(energy_min, 1-p));
        for(i = 0; i < size; i++)
            population[i] = norm * pow(energy[i], -p);
    }
}

void generate_broken_power_law(double *population, double *energy,
        double gc, double p1, double p2, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm1 = p1 == 1 ?
                    log(energy_min / gc) :
                    (pow(gc, 1-p1) - pow(energy_min, 1-p1)) / (1 - p1);

    double norm2 = p2 == 1 ?
                    log(gc / energy_max) :
                    (pow(energy_max, 1-p2) - pow(gc, 1-p2)) / (1 - p2);

    double aux  = pow(gc, p2 - p1);
    double norm = 1/(norm1 + aux * norm2);

    for(i = 0; i < size && energy[i] < gc; i++)
        population[i] = norm * pow(energy[i], -p1);
    for(     ; i < size                  ; i++)
        population[i] = aux * norm * pow(energy[i], -p2);
}

void generate_power_law_with_exponential_cutoff(double *population, double *energy,
        double p, double e, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm = pow(e, 1 - p) * gsl_sf_gamma_inc(1 - p, energy_min / e);

    for(i = 0; i < size; i++)
        population[i] = norm * pow(energy[i], -p) * exp(-energy[i] / e);
}

struct hybrid_norm_mj_params { double theta; };
static double hybrid_norm_mj(double x, void *params)
{
    struct hybrid_norm_mj_params p = *(struct hybrid_norm_mj_params *)params;

    return x * sqrt(x*x - 1) * exp(-x/p.theta);
}

void generate_hybrid(double *population, double *energy,
        double theta, unsigned int size)
{
    unsigned int i;
    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);

    gsl_function F;
    F.function = &hybrid_norm_mj;
    F.params = &theta;

    double g_c = theta < 0.73262 ?
                    1 + 1.2 * theta + 1.9 * theta*theta
                 :
                    4 * theta - 1/(8*theta) + 7/(128*theta*theta*theta);

    double p = (g_c*g_c*g_c - 2 * theta * g_c*g_c - g_c + theta) /
               (theta * (g_c*g_c - 1));

    double norm_mj;
    gsl_integration_qags(&F, 1, g_c, 0, 1e-7, 256, w, &norm_mj, &e);
    norm_mj = 1 / norm_mj;

    double norm_pl = (1 - p) / (pow(energy[size - 1], 1-p) - pow(g_c, 1-p));

    double a_mj = 1 / ((1 - pow(g_c / energy[size - 1], p-1)) * g_c * g_c * sqrt(g_c * g_c - 1) * norm_mj/((p-1)*exp(g_c/theta)) + 1);
    double a_pl = 1 - a_mj;

    for(i = 0; i < size && energy[i] < g_c; i++)
    {
        double g = energy[i];
        population[i] = a_mj * g * sqrt(g*g - 1) * exp(-g / theta) * norm_mj;
    }

    for(     ; i < size; i++)
    {
        double g = energy[i];
        population[i] = a_pl * norm_pl * pow(g, -p);
    }

    gsl_integration_workspace_free(w);
}

struct cpl_norm_params { double p1, p2; };
static double cpl_norm_lo(double x, void *params)
{
    struct cpl_norm_params p = *(struct cpl_norm_params *)params;

    return pow(x, p.p1) / (1 + pow(x, p.p1 + p.p2));
}
static double cpl_norm_hi(double x, void *params)
{
    struct cpl_norm_params p = *(struct cpl_norm_params *)params;

    return pow(x, p.p2 - 2) / (1 + pow(x, p.p1 + p.p2));
}

void generate_connected_power_law(double *population, double *energy,
        double g_c, double p1, double p2, unsigned int size)
{
    unsigned int i;

    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);
    struct cpl_norm_params params = {p1, p2};

    gsl_function F;
    F.function = &cpl_norm_lo;
    F.params = &params;

    double norm_lo, norm_hi;

    gsl_integration_qags(&F, energy[0] / g_c,        1, 0, 1e-7, 256, w, &norm_lo, &e);
    F.function = &cpl_norm_hi;
    gsl_integration_qags(&F, g_c / energy[size - 1], 1, 0, 1e-7, 256, w, &norm_hi, &e);

    double norm = g_c * (norm_lo + norm_hi);

    for(i = 0; i < size; i++)
    {
        double g = energy[i] / g_c;
        population[i] = pow(g, p1) / (norm * (1 + pow(g, p1 + p2)));
    }

    gsl_integration_workspace_free(w);
}


double distribution_average(distribution_metadata_t *dm)
{
    switch(dm->dt)
    {
        case maxwell_juttner:
            return maxwell_juttner_average(dm->t);

        case power_law:
            return power_law_average(dm->min, dm->max, dm->p);

        case broken_power_law:
            return broken_power_law_average(dm->min, dm->max, dm->gc, dm->p1, dm->p2);

        case power_law_with_exponential_cutoff:
            assert(0);
            break;

        case hybrid:
            assert(0);
            break;

        case connected_power_law:
            assert(0);
            break;

        default:
            assert(0);
            break;
    }
}


struct mj_average_params { double theta; };
static double mj_average(double x, void *params)
{
    struct hybrid_norm_mj_params p = *(struct hybrid_norm_mj_params *)params;

    return x*x * sqrt(x*x - 1) * exp(-x/p.theta);
}

double maxwell_juttner_average(double theta)
{
    double norm = 1 / (theta * gsl_sf_bessel_Kn(2, 1 / theta));

    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);

    gsl_function F;
    F.function = &mj_average;
    F.params = &theta;

    double average;
    gsl_integration_qags(&F, 1, 1e10, 0, 1e-7, 256, w, &average, &e);

    return norm * average;
}

double power_law_average(double energy_min, double energy_max, double p)
{
    double norm = p == 1 ?
                    log(energy_min / energy_max) :
                    (1 - p) / (pow(energy_max, 1-p) - pow(energy_min, 1-p));

    double average = p == 2 ?
                        log(energy_max / energy_min) :
                        (pow(energy_max, 2-p) - pow(energy_min, 2-p)) / (2 - p);

    return average * norm;
}

double broken_power_law_average(double energy_min, double energy_max,
        double energy_break, double p1, double p2)
{
    double norm1 = p1 == 1 ?
                    log(energy_min / energy_break) :
                    (pow(energy_break, 1-p1) - pow(energy_min, 1-p1)) / (1 - p1);

    double norm2 = p2 == 1 ?
                    log(energy_break / energy_max) :
                    (pow(energy_max, 1-p2) - pow(energy_break, 1-p2)) / (1 - p2);

    double avg1 = p1 == 2 ?
                    log(energy_break / energy_min) :
                    (pow(energy_break, 2-p1) - pow(energy_min, 2-p1)) / (2 - p1);

    double avg2 = p2 == 2 ?
                    log(energy_max / energy_break) :
                    (pow(energy_max, 2-p2) - pow(energy_break, 2-p2)) / (2 - p2);

    double aux  = pow(energy_break, p2 - p1);
    double norm = 1/(norm1 + aux * norm2);

    return norm * (avg1 + aux * avg2);
}
