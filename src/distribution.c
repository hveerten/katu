#include "distribution.h"

#include <assert.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

/* >>> Maxwell-Jüttner Functions */
static inline double maxwell_juttner_norm(double theta)
{
    return 1 / (theta * gsl_sf_bessel_Kn(2, 1 / theta));
}

static inline double maxwell_juttner_average(double theta)
{
    return gsl_sf_bessel_Kn(3, 1 / theta) / gsl_sf_bessel_Kn(2, 1 / theta) - theta;
}

static inline void generate_maxwell_juttner(double *population, double *energy,
        double theta, unsigned int size)
{
    unsigned int i;
    double norm = maxwell_juttner_norm(theta);

    for(i = 0; i < size; i++)
    {
        double aux0 = energy[i] * sqrt(energy[i] * energy[i] - 1);

        population[i] = aux0 * norm * exp(-energy[i] / theta);
    }
}
/* <<< Maxwell-Jüttner Functions */

/* >>> Hybrid Functions */
/* The hybrid function is defined as
 * f(x, t, gc) =
 *      N *     x * sqrt(x^2 - 1) * exp(-x / t)     for x < gc
 *      N * C * x^-p                                for x > gc
 * with
 *  p = (gc^3 - 2 * theta * gc^2 - gc + theta) / (theta * (gc*gc - 1));
 *  C = gc^(1 + p) * sqrt(gc^2 - 1) * exp(-gc / theta)
 *
 *  NOTE: gc should be bigger than the maximum of the distribution so that p
 *          is negative.
 *  So, gc should be bigger than gM with
 *  gM =
 *      1 + t/2 + 5/8*t^2 + t^3/2               for t < 0.5
 *      2*t + 1 / (4*t) - 1/ (256*t**3)         for t > 0.5
 *
 *  NOTE: It is the user's responsibility to make sure that t and gc are correct
 */

struct hybrid_params { double theta, gc, p, C;};
static double h_norm(double x, void *params)
{
    struct hybrid_params p = *(struct hybrid_params *)params;

    if (x < p.gc)
        return x * sqrt(x*x - 1) * exp(-x/p.theta);
    else
        return p.C * pow(x, -p.p);
}

static double h_avg(double x, void *params)
{
    struct hybrid_params p = *(struct hybrid_params *)params;

    if (x < p.gc)
        return x*x * sqrt(x*x - 1) * exp(-x/p.theta);
    else
        return p.C * pow(x, 1-p.p);
}

static inline double hybrid_norm(double min, double max, double theta, double gc)
{
    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);

    double p = (gc*gc*gc - 2 * theta * gc*gc - gc + theta) /
               (theta * (gc*gc - 1));

    double C = pow(gc, 1 + p) * sqrt(gc*gc - 1) * exp(-gc / theta);

    struct hybrid_params params = {theta, gc, p, C};

    gsl_function F;
    F.function = &h_norm;
    F.params = &params;

    double n1, n2;

    gsl_integration_qag(&F, min,  gc, 0, 1e-7, 256, GSL_INTEG_GAUSS61, w, &n1, &e);
    n2 = p == 1 ?
            log(max / gc) : (pow(max, 1 - p) - pow(gc, 1 - p)) / (1 - p);

    gsl_integration_workspace_free(w);

    return 1 / (n1 + C * n2);
}

static inline double hybrid_average(double min, double max, double theta, double gc)
{
    double e;
    double norm;
    double avg1, avg2;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);

    double p = (gc*gc*gc - 2 * theta * gc*gc - gc + theta) /
               (theta * (gc*gc - 1));

    double C = pow(gc, 1 + p) * sqrt(gc*gc - 1) * exp(-gc / theta);

    struct hybrid_params params = {theta, gc, p, C};

    norm = hybrid_norm(min, max, theta, gc);

    gsl_function F;
    F.params = &params;

    F.function = &h_avg;
    gsl_integration_qag(&F, min, max, 0, 1e-7, 256, GSL_INTEG_GAUSS61, w, &avg1, &e);
    avg2 = p == 2 ?
        log(max / gc) : (pow(max, 2 - p) - pow(gc, 2 - p)) / (2 - p);

    gsl_integration_workspace_free(w);

    return norm * (avg1 + C * avg2);
}


void generate_hybrid(double *population, double *energy,
        double theta, double gc, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double p = (gc*gc*gc - 2 * theta * gc*gc - gc + theta) /
               (theta * (gc*gc - 1));

    double C = pow(gc, 1 + p) * sqrt(gc*gc - 1) * exp(-gc / theta);

    double norm = hybrid_norm(energy_min, energy_max, theta, gc);

    for(i = 0; i < size && energy[i] < gc; i++)
    {
        double g = energy[i];
        population[i] = norm * g * sqrt(g*g - 1) * exp(-g / theta);
    }

    for(     ; i < size; i++)
    {
        double g = energy[i];
        population[i] = norm * C * pow(g, -p);
    }
}
/* <<< Hybrid Functions */

/* >>> Power-Law Functions */
static inline double power_law_norm(double min, double max, double p)
{
    return p == 1 ?
        log(min / max) : (1 - p) / (pow(max, 1 - p) - pow(min, 1 - p));
}

static inline double power_law_average(double min, double max, double p)
{
    double norm = power_law_norm(min, max, p);

    double avg = p == 2 ?
        log(max / min) : (pow(max, 2 - p) - pow(min, 2 - p)) / (2 - p);

    return norm * avg;
}

static inline void generate_power_law(double *population, double *energy,
        double p, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm = power_law_norm(energy_min, energy_max, p);

    for(i = 0; i < size; i++)
        population[i] = norm * pow(energy[i], -p);
}
/* <<< Power-Law Functions */

/* >>> Power-Law with Exponential Cutoff Functions */
/* For now, we will use the simpler formulas of the power law to calculate the
 * normalization constant and the average
 * TODO: Improve this (if neccessary)
 */
static inline double power_law_with_exponential_cutoff_norm(
        double min, double max, double p, double e)
{
    (void) e;
    return power_law_norm(min, max, p);
}

static inline double power_law_with_exponential_cutoff_average(
        double min, double max, double p, double e)
{
    (void) e;
    return power_law_average(min, max, p);
}

void generate_power_law_with_exponential_cutoff(double *population, double *energy,
        double p, double e, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm = power_law_with_exponential_cutoff_norm(energy_min, energy_max, p, e);

    for(i = 0; i < size; i++)
        population[i] = norm * pow(energy[i], -p) * exp(-energy[i] / e);
}
/* <<< Power-Law with Exponential Cutoff Functions */

/* >>> Broken Power Law Functions */
static inline double broken_power_law_norm(double min, double max,
        double gc, double p1, double p2)
{
    double n1 = 1 / power_law_norm(min,  gc, p1);
    double n2 = 1 / power_law_norm(gc,  max, p2);

    return 1 / (n1 + pow(gc, p2 - p1) * n2);
}

static inline double broken_power_law_average(double min, double max,
        double gc, double p1, double p2)
{
    double norm = broken_power_law_norm(min, max, gc, p1, p2);

    double a1 = p1 == 2 ?
        log(gc / min) : (pow(gc,  2 - p1) - pow(min, 2 - p1)) / (2 - p1);
    double a2 = p2 == 2 ?
        log(max / gc) : (pow(max, 2 - p2) - pow(gc,  2 - p2)) / (2 - p2);

    return norm * (a1 + pow(gc, p2 - p1) * a2);
}

static void generate_broken_power_law(double *population, double *energy,
        double gc, double p1, double p2, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm = broken_power_law_norm(energy_min, energy_max, gc, p1, p2);
    double aux  = pow(gc, p2 - p1);

    for(i = 0; i < size && energy[i] < gc; i++)
        population[i] =       norm * pow(energy[i], -p1);
    for(     ; i < size                  ; i++)
        population[i] = aux * norm * pow(energy[i], -p2);
}
/* <<< Broken Power Law Functions */

/* >>> Connected Power Law Functions */
struct cpl_params { double p1, p2; };
static double cpl_norm_lo(double x, void *params)
{
    struct cpl_params p = *(struct cpl_params *)params;

    return pow(x, p.p1) / (1 + pow(x, p.p1 + p.p2));
}
static double cpl_norm_hi(double x, void *params)
{
    struct cpl_params p = *(struct cpl_params *)params;

    return pow(x, p.p2 - 2) / (1 + pow(x, p.p1 + p.p2));
}
static inline double connected_power_law_norm(double min, double max,
        double gc, double p1, double p2)
{
    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);
    struct cpl_params params = {p1, p2};

    gsl_function F;
    F.params = &params;

    double norm_lo, norm_hi;

    F.function = &cpl_norm_lo;
    gsl_integration_qags(&F, min / gc, 1, 0, 1e-7, 256, w, &norm_lo, &e);
    F.function = &cpl_norm_hi;
    gsl_integration_qags(&F, gc / max, 1, 0, 1e-7, 256, w, &norm_hi, &e);

    gsl_integration_workspace_free(w);

    return 1 / gc * (norm_lo + norm_hi);
}

static double cpl_avg_lo(double x, void *params)
{
    struct cpl_params p = *(struct cpl_params *)params;

    return pow(x, p.p1 + 1) / (1 + pow(x, p.p1 + p.p2));
}
static double cpl_avg_hi(double x, void *params)
{
    struct cpl_params p = *(struct cpl_params *)params;

    return pow(x, p.p2 - 1) / (1 + pow(x, p.p1 + p.p2));
}
static inline double connected_power_law_average(double min, double max,
        double gc, double p1, double p2)
{
    double e;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);
    struct cpl_params params = {p1, p2};

    gsl_function F;
    F.params = &params;

    double avg_lo, avg_hi;

    F.function = &cpl_avg_lo;
    gsl_integration_qags(&F, min / gc, 1, 0, 1e-7, 256, w, &avg_lo, &e);
    F.function = &cpl_avg_hi;
    gsl_integration_qags(&F, gc / max, 1, 0, 1e-7, 256, w, &avg_hi, &e);

    gsl_integration_workspace_free(w);

    double norm = connected_power_law_norm(min, max, gc, p1, p2);

    return norm * gc * gc * (avg_lo + avg_hi);
}

static void generate_connected_power_law(double *population, double *energy,
        double gc, double p1, double p2, unsigned int size)
{
    unsigned int i;
    double energy_min = energy[0];
    double energy_max = energy[size - 1];

    double norm = connected_power_law_norm(energy_min, energy_max, gc, p1, p2);

    for(i = 0; i < size; i++)
    {
        double g = energy[i] / gc;
        population[i] = pow(g, p1) / (norm * (1 + pow(g, p1 + p2)));
    }
}
/* <<< Connected Power Law Functions */

/* >>> Black Body Functions */
static inline double black_body_norm(double theta)
{
    return 1 / (theta * theta * 2.40411380631918); // zeta(3) * Gamma(3)
}

static inline double black_body_average(double theta)
{
    return theta * pow(M_PI, 4) / 90 / 2.40411380631918; // zeta(3) * Gamma(3)
}

void generate_black_body(double *population, double *energy,
        double theta, unsigned int size)
{
    unsigned int i;
    double norm = black_body_norm(theta);

    for(i = 0; i < size; i++)
    {
        double aux0 = energy[i] * energy[i];

        population[i] = aux0 * norm / expm1(energy[i] / theta);
    }
}
/* <<< Black Body Functions */


void generate_distribution(double *population, double *energy,
        distribution_metadata_t *dm, unsigned int size)
{
    switch(dm->dt)
    {
        case maxwell_juttner:
            generate_maxwell_juttner(population, energy, dm->t, size);
            break;

        case hybrid:
            generate_hybrid(population, energy, dm->t, dm->gc, size);
            break;

        case power_law:
            generate_power_law(population, energy, dm->p, size);
            break;

        case power_law_with_exponential_cutoff:
            generate_power_law_with_exponential_cutoff(population, energy, dm->p, dm->e, size);
            break;

        case broken_power_law:
            generate_broken_power_law(population, energy, dm->gc, dm->p1, dm->p2, size);
            break;

        case connected_power_law:
            generate_connected_power_law(population, energy, dm->gc, dm->p1, dm->p2, size);
            break;

        case black_body:
            generate_black_body(population, energy, dm->t, size);
            break;


        default:
            break;
    }
}

/* Returns the normalization factor the corresponding distributions assuming
 * that only the dependencies with energy and so are left.
 *
 * NOTE: For the case of maxwell-jüttner, black_body and hybrid distributions
 *       the norm is calculated as integrated from 1 to infinity and from 0 to
 *       infinity, respectively. In the other cases from dm->min to dm->max
 */
double distribution_norm(distribution_metadata_t *dm)
{
    switch(dm->dt)
    {
        case maxwell_juttner:
            return maxwell_juttner_norm(dm->t);

        case hybrid:
            return hybrid_norm(dm->min, dm->max, dm->t, dm->gc);

        case power_law:
            return power_law_norm(dm->min, dm->max, dm->p);

        case power_law_with_exponential_cutoff:
            return power_law_with_exponential_cutoff_norm(dm->min, dm->max, dm->p, dm->e);

        case broken_power_law:
            return broken_power_law_norm(dm->min, dm->max, dm->gc, dm->p1, dm->p2);

        case connected_power_law:
            return connected_power_law_norm(dm->min, dm->max, dm->gc, dm->p1, dm->p2);

        case black_body:
            return black_body_norm(dm->t);

        default:
            assert(0 && "Non supported distribution");
            break;
    }
}

double distribution_average(distribution_metadata_t *dm)
{
    switch(dm->dt)
    {
        case maxwell_juttner:
            return maxwell_juttner_average(dm->t);

        case hybrid:
            return hybrid_average(dm->min, dm->max, dm->t, dm->gc);

        case power_law:
            return power_law_average(dm->min, dm->max, dm->p);

        case power_law_with_exponential_cutoff:
            return power_law_with_exponential_cutoff_average(dm->min, dm->max, dm->p, dm->e);

        case broken_power_law:
            return broken_power_law_average(dm->min, dm->max, dm->gc, dm->p1, dm->p2);

        case connected_power_law:
            return connected_power_law_average(dm->min, dm->max, dm->gc, dm->p1, dm->p2);

        case black_body:
            return black_body_average(dm->t);

        default:
            assert(0 && "Non supported distribution");
            break;
    }
}
