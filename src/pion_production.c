#include "constants.h"
#include "state.h"
#include "pion_production.h"

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#include <float.h>
#include <stdlib.h>

const double multi_resonances_constants[][10] =
{
    //                                                                                                   FROM PROTONS            FROM NEUTRONS
    // Low threshold                High threshold              Cross section      energy ratio      M_pi0   M_pi+   M_pi-   M_pi0   M_pi+   M_pi-
    {GEV_TO_EPSILON( 0.5) / 2,      GEV_TO_EPSILON( 0.9) / 2,    60 * 1e-30,       0.10 ,            0.32,   0.34,   0.04,   0.32,   0.04,   0.34}, // Multi resonance 1 Low
    {GEV_TO_EPSILON( 0.5) / 2,      GEV_TO_EPSILON( 0.9) / 2,    60 * 1e-30,       0.40 ,            0.17,   0.29,   0.05,   0.17,   0.05,   0.29}, // Multi resonance 1 High
    {GEV_TO_EPSILON( 0.9) / 2,      GEV_TO_EPSILON( 1.5) / 2,    85 * 1e-30,       0.15 ,            0.42,   0.31,   0.07,   0.42,   0.07,   0.31}, // Multi resonance 2 Low
    {GEV_TO_EPSILON( 0.9) / 2,      GEV_TO_EPSILON( 1.5) / 2,    85 * 1e-30,       0.35 ,            0.19,   0.35,   0.08,   0.19,   0.08,   0.35}, // Multi resonance 2 High
    {GEV_TO_EPSILON( 1.5) / 2,      GEV_TO_EPSILON( 5.0) / 2,   120 * 1e-30,       0.15 ,            0.59,   0.57,   0.30,   0.59,   0.30,   0.57}, // Multi resonance 3 Low
    {GEV_TO_EPSILON( 1.5) / 2,      GEV_TO_EPSILON( 5.0) / 2,   120 * 1e-30,       0.35 ,            0.16,   0.21,   0.13,   0.16,   0.13,   0.21}, // Multi resonance 3 High
    {GEV_TO_EPSILON( 5.0) / 2,      GEV_TO_EPSILON(  50) / 2,   120 * 1e-30,       0.07 ,            1.38,   1.37,   1.11,   1.38,   1.11,   1.37}, // Multi resonance 4 Low
    {GEV_TO_EPSILON( 5.0) / 2,      GEV_TO_EPSILON(  50) / 2,   120 * 1e-30,       0.35 ,            0.16,   0.25,   0.23,   0.16,   0.23,   0.25}, // Multi resonance 4 High
    {GEV_TO_EPSILON(  50) / 2,      GEV_TO_EPSILON( 500) / 2,   120 * 1e-30,       0.02 ,            3.01,   2.86,   2.64,   3.01,   2.64,   2.86}, // Multi resonance 5 Low
    {GEV_TO_EPSILON(  50) / 2,      GEV_TO_EPSILON( 500) / 2,   120 * 1e-30,       0.50 ,            0.20,   0.21,   0.14,   0.20,   0.14,   0.21}, // Multi resonance 5 High
    {GEV_TO_EPSILON( 500) / 2,      GEV_TO_EPSILON(5000) / 2,   120 * 1e-30,       0.007,            5.13,   4.68,   4.57,   5.13,   4.57,   4.68}, // Multi resonance 6 Low
    {GEV_TO_EPSILON( 500) / 2,      GEV_TO_EPSILON(5000) / 2,   120 * 1e-30,       0.50 ,            0.27,   0.29,   0.12,   0.27,   0.12,   0.29}, // Multi resonance 6 High
    {GEV_TO_EPSILON(5000) / 2,      INFINITY,                   120 * 1e-30,       0.002,            7.59,   6.80,   6.65,   7.59,   6.65,   6.80}, // Multi resonance 7 Low
    {GEV_TO_EPSILON(5000) / 2,      INFINITY,                   120 * 1e-30,       0.60 ,            0.26,   0.27,   0.13,   0.26,   0.13,   0.27}, // Multi resonance 7 High

    {GEV_TO_EPSILON( 0.2) / 2,      GEV_TO_EPSILON( 0.5) / 2,   200 * 1e-30,       0.22,             2./3.,  1./3.,  0.00,   2./3.,  0.00,   1./3.}, // Low resonance
    {GEV_TO_EPSILON( 0.5) / 2,      GEV_TO_EPSILON( 1.2) / 2,    90 * 1e-30,       0.22,             0.00,   0.00,   0.34,   0.00,   0.34,   0.00}, // High resonance 1
    {GEV_TO_EPSILON( 0.5) / 2,      GEV_TO_EPSILON( 1.2) / 2,    90 * 1e-30,       0.25,             0.00,   0.77,   0.00,   0.00,   0.00,   0.77}, // High resonance 2
    {GEV_TO_EPSILON( 0.5) / 2,      GEV_TO_EPSILON( 1.2) / 2,    90 * 1e-30,       0.26,             0.47,   0.00,   0.00,   0.47,   0.00,   0.00}, // High resonance 3
};

const double multi_resonances_hadron_losses_constants[][8] =
{
    //                                                                                              FROM PROTONS     FROM NEUTRONS
    // Low threshold                      High threshold        Cross section      Inelasticity     M_p   M_n        M_p   M_n
    {GEV_TO_EPSILON(0.5) / 2,       GEV_TO_EPSILON( 14) / 2,    NAN,               0.60,            0.69, 0.31,      0.31, 0.69}, // Multi resonances

    {GEV_TO_EPSILON(0.2) / 2,       GEV_TO_EPSILON(0.5) / 2,    200 * 1e-30,       0.22,            2./3, 1./3,      1./3, 2./3}, // Low  resonance
    {GEV_TO_EPSILON(0.5) / 2,       GEV_TO_EPSILON(1.2) / 2,     90 * 1e-30,       0.39,            0.57, 0.43,      0.43, 0.57}, // High resonance
};

const double direct_pion_production_constants[][9] =
{
    //                                                                                FROM PROTONS            FROM NEUTRONS
    // Low threshold                        High threshold       energy ratio     M_pi0   M_pi+   M_pi-   M_pi0   M_pi+   M_pi-
    {GEV_TO_EPSILON(0.17) / 2,      GEV_TO_EPSILON(0.56) / 2,    0.13,            0.00,   1.00,   0.00,   0.00,   0.00,   1.00}, // Direct One Pion Low
    {GEV_TO_EPSILON(0.56) / 2,      GEV_TO_EPSILON(  10) / 2,    0.05,            0.00,   1.00,   0.00,   0.00,   0.00,   1.00}, // Direct One Pion Mid
    {GEV_TO_EPSILON(  10) / 2,      INFINITY,                    0.001,           0.00,   1.00,   0.00,   0.00,   0.00,   1.00}, // Direct One Pion High

    {GEV_TO_EPSILON(0.40) / 2,      GEV_TO_EPSILON(1.58) / 2,    0.08,            0.00,   0.25,   0.75,   0.00,   0.75,   0.25}, // Direct Two Pion Low
    {GEV_TO_EPSILON(1.58) / 2,      GEV_TO_EPSILON(  10) / 2,    0.02,            0.00,   0.25,   0.75,   0.00,   0.75,   0.25}, // Direct Two Pion Mid
    {GEV_TO_EPSILON(  10) / 2,      INFINITY,                    0.001,           0.00,   0.25,   0.75,   0.00,   0.75,   0.25}, // Direct Two Pion High

    {GEV_TO_EPSILON(0.40) / 2,      INFINITY,                    0.20,            1./6.,  0.75,   1./12., 1./6.,  1./12., 0.75}, // Direct Two Pion ---
};

const double direct_pion_production_hadron_losses_constants[][7] =
{
    //                                                                            FROM PROTONS      FROM NEUTRONS
    // Low threshold                        High threshold       Inelasticity     M_p    M_n        M_p    M_n
    {GEV_TO_EPSILON(0.17) / 2,      GEV_TO_EPSILON(0.56) / 2,    0.13,            0.00,  1.00,      1.00,  0.00}, // Direct One Pion Low
    {GEV_TO_EPSILON(0.56) / 2,      GEV_TO_EPSILON(  10) / 2,    0.05,            0.00,  1.00,      1.00,  0.00}, // Direct One Pion Mid
    {GEV_TO_EPSILON(  10) / 2,      INFINITY,                    0.001,           0.00,  1.00,      1.00,  0.00}, // Direct One Pion High

    {GEV_TO_EPSILON(0.40) / 2,      GEV_TO_EPSILON(1.58) / 2,    0.28,            5./6., 1./6.,     1./6., 5./6.}, // Direct Two Pion Low
    {GEV_TO_EPSILON(1.58) / 2,      GEV_TO_EPSILON(  10) / 2,    0.22,            5./6., 1./6.,     1./6., 5./6.}, // Direct Two Pion Mid
    {GEV_TO_EPSILON(  10) / 2,      INFINITY,                    0.201,           5./6., 1./6.,     1./6., 5./6.}, // Direct Two Pion High
};


static double f_function(double y, double s, double e_min, double e_max)
{
    if(y < e_min)
        return 0;
    else
    {
        double aux = s / (y*y);
        if(y < e_max)
            return (y * y - e_min * e_min) * aux;
        else
            return (e_max * e_max - e_min * e_min) * aux;
    }
}

static double f_function_multi_resonance(double y, double e_min, double e_max)
{
    if(y < e_min)
        return 0;

    // the x functions are parametrized in
    // log_10(y / GeV)
    double x = log10(EPSILON_TO_GEV(y));

    if(y < e_max)
        return 87.5538 +
               120.894 * x -
               98.4187 * x * x -
               59.6965 * x * x * x +
               67.2251 * x * x * x * x;
    else
        return 131.839   -
               25.3296   * x +
               10.612    * x * x -
               0.858307  * x * x * x +
               0.0493614 * x * x * x * x;

}

static double I_one_pion(double y)
{
    const double e_min = GEV_TO_EPSILON(0.17) / 2;
    const double e_max = GEV_TO_EPSILON(0.96) / 2;

    if(y < e_min)
        return 0;

    // the x functions are parametrized in
    // log_10(y / GeV)
    double x = log10(EPSILON_TO_GEV(y));

    if(y < e_max)
        return 35.9533 +
               84.0859 * x +
               110.765 * x * x +
               102.728 * x * x * x +
               40.4699 * x * x * x * x;
    else
        return 30.2004  +
               40.5478  * x +
               2.03074  * x * x -
               0.387884 * x * x * x +
               0.025044 * x * x * x * x;
}

static double I_two_pions(double y)
{
    const double e_min = GEV_TO_EPSILON(0.4) / 2;

    if(y < e_min)
        return 0;

    // y in epsilon, x in GeV, the 2 factor is because of the equations
    double x = 2 * EPSILON_TO_GEV(y);

    return -3.4083 + 16.2864 / x + 40.7160 * log(x);
}

static double f_function_direct_one_pion(double y, double e_min, double e_max)
{
    if(y < e_min)
        return 0;

    double aux = pow(EPSILON_TO_GEV(y), -2);
    double aux1 = I_one_pion(y);
    double aux2 = I_one_pion(e_min);
    double aux3 = I_one_pion(e_max);

    if(y < e_max)
        return (aux1 - aux2) * aux;
    else
        return (aux3 - aux2) * aux;
}

static double f_function_direct_two_pions(double y, double e_min, double e_max)
{
    if(y < e_min)
        return 0;

    double aux = pow(EPSILON_TO_GEV(y), -2);
    double aux1 = I_two_pions(y);
    double aux2 = I_two_pions(e_min);
    double aux3 = I_two_pions(e_max);

    if(y < e_max)
        return (aux1 - aux2) * aux;
    else
        return (aux3 - aux2) * aux;
}

void multi_resonance_pion_production(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    double *proton_interpolated;
    /*gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_cspline, st->protons.size);*/
    gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_steffen, st->protons.size);
    gsl_interp_accel *proton_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, st->protons.log_energy, st->protons.log_population, st->protons.size);

    double *neutron_interpolated;
    /*gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_cspline, st->neutrons.size);*/
    gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_steffen, st->neutrons.size);
    gsl_interp_accel *neutron_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(neutron_spline, st->neutrons.log_energy, st->neutrons.log_population, st->neutrons.size);

    posix_memalign((void **) &proton_interpolated,  32, st->neutral_pions.size * sizeof(double));
    posix_memalign((void **) &neutron_interpolated, 32, st->neutral_pions.size * sizeof(double));

    memset(st->multi_resonances_neutral_pion_gains,  '\0', st->neutral_pions.size  * sizeof(double));
    memset(st->multi_resonances_positive_pion_gains, '\0', st->positive_pions.size * sizeof(double));
    memset(st->multi_resonances_negative_pion_gains, '\0', st->negative_pions.size * sizeof(double));

    for(process = 0; process < 18; process++)
    {
        unsigned int index_base1 = process * st->neutral_pions.size;

        double energy_ratio                       = multi_resonances_constants[process][3];

        double proton_neutral_pion_multiplicity   = multi_resonances_constants[process][4];
        double proton_positive_pion_multiplicity  = multi_resonances_constants[process][5];
        double proton_negative_pion_multiplicity  = multi_resonances_constants[process][6];

        double neutron_neutral_pion_multiplicity  = multi_resonances_constants[process][7];
        double neutron_positive_pion_multiplicity = multi_resonances_constants[process][8];
        double neutron_negative_pion_multiplicity = multi_resonances_constants[process][9];


        for(i = 0; i < st->neutral_pions.size; i++)
        {
            double g_pion = st->neutral_pions.energy[i];
            double g_proton = g_pion * PION_MASS / PROTON_MASS / energy_ratio;

            if(g_proton < st->protons.energy[0] ||
               g_proton > st->protons.energy[st->protons.size - 1])
            {
                proton_interpolated[i] = DBL_MIN;
                neutron_interpolated[i] = DBL_MIN;
            }
            else
            {
                proton_interpolated[i] = exp(gsl_spline_eval(proton_spline, log(g_proton), proton_accelerator));
                neutron_interpolated[i] = exp(gsl_spline_eval(neutron_spline, log(g_proton), neutron_accelerator));
            }
        }

        // We assume that all pions have the same mass, so their
        // gamma factor is also the same
        for(i = 0; i < st->neutral_pions.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            unsigned int index_min = st->multi_resonances_LUT_pion_gains_e_min_index[index_base1 + i];

            double production = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                double n = st->photons.population[j];

                unsigned int index = index_base2 + j;

                production += 2 * n * st->multi_resonances_LUT_pion_gains[index];
            }
            production *= dlne;

            double neutral_pion_gains_from_protons  = proton_interpolated[i] * proton_neutral_pion_multiplicity  * production;
            double positive_pion_gains_from_protons = proton_interpolated[i] * proton_positive_pion_multiplicity * production;
            double negative_pion_gains_from_protons = proton_interpolated[i] * proton_negative_pion_multiplicity * production;

            double neutral_pion_gains_from_neutrons  = neutron_interpolated[i] * neutron_neutral_pion_multiplicity  * production;
            double positive_pion_gains_from_neutrons = neutron_interpolated[i] * neutron_positive_pion_multiplicity * production;
            double negative_pion_gains_from_neutrons = neutron_interpolated[i] * neutron_negative_pion_multiplicity * production;

            st->multi_resonances_neutral_pion_gains[i]  += (neutral_pion_gains_from_protons  + neutral_pion_gains_from_neutrons);
            st->multi_resonances_positive_pion_gains[i] += (positive_pion_gains_from_protons + positive_pion_gains_from_neutrons);
            st->multi_resonances_negative_pion_gains[i] += (negative_pion_gains_from_protons + negative_pion_gains_from_neutrons);
        }
    }

    free(proton_interpolated);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(proton_accelerator);

    free(neutron_interpolated);
    gsl_spline_free(neutron_spline);
    gsl_interp_accel_free(neutron_accelerator);
}

void multi_resonance_pion_production_hadron_losses(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    memset(st->multi_resonances_proton_losses,  '\0', st->protons.size  * sizeof(double));
    memset(st->multi_resonances_neutron_losses, '\0', st->neutrons.size * sizeof(double));

    for(process = 0; process < 3; process++)
    {
        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            unsigned int index_min = st->multi_resonances_LUT_hadron_losses_e_min_index[index_base1 + i];

            double losses = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                double n = st->photons.population[j];
                losses += 2 * n * st->multi_resonances_LUT_hadron_losses[index_base2 + j];
            }

            losses *= dlne/2;

            // Model all the losses as an scape term.
            // If we don't do this, then this is only valid for proton->neutron
            // and neutron -> proton and we have to handle the case
            // proton -> proton and neutron -> neutron as cooling, which is
            // more complicated
            st->multi_resonances_proton_losses[i]  += -st->protons.population[i]  * losses;
            st->multi_resonances_neutron_losses[i] += -st->neutrons.population[i] * losses;
        }
    }
}

void multi_resonance_pion_production_hadron_gains(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    double *proton_interpolated;
    /*gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_cspline, st->protons.size);*/
    gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_steffen, st->protons.size);
    gsl_interp_accel *proton_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, st->protons.log_energy, st->protons.log_population, st->protons.size);

    double *neutron_interpolated;
    /*gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_cspline, st->neutrons.size);*/
    gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_steffen, st->neutrons.size);
    gsl_interp_accel *neutron_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(neutron_spline, st->neutrons.log_energy, st->neutrons.log_population, st->neutrons.size);

    posix_memalign((void **) &proton_interpolated,  32, st->protons.size  * sizeof(double));
    posix_memalign((void **) &neutron_interpolated, 32, st->neutrons.size * sizeof(double));

    memset(st->multi_resonances_proton_gains,  '\0', st->protons.size  * sizeof(double));
    memset(st->multi_resonances_neutron_gains, '\0', st->neutrons.size * sizeof(double));

    for(process = 0; process < 3; process++)
    {
        double inelasticity                 = multi_resonances_hadron_losses_constants[process][3];

        unsigned int index_base1 = process * st->protons.size;

        double proton_proton_multiplicity   = multi_resonances_hadron_losses_constants[process][4];
        double proton_neutron_multiplicity  = multi_resonances_hadron_losses_constants[process][5];
        double neutron_proton_multiplicity  = multi_resonances_hadron_losses_constants[process][6];
        double neutron_neutron_multiplicity = multi_resonances_hadron_losses_constants[process][7];

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            double g_proton_high = st->protons.energy[i] / (1 - inelasticity);

            if(g_proton_high < st->protons.energy[0] ||
               g_proton_high > st->protons.energy[st->protons.size - 1])
            {
                proton_interpolated[i] = DBL_MIN;
                neutron_interpolated[i] = DBL_MIN;
            }
            else
            {
                proton_interpolated[i]  = exp(gsl_spline_eval(proton_spline,  log(g_proton_high), proton_accelerator));
                neutron_interpolated[i] = exp(gsl_spline_eval(neutron_spline, log(g_proton_high), neutron_accelerator));
            }
        }

        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            unsigned int index_min = st->multi_resonances_LUT_hadron_gains_e_min_index[index_base1 + i];

            double gains = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                double n = st->photons.population[j];
                gains += 2 * n * st->multi_resonances_LUT_hadron_gains[index_base2 + j];
            }

            gains *= dlne/2;

            // -------- NOTE ---------
            // Model all the gains through the gains term including
            // proton -> proton and neutron -> neutron
            // Remember that, as we have NOT used a cooling term and
            // we have included the process p->p and n->n as an escape
            // we need to reintroduce those particles as a gains term

            double proton_gains_from_protons  = proton_interpolated[i]  * proton_proton_multiplicity  * gains;
            double proton_gains_from_neutrons = neutron_interpolated[i] * neutron_proton_multiplicity * gains;

            double neutron_gains_from_protons  = proton_interpolated[i]  * proton_neutron_multiplicity  * gains;
            double neutron_gains_from_neutrons = neutron_interpolated[i] * neutron_neutron_multiplicity * gains;

            st->multi_resonances_proton_gains[i]  += (proton_gains_from_protons + proton_gains_from_neutrons);
            st->multi_resonances_neutron_gains[i] += (neutron_gains_from_protons + neutron_gains_from_neutrons);
        }
    }

    free(proton_interpolated);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(proton_accelerator);

    free(neutron_interpolated);
    gsl_spline_free(neutron_spline);
    gsl_interp_accel_free(neutron_accelerator);
}

void direct_pion_production(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    double *proton_interpolated;
    /*gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_cspline, st->protons.size);*/
    gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_steffen, st->protons.size);
    gsl_interp_accel *proton_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, st->protons.log_energy, st->protons.log_population, st->protons.size);

    double *neutron_interpolated;
    /*gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_cspline, st->neutrons.size);*/
    gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_steffen, st->neutrons.size);
    gsl_interp_accel *neutron_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(neutron_spline, st->neutrons.log_energy, st->neutrons.log_population, st->neutrons.size);

    posix_memalign((void **) &proton_interpolated,  32, st->neutral_pions.size * sizeof(double));
    posix_memalign((void **) &neutron_interpolated, 32, st->neutral_pions.size * sizeof(double));

    memset(st->direct_neutral_pion_gains,  '\0', st->neutral_pions.size  * sizeof(double));
    memset(st->direct_positive_pion_gains, '\0', st->positive_pions.size * sizeof(double));
    memset(st->direct_negative_pion_gains, '\0', st->negative_pions.size * sizeof(double));

    for(process = 0; process < 7; process++)
    {
        unsigned int index_base1 = process * st->neutral_pions.size;

        double energy_ratio                       = direct_pion_production_constants[process][2];

        double proton_neutral_pion_multiplicity   = direct_pion_production_constants[process][3];
        double proton_positive_pion_multiplicity  = direct_pion_production_constants[process][4];
        double proton_negative_pion_multiplicity  = direct_pion_production_constants[process][5];

        double neutron_neutral_pion_multiplicity  = direct_pion_production_constants[process][6];
        double neutron_positive_pion_multiplicity = direct_pion_production_constants[process][7];
        double neutron_negative_pion_multiplicity = direct_pion_production_constants[process][8];

        for(i = 0; i < st->neutral_pions.size; i++)
        {
            double g_pion = st->neutral_pions.energy[i];
            double g_proton = g_pion * PION_MASS / PROTON_MASS / energy_ratio;

            if(g_proton < st->protons.energy[0] ||
               g_proton > st->protons.energy[st->protons.size - 1])
            {
                proton_interpolated[i] = DBL_MIN;
                neutron_interpolated[i] = DBL_MIN;
            }
            else
            {
                proton_interpolated[i] = exp(gsl_spline_eval(proton_spline, log(g_proton), proton_accelerator));
                neutron_interpolated[i] = exp(gsl_spline_eval(neutron_spline, log(g_proton), neutron_accelerator));
            }
        }

        // We assume that all pions have the same mass, so their
        // gamma factor is also the same
        for(i = 0; i < st->neutral_pions.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;
            unsigned int index_min = st->direct_LUT_pion_gains_e_min_index[index_base1 + i];

            double production = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double n = st->photons.population[j];

                production += 2 * n * st->direct_LUT_pion_gains[index];
            }
            production *= dlne;

            double neutral_pion_gains_from_protons  = proton_interpolated[i] * proton_neutral_pion_multiplicity  * production;
            double positive_pion_gains_from_protons = proton_interpolated[i] * proton_positive_pion_multiplicity * production;
            double negative_pion_gains_from_protons = proton_interpolated[i] * proton_negative_pion_multiplicity * production;

            double neutral_pion_gains_from_neutrons  = neutron_interpolated[i] * neutron_neutral_pion_multiplicity  * production;
            double positive_pion_gains_from_neutrons = neutron_interpolated[i] * neutron_positive_pion_multiplicity * production;
            double negative_pion_gains_from_neutrons = neutron_interpolated[i] * neutron_negative_pion_multiplicity * production;

            st->direct_neutral_pion_gains[i]  += (neutral_pion_gains_from_protons  + neutral_pion_gains_from_neutrons);
            st->direct_positive_pion_gains[i] += (positive_pion_gains_from_protons + positive_pion_gains_from_neutrons);
            st->direct_negative_pion_gains[i] += (negative_pion_gains_from_protons + negative_pion_gains_from_neutrons);
        }
    }

    free(proton_interpolated);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(proton_accelerator);

    free(neutron_interpolated);
    gsl_spline_free(neutron_spline);
    gsl_interp_accel_free(neutron_accelerator);
}

void direct_pion_production_hadron_losses(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    memset(st->direct_pion_production_proton_losses,  '\0', st->protons.size  * sizeof(double));
    memset(st->direct_pion_production_neutron_losses, '\0', st->neutrons.size * sizeof(double));

    for(process = 0; process < 6; process++)
    {
        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            unsigned int index_min = st->direct_pion_production_LUT_hadron_losses_e_min_index[index_base1 + i];

            double losses = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double n = st->photons.population[j];

                losses += 2 * n * st->direct_pion_production_LUT_hadron_losses[index];
            }
            losses *= dlne/2;

            // Model all the losses as an scape term.
            // If we don't do this, then this is only valid for proton->neutron
            // and neutron -> proton and we have to handle the case
            // proton -> proton and neutron -> neutron as cooling, which is
            // more complicated
            st->direct_pion_production_proton_losses[i]  += -st->protons.population[i]  * losses;
            st->direct_pion_production_neutron_losses[i] += -st->neutrons.population[i] * losses;
        }
    }
}

void direct_pion_production_hadron_gains(state_t *st)
{
    unsigned int i, j;
    unsigned int process;

    double dlne = st->photons.log_energy[1] - st->photons.log_energy[0];

    double *proton_interpolated;
    /*gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_cspline, st->protons.size);*/
    gsl_spline       *proton_spline       = gsl_spline_alloc(gsl_interp_steffen, st->protons.size);
    gsl_interp_accel *proton_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(proton_spline, st->protons.log_energy, st->protons.log_population, st->protons.size);

    double *neutron_interpolated;
    /*gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_cspline, st->neutrons.size);*/
    gsl_spline       *neutron_spline       = gsl_spline_alloc(gsl_interp_steffen, st->neutrons.size);
    gsl_interp_accel *neutron_accelerator  = gsl_interp_accel_alloc();
    gsl_spline_init(neutron_spline, st->neutrons.log_energy, st->neutrons.log_population, st->neutrons.size);

    posix_memalign((void **) &proton_interpolated,  32, st->protons.size  * sizeof(double));
    posix_memalign((void **) &neutron_interpolated, 32, st->neutrons.size * sizeof(double));

    memset(st->direct_pion_production_proton_gains,  '\0', st->protons.size  * sizeof(double));
    memset(st->direct_pion_production_neutron_gains, '\0', st->neutrons.size * sizeof(double));

    for(process = 0; process < 6; process++)
    {
        unsigned int index_base1 = process * st->protons.size;

        double inelasticity                 = direct_pion_production_hadron_losses_constants[process][2];

        double proton_proton_multiplicity   = direct_pion_production_hadron_losses_constants[process][3];
        double proton_neutron_multiplicity  = direct_pion_production_hadron_losses_constants[process][4];
        double neutron_proton_multiplicity  = direct_pion_production_hadron_losses_constants[process][5];
        double neutron_neutron_multiplicity = direct_pion_production_hadron_losses_constants[process][6];

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            double g_proton_high = st->protons.energy[i] / (1 - inelasticity);

            if(g_proton_high < st->protons.energy[0] ||
               g_proton_high > st->protons.energy[st->protons.size - 1])
            {
                proton_interpolated[i] = DBL_MIN;
                neutron_interpolated[i] = DBL_MIN;
            }
            else
            {
                proton_interpolated[i]  = exp(gsl_spline_eval(proton_spline,  log(g_proton_high), proton_accelerator));
                neutron_interpolated[i] = exp(gsl_spline_eval(neutron_spline, log(g_proton_high), neutron_accelerator));
            }
        }

        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            unsigned int index_min = st->direct_pion_production_LUT_hadron_gains_e_min_index[index_base1 + i];

            double gains = 0;
            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double n = st->photons.population[j];

                gains += 2 * n * st->direct_pion_production_LUT_hadron_gains[index];
            }
            gains *= dlne/2;

            // -------- NOTE ---------
            // Model all the gains through the gains term including
            // proton -> proton and neutron -> neutron
            // Remember that, as we have NOT used a cooling term and
            // we have included the process p->p and n->n as an escape
            // we need to reintroduce those particles as a gains term

            double proton_gains_from_protons  = proton_interpolated[i]  * proton_proton_multiplicity  * gains;
            double proton_gains_from_neutrons = neutron_interpolated[i] * neutron_proton_multiplicity * gains;

            double neutron_gains_from_protons  = proton_interpolated[i]  * proton_neutron_multiplicity  * gains;
            double neutron_gains_from_neutrons = neutron_interpolated[i] * neutron_neutron_multiplicity * gains;

            st->direct_pion_production_proton_gains[i]  += (proton_gains_from_protons + proton_gains_from_neutrons);
            st->direct_pion_production_neutron_gains[i] += (neutron_gains_from_protons + neutron_gains_from_neutrons);
        }
    }

    free(proton_interpolated);
    gsl_spline_free(proton_spline);
    gsl_interp_accel_free(proton_accelerator);

    free(neutron_interpolated);
    gsl_spline_free(neutron_spline);
    gsl_interp_accel_free(neutron_accelerator);
}

void init_multi_resonances_LUT_pion_gains(state_t *st)
{
    size_t size = 18 * st->neutral_pions.size * st->photons.size;

    posix_memalign((void **) &st->multi_resonances_LUT_pion_gains, 32, sizeof(double) * size);

    posix_memalign((void **) &st->multi_resonances_LUT_pion_gains_e_min_index, 32, sizeof(double) * 18 * st->neutral_pions.size);
}

void init_multi_resonances_LUT_hadron_losses(state_t *st)
{
    size_t size = 3 * st->protons.size * st->photons.size;

    posix_memalign((void **) &st->multi_resonances_LUT_hadron_losses, 32, sizeof(double) * size);

    posix_memalign((void **) &st->multi_resonances_LUT_hadron_losses_e_min_index, 32, sizeof(double) * 3 * st->protons.size);
}

void init_multi_resonances_LUT_hadron_gains(state_t *st)
{
    size_t size = 3 * st->protons.size * st->photons.size;

    posix_memalign((void **) &st->multi_resonances_LUT_hadron_gains, 32, sizeof(double) * size);

    posix_memalign((void **) &st->multi_resonances_LUT_hadron_gains_e_min_index, 32, sizeof(double) * 3 * st->protons.size);
}

void init_direct_LUT_pion_gains(state_t *st)
{
    size_t size = 7 * st->neutral_pions.size * st->photons.size;

    posix_memalign((void **) &st->direct_LUT_pion_gains, 32, sizeof(double) * size);

    posix_memalign((void **) &st->direct_LUT_pion_gains_e_min_index, 32, sizeof(double) * 18 * st->neutral_pions.size);
}

void init_direct_pion_production_LUT_hadron_losses(state_t *st)
{
    size_t size = 6 * st->protons.size * st->photons.size;

    posix_memalign((void **) &st->direct_pion_production_LUT_hadron_losses, 32, sizeof(double) * size);

    posix_memalign((void **) &st->direct_pion_production_LUT_hadron_losses_e_min_index, 32, sizeof(double) * 6 * st->protons.size);
}

void init_direct_pion_production_LUT_hadron_gains(state_t *st)
{
    size_t size = 6 * st->protons.size * st->photons.size;

    posix_memalign((void **) &st->direct_pion_production_LUT_hadron_gains, 32, sizeof(double) * size);

    posix_memalign((void **) &st->direct_pion_production_LUT_hadron_gains_e_min_index, 32, sizeof(double) * 6 * st->protons.size);
}

void calculate_multi_resonances_LUT_pion_gains(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 18; process++)
    {
        double low_threshold  = multi_resonances_constants[process][0];
        double high_threshold = multi_resonances_constants[process][1];
        double cross_section  = multi_resonances_constants[process][2];
        double energy_ratio   = multi_resonances_constants[process][3];

        unsigned int index_base1 = process * st->neutral_pions.size;

        // We assume that all pions have the same mass, so their
        // gamma factor is also the same
        for(i = 0; i < st->neutral_pions.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_pion = st->neutral_pions.energy[i];
            double g_proton = g_pion * PION_MASS / PROTON_MASS / energy_ratio;

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->multi_resonances_LUT_pion_gains_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;
                double f = f_function(y, cross_section, low_threshold, high_threshold);

                st->multi_resonances_LUT_pion_gains[index] = f * e / energy_ratio;
            }
        }
    }
}

void calculate_multi_resonances_LUT_hadron_losses(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 3; process++)
    {
        double low_threshold                = multi_resonances_hadron_losses_constants[process][0];
        double high_threshold               = multi_resonances_hadron_losses_constants[process][1];
        double cross_section                = multi_resonances_hadron_losses_constants[process][2];

        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_proton = st->protons.energy[i];

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->multi_resonances_LUT_hadron_losses_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;

                double f = (process == 0 ?
                            f_function_multi_resonance(y, low_threshold, high_threshold) * 1e-30:
                            f_function(y, cross_section, low_threshold, high_threshold));

                st->multi_resonances_LUT_hadron_losses[index] = f * e;
            }
        }
    }
}

void calculate_multi_resonances_LUT_hadron_gains(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 3; process++)
    {
        double low_threshold                = multi_resonances_hadron_losses_constants[process][0];
        double high_threshold               = multi_resonances_hadron_losses_constants[process][1];
        double cross_section                = multi_resonances_hadron_losses_constants[process][2];
        double inelasticity                 = multi_resonances_hadron_losses_constants[process][3];

        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_proton = st->protons.energy[i] / (1 - inelasticity);

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->multi_resonances_LUT_hadron_gains_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;

                double f = (process == 0 ?
                            f_function_multi_resonance(y, low_threshold, high_threshold) * 1e-30:
                            f_function(y, cross_section, low_threshold, high_threshold));

                st->multi_resonances_LUT_hadron_gains[index] = f * e / (1 - inelasticity);
            }
        }
    }
}

void calculate_direct_LUT_pion_gains(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 7; process++)
    {
        double low_threshold  = direct_pion_production_constants[process][0];
        double high_threshold = direct_pion_production_constants[process][1];
        double energy_ratio   = direct_pion_production_constants[process][2];

        unsigned int index_base1 = process * st->neutral_pions.size;

        // We assume that all pions have the same mass, so their
        // gamma factor is also the same
        for(i = 0; i < st->neutral_pions.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_pion = st->neutral_pions.energy[i];
            double g_proton = g_pion * PION_MASS / PROTON_MASS / energy_ratio;

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->direct_LUT_pion_gains_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;
                double f = (process < 3 ?
                                f_function_direct_one_pion (y, low_threshold, high_threshold) :
                                f_function_direct_two_pions(y, low_threshold, high_threshold));

                st->direct_LUT_pion_gains[index] = f * e * 1e-30 / energy_ratio;
            }
        }
    }
}

void calculate_direct_pion_production_LUT_hadron_losses(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 6; process++)
    {
        double low_threshold  = direct_pion_production_hadron_losses_constants[process][0];
        double high_threshold = direct_pion_production_hadron_losses_constants[process][1];

        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_proton = st->protons.energy[i];

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->direct_pion_production_LUT_hadron_losses_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;

                double f = (process < 3 ?
                                f_function_direct_one_pion (y, low_threshold, high_threshold) :
                                f_function_direct_two_pions(y, low_threshold, high_threshold));

                st->direct_pion_production_LUT_hadron_losses[index] = f * e * 1e-30;
            }
        }
    }
}

void calculate_direct_pion_production_LUT_hadron_gains(state_t *st)
{
    unsigned int process, i, j;

    for(process = 0; process < 6; process++)
    {
        double low_threshold  = direct_pion_production_hadron_losses_constants[process][0];
        double high_threshold = direct_pion_production_hadron_losses_constants[process][1];
        double inelasticity   = direct_pion_production_hadron_losses_constants[process][2];

        unsigned int index_base1 = process * st->protons.size;

        // We assume that protons and neutrons have the same mass,
        // same energy limits and same size
        for(i = 0; i < st->protons.size; i++)
        {
            unsigned int index_base2 = (i + index_base1) * st->photons.size;

            double g_proton = st->protons.energy[i] / (1 - inelasticity);

            double e_min = low_threshold / g_proton;
            unsigned int index_min;
            for(index_min = 0; index_min < st->photons.size && st->photons.energy[index_min] < e_min; index_min++)
            { }

            st->direct_pion_production_LUT_hadron_gains_e_min_index[index_base1 + i] = index_min;

            for(j = index_min; j < st->photons.size; j++)
            {
                unsigned int index = index_base2 + j;

                double e = st->photons.energy[j];
                double y = e * g_proton;

                double f = (process < 3 ?
                                f_function_direct_one_pion (y, low_threshold, high_threshold) :
                                f_function_direct_two_pions(y, low_threshold, high_threshold));

                st->direct_pion_production_LUT_hadron_gains[index] = f * e * 1e-30;
            }
        }
    }
}
