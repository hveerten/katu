#include "constants.h"
#include "synchrotron.h"
#include "state.h"

#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#include <stdio.h>
#include <string.h>

#include <float.h>

#include <assert.h>

#define MAX(x,y) ((x) < (y) ? (y) : (x))

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);


static void synchrotron_process_full(synchrotron_t *synchro);
static void init_synchrotron_LUT_xCS(synchrotron_t *synchro);
static void calculate_synchrotron_LUT_xCS(synchrotron_t *synchro);

#if USE_THREADS
void *synchrotron_process_full_wrapper(void *args)
{
    synchrotron_t *s = (synchrotron_t *)args;
    synchrotron_process_full(s);
    return NULL;
}
#endif

void init_synchrotron(state_t *st, synchrotron_t *synchro, enum particle_type pt)
{
    synchro->photons = &st->photons;

    switch(pt)
    {
        case electron:      synchro->particles = &st->electrons;      break;
        case positron:      synchro->particles = &st->positrons;      break;
        case proton:        synchro->particles = &st->protons;        break;
        case positive_pion: synchro->particles = &st->positive_pions; break;
        case negative_pion: synchro->particles = &st->negative_pions; break;

        case positive_left_muon:  synchro->particles = &st->positive_left_muons;  break;
        case positive_right_muon: synchro->particles = &st->positive_right_muons; break;
        case negative_left_muon:  synchro->particles = &st->negative_left_muons;  break;
        case negative_right_muon: synchro->particles = &st->negative_right_muons; break;

        // TODO: ADD FAIL MESSAGES
        case photon:   assert(1);   break;
        case neutron:   assert(1);   break;
        case neutral_pion:   assert(1);   break;
        case electron_neutrino:   assert(1);   break;
        case electron_antineutrino:   assert(1);   break;
        case muon_neutrino:   assert(1);   break;
        case muon_antineutrino:   assert(1);   break;
        default: assert(1); break;
    }

    synchro->nu_0 = 3 / (4 * M_PI) * ELECTRON_CHARGE * st->B / (synchro->particles->mass * LIGHT_SPEED);
    double magnetic_energy = st->B * st->B / (8 * M_PI);

    synchro->photon_gains_factor    = sqrt(3) * M_PI /  3 * FINE_STRUCTURE_CONSTANT * synchro->nu_0;
    synchro->photon_losses_factor   = sqrt(3) * M_PI / 12 * ELECTRON_RADIUS * ELECTRON_ENERGY * synchro->nu_0 / synchro->particles->mass;
    synchro->particle_losses_factor = 4 * LIGHT_SPEED * THOMSON_CROSS_SECTION * magnetic_energy / (3 * ELECTRON_ENERGY) * pow(ELECTRON_MASS / synchro->particles->mass, 3);

    synchro->dt = st->dt;

    MEM_PREPARE(synchro->temp1,           synchro->particles->size);
    MEM_PREPARE(synchro->temp2,           synchro->particles->size);
    MEM_PREPARE(synchro->photon_gains,    synchro->photons->size);
    MEM_PREPARE(synchro->photon_losses,   synchro->photons->size);
    MEM_PREPARE(synchro->particle_losses, synchro->particles->size);

    init_synchrotron_LUT_xCS(synchro);
    calculate_synchrotron_LUT_xCS(synchro);
#if USE_THREADS
    synchro->synchrotron_function = &synchrotron_process_full_wrapper;
#else
    synchro->synchrotron_function = &synchrotron_process_full;
#endif
}

void synchrotron_update_dt(synchrotron_t *synchro, double dt)
{
    synchro->dt = dt;
}

void synchrotron_update_B(synchrotron_t *synchro, double B)
{
    synchro->nu_0 = 3 / (4 * M_PI) * ELECTRON_CHARGE * B / (synchro->particles->mass * LIGHT_SPEED);
    double magnetic_energy = B * B / (8 * M_PI);

    synchro->photon_gains_factor    = sqrt(3) * M_PI /  3 * FINE_STRUCTURE_CONSTANT * synchro->nu_0;
    synchro->photon_losses_factor   = sqrt(3) * M_PI / 12 * ELECTRON_RADIUS * ELECTRON_ENERGY * synchro->nu_0 / synchro->particles->mass;
    synchro->particle_losses_factor = 4 * LIGHT_SPEED * THOMSON_CROSS_SECTION * magnetic_energy / (3 * ELECTRON_ENERGY) * pow(ELECTRON_MASS / synchro->particles->mass, 3);

    calculate_synchrotron_LUT_xCS(synchro);
}


static double xCS(double x)
{
    return x * x / M_PI * gsl_sf_bessel_Knu(4./3, x / 2) * gsl_sf_bessel_Knu(1./3, x / 2) - \
           pow(x, 11./3) * exp(-x) * gsl_sf_hyperg_U(5./6, 8./3, x) * gsl_sf_hyperg_U(11./6, 8./3, x);
}

static double xCS_approx(double x)
{
    if(x < 0.291)
        return (29 * cbrt(x) - 27 * x) / 25;
    else if(x < 2.7675)
        return (135 - 77 * x + 12 * x * x) / 250;
    else
        return exp(-x) * (1 - 1 / (3 * x));
}

static void synchrotron_process_full(synchrotron_t *synchro)
{
    population_t *part = synchro->particles;
    population_t *phot = synchro->photons;
    unsigned int i, j;

    double dlng = part->log_energy[1] - part->log_energy[0];
    for(i = 0; i < part->size - 1; i++)
        synchro->temp1[i] = (part->log_population[i + 1] - part->log_population[i]) / dlng;

    for(i = 0; i < part->size - 1; i++)
        synchro->temp2[i] = part->population[i] * (synchro->temp1[i] - 2);
    synchro->temp2[i] = synchro->temp2[i - 1];


    for(i = 0; i < phot->size; i++)
    {
        double e = phot->energy[i];
        double nu = ELECTRON_ENERGY * e / PLANCK_CONSTANT;
        double nu_nu0 = nu / synchro->nu_0;

        double g_min = MAX(part->energy[0], 0.53 / nu_nu0);
        double g_max = part->energy[part->size - 1];

        /*fprintf(stderr,"%u:\t%lg\t%lg\n",i,g_min,g_max);*/

        unsigned int index_of_g_min = 0; // TODO
        unsigned int index_of_g_max = part->size - 1;

        if(g_min >= g_max) continue;
        for(/* START */; part->energy[index_of_g_min] < g_min; index_of_g_min++)
        {
        }

        double emission = 0;
        double absorption = 0;
        { // Integral, to be separated into a function
            double g0 = part->energy[index_of_g_min];
            double gN = part->energy[index_of_g_max];
            double n0 = part->population[index_of_g_min];
            double nN = part->population[index_of_g_max];

            double dlng = part->log_energy[1] - part->log_energy[0];

/*
 *            emission += g0 * n0 * xCS_approx(nu_nu0 / (g0 * g0));
 *            emission += g0 * nN * xCS_approx(nu_nu0 / (gN * gN));
 *
 *            absorption += xCS_approx(nu_nu0 / (g0 * g0)) * synchro->temp2[index_of_g_min];
 *            absorption += xCS_approx(nu_nu0 / (gN * gN)) * synchro->temp2[index_of_g_max];
 */
            unsigned int index_min = i * part->size + index_of_g_min;
            unsigned int index_max = i * part->size + index_of_g_max;
            emission += g0 * n0 * synchro->LUT_xCS[index_min];
            emission += gN * nN * synchro->LUT_xCS[index_max];

            absorption += synchro->LUT_xCS[index_min] * synchro->temp2[index_of_g_min];
            absorption += synchro->LUT_xCS[index_max] * synchro->temp2[index_of_g_max];

            for(j = index_of_g_min + 1; j < index_of_g_max; j++)
            {
                double g = part->energy[j];
                double n = part->population[j];

                unsigned int index = i * part->size + j;

                /*
                 *emission   += 2 * g * n * xCS_approx(nu_nu0 / (g * g));
                 *absorption += 2 *         xCS_approx(nu_nu0 / (g * g)) * synchro->temp2[j];
                 */
                emission   += 2 * g * n * synchro->LUT_xCS[index];
                absorption += 2 *         synchro->LUT_xCS[index] * synchro->temp2[j];
            }

            emission   *= dlng / 2;
            absorption *= dlng / 2;
        }

        synchro->photon_gains[i]  = synchro->photon_gains_factor  * emission   / e;
        synchro->photon_losses[i] = synchro->photon_losses_factor * absorption / (nu * nu) * phot->population[i];
    }

    for(i = 0; i < part->size; i++)
    {
        synchro->particle_losses[i] = synchro->particle_losses_factor *
            part->population[i] * part->energy[i] *
            ((1 - 1 / (part->energy[i] * part->energy[i])) * synchro->temp1[i] + 2);
    }
    synchro->particle_losses[part->size - 1] = 0;

    double log_n_new = part->log_population[part->size - 1];
    for(i = part->size - 2; i < part->size; i--)
    {
        double aux1 = synchro->dt * part->energy[i] * synchro->particle_losses_factor / dlng;
        double aux2 = 1 - 1 / (part->energy[i] * part->energy[i]);

        log_n_new = (part->log_population[i] + aux1 * (aux2 * log_n_new + 2 * dlng)) / (1 + aux1 * aux2);

        if(log_n_new < log(DBL_MIN))
        {
            log_n_new = log(DBL_MIN);
            synchro->particle_losses[i] = (DBL_MIN - part->population[i]) / synchro->dt;
        }
        else
            synchro->particle_losses[i] = (exp(log_n_new) - part->population[i]) / synchro->dt;
    }
    synchro->particle_losses[part->size - 1] = 0;
}

static void init_synchrotron_LUT_xCS(synchrotron_t *synchro)
{
    size_t size = synchro->particles->size * synchro->photons->size;

    MEM_PREPARE(synchro->LUT_xCS, size);
}

static void calculate_synchrotron_LUT_xCS(synchrotron_t *synchro)
{
    unsigned int i, j;

    gsl_set_error_handler_off();
    for(i = 0; i < synchro->photons->size; i++)
    {
        double e = synchro->photons->energy[i];
        double nu  = ELECTRON_ENERGY * e / PLANCK_CONSTANT;

        for(j = 0; j < synchro->particles->size; j++)
        {
            unsigned int index = i * synchro->particles->size + j;
            double g = synchro->particles->energy[j];

            synchro->LUT_xCS[index] = xCS(nu / (synchro->nu_0 * g * g));
        }
    }
}
