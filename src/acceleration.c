#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "acceleration.h"
#include "state.h"

#include <stdio.h>

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void linear_acceleration_implicit(acceleration_t *acc);
void quadratic_acceleration_implicit(acceleration_t *acc);

void init_acceleration(state_t *st, acceleration_t *acc, enum particle_type pt, double t)
{
    switch(pt)
    {
        case electron:      acc->particles = &st->electrons;      break;
        case photon:        acc->particles = &st->photons;        break;
        case proton:        acc->particles = &st->protons;        break;
        case neutron:       acc->particles = &st->neutrons;       break;

        case positive_pion: acc->particles = &st->positive_pions; break;
        case neutral_pion:  acc->particles = &st->neutral_pions;  break;
        case negative_pion: acc->particles = &st->negative_pions; break;

        case positive_left_muon:  acc->particles = &st->positive_left_muons;  break;
        case positive_right_muon: acc->particles = &st->positive_right_muons; break;
        case negative_left_muon:  acc->particles = &st->negative_left_muons;  break;
        case negative_right_muon: acc->particles = &st->negative_right_muons; break;

        case electron_neutrino:     acc->particles = &st->electron_neutrinos;       break;
        case electron_antineutrino: acc->particles = &st->electron_antineutrinos;   break;
        case muon_neutrino:         acc->particles = &st->muon_neutrinos;           break;
        case muon_antineutrino:     acc->particles = &st->muon_antineutrinos;       break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }
    MEM_PREPARE(acc->gains, acc->particles->size);
    MEM_PREPARE(acc->temp1, acc->particles->size);
    MEM_PREPARE(acc->temp2, acc->particles->size);

    acc->t = t;
    acc->dt = st->dt;

    acc->acceleration_function = linear_acceleration_implicit;
    /*acc->acceleration_function = quadratic_acceleration_implicit;*/
}

void update_acceleration(state_t *st, acceleration_t *acc, double t)
{
    acc->t = t;
    acc->dt = st->dt;
}

void linear_acceleration(acceleration_t *acc)
{
    unsigned int i;
    population_t *p = acc->particles;
    double dlng = p->log_energy[1] - p->log_energy[0];

    for(i = 1; i < p->size; i++)
        acc->temp1[i] = (p->log_population[i] - p->log_population[i - 1]) / dlng;

    for(i = 1; i < p->size; i++)
        acc->temp2[i] = p->population[i] * (acc->temp1[i] + 1);

    // Temp2 now holds partial(g*n)/partial(g)

    for(i = 0; i < p->size; i++)
        acc->gains[i] = -acc->temp2[i] / acc->t;
    acc->gains[0] = 0;
}

void linear_acceleration_implicit(acceleration_t *acc)
{
    unsigned int i;
    population_t *p = acc->particles;
    double dlng = p->log_energy[1] - p->log_energy[0];
    double aux = acc->dt / dlng / acc->t;

    double ln_n_new = p->log_population[0];
    /*
     *fprintf(stderr,"%u:\t%lg\n",0,ln_n_new);
     *fprintf(stderr,"%u:\t%lg\n",0,aux);
     *fprintf(stderr,"%u:\t%lg\n",0,acc->t);
     *fprintf(stderr,"%u:\t%lg\n",0,dlng);
     */
    for(i = 1; i < p->size; i++)
    {
        ln_n_new = (p->log_population[i] - aux * (dlng - ln_n_new)) /
                   (1 + aux);

        acc->gains[i] = (exp(ln_n_new) - p->population[i]) / acc->dt;

        /*
         *fprintf(stderr,"%u:\t%lg\n",i,ln_n_new);
         *fprintf(stderr,"%u:\t%lg\n",i,p->energy[i]);
         */
    }

    acc->gains[0] = 0;
}

void quadratic_acceleration(acceleration_t *acc)
{
    unsigned int i;
    population_t *p = acc->particles;
    double dlng = p->log_energy[1] - p->log_energy[0];

    for(i = 1; i < p->size; i++)
        acc->temp1[i] = (p->log_population[i] - p->log_population[i - 1]) / dlng;

    for(i = 1; i < p->size; i++)
        acc->temp2[i] = p->population[i] * p->energy[i] * (acc->temp1[i] + 2);

    // Temp2 now holds partial(g^2*n)/partial(g)

    for(i = 0; i < p->size; i++)
        acc->gains[i] = -acc->temp2[i] / acc->t;
    acc->gains[0] = 0;
}

void quadratic_acceleration_implicit(acceleration_t *acc)
{
    unsigned int i;
    population_t *p = acc->particles;
    double dlng = p->log_energy[1] - p->log_energy[0];
    double aux = acc->dt / dlng / acc->t;

    double ln_n_new = p->log_population[0];
    for(i = 1; i < p->size; i++)
    {
        double aux2 = aux * p->energy[i];
        ln_n_new = (p->log_population[i] - aux2 * (2 * dlng - ln_n_new)) /
                   (1 + aux2);

        acc->gains[i] = (exp(ln_n_new) - p->population[i]) / acc->dt;
    }

    acc->gains[0] = 0;
}
