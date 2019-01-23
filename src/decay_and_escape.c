#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "constants.h"
#include "decay_and_escape.h"
#include "state.h"

static void free_decay_and_escape(decay_and_escape_t *dea);

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void init_decay_and_escape(struct state_t *st, decay_and_escape_t *dea, enum particle_type pt, double escape_lifetime)
{
    switch(pt)
    {
        case neutron:
            dea->particles = &st->neutrons;
            assert(0); // Not implemented
            break;

        case positive_pion: dea->particles = &st->positive_pions; break;
        case neutral_pion:  dea->particles = &st->neutral_pions;  break;
        case negative_pion: dea->particles = &st->negative_pions; break;

        case positive_left_muon:  dea->particles = &st->positive_left_muons;  break;
        case positive_right_muon: dea->particles = &st->positive_right_muons; break;
        case negative_left_muon:  dea->particles = &st->negative_left_muons;  break;
        case negative_right_muon: dea->particles = &st->negative_right_muons; break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    switch(pt)
    {
        case neutron:
            dea->decay_lifetime = 0;
            assert(0); // Not implemented
            break;

        case neutral_pion:
            dea->decay_lifetime = NEUTRAL_PION_LIFETIME;

        case positive_pion:
        case negative_pion:
            dea->decay_lifetime = CHARGED_PION_LIFETIME;
            break;

        case positive_left_muon:
        case positive_right_muon:
        case negative_left_muon:
        case negative_right_muon:
            dea->decay_lifetime = MUON_LIFETIME;
            break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    MEM_PREPARE(dea->losses,        dea->particles->size);
    MEM_PREPARE(dea->losses_factor, dea->particles->size);

    dea->losses_function = free_decay_and_escape;

    update_decay_and_escape(st, dea, escape_lifetime);
}

void update_decay_and_escape(struct state_t *st, decay_and_escape_t *dea, double escape_lifetime)
{
    dea->escape_lifetime = escape_lifetime;
    dea->dt = st->dt;

    unsigned int i;
    for(i = 0; i < dea->particles->size; i++)
    {
        double t = (dea->particles->energy[i] * dea->decay_lifetime + dea->escape_lifetime) / 
                   (dea->particles->energy[i] * dea->decay_lifetime * dea->escape_lifetime);
        dea->losses_factor[i] = expm1(-dea->dt / t) / dea->dt;
    }
}

static void free_decay_and_escape(decay_and_escape_t *dea)
{
    unsigned int i;
    population_t *p = dea->particles;

    for(i = 0; i < p->size; i++)
        dea->losses[i] = p->population[i] * dea->losses_factor[i];
}
