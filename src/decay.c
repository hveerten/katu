#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "constants.h"
#include "decay.h"
#include "state.h"

static void decay(decay_t *d);

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void init_decay(struct state_t *st, decay_t *d, enum particle_type pt)
{
    switch(pt)
    {
        case neutron: d->particles = &st->neutrons; break;

        case positive_pion: d->particles = &st->positive_pions; break;
        case neutral_pion:  d->particles = &st->neutral_pions;  break;
        case negative_pion: d->particles = &st->negative_pions; break;

        case positive_left_muon:  d->particles = &st->positive_left_muons;  break;
        case positive_right_muon: d->particles = &st->positive_right_muons; break;
        case negative_left_muon:  d->particles = &st->negative_left_muons;  break;
        case negative_right_muon: d->particles = &st->negative_right_muons; break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    switch(pt)
    {
        case neutron: d->t = NEUTRON_LIFETIME; break;

        case positive_pion:
        case negative_pion:
            d->t = CHARGED_PION_LIFETIME;
            break;

        case neutral_pion: d->t = NEUTRAL_PION_LIFETIME; break;

        case positive_left_muon:
        case positive_right_muon:
        case negative_left_muon:
        case negative_right_muon:
            d->t = MUON_LIFETIME;
            break;

        // Unreachable due to previous switch
        default: assert(0); break;
    }

    MEM_PREPARE(d->losses,       d->particles->size);
    MEM_PREPARE(d->decay_factor, d->particles->size);

    d->decay_function = decay;

    update_decay(st, d);
}

void update_decay(struct state_t *st, decay_t *d)
{
    d->dt = st->dt;

    unsigned int i;
    population_t *p = d->particles;

    for(i = 0; i < p->size; i++)
        d->decay_factor[i] = expm1(-d->dt / (p->energy[i] * d->t)) / d->dt;
}

static void decay(decay_t *d)
{
    unsigned int i;
    population_t *p = d->particles;

    for(i = 0; i < p->size; i++)
        d->losses[i] = p->population[i] * d->decay_factor[i];
}
