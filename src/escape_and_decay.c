#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "escape_and_decay.h"
#include "state.h"
#include "constants.h"

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void init_escape_and_decay(struct state_t *st, escape_and_decay_t *ead, enum particle_type pt, double escape_lifetime)
{
    switch(pt)
    {
        case neutron:
            ead->particles = &st->neutrons;
            assert(0); // Not implemented
            break;

        case positive_pion: ead->particles = &st->positive_pions; break;
        case neutral_pion:  ead->particles = &st->neutral_pions;  break;
        case negative_pion: ead->particles = &st->negative_pions; break;

        case positive_left_muon:  ead->particles = &st->positive_left_muons;  break;
        case positive_right_muon: ead->particles = &st->positive_right_muons; break;
        case negative_left_muon:  ead->particles = &st->negative_left_muons;  break;
        case negative_right_muon: ead->particles = &st->negative_right_muons; break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    switch(pt)
    {
        case neutron:
            ead->decay_lifetime = 0;
            assert(0); // Not implemented
            break;

        case positive_pion:
        case neutral_pion:
        case negative_pion:
            ead->decay_lifetime = CHARGED_PION_LIFETIME;
            break;

        case positive_left_muon:
        case positive_right_muon:
        case negative_left_muon:
        case negative_right_muon:
            ead->decay_lifetime = MUON_LIFETIME;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    MEM_PREPARE(ead->losses,        ead->particles->size);
    MEM_PREPARE(ead->losses_factor, ead->particles->size);

    ead->losses_function = free_escape_and_decay;

    update_escape_and_decay(st, ead, escape_lifetime);
}

void update_escape_and_decay(struct state_t *st, escape_and_decay_t *ead, double escape_lifetime)
{
    ead->escape_lifetime = escape_lifetime;
    ead->dt = st->dt;

    unsigned int i;
    for(i = 0; i < ead->particles->size; i++)
    {
        double t = (ead->particles->energy[i] * ead->decay_lifetime + ead->escape_lifetime) / 
                   (ead->particles->energy[i] * ead->decay_lifetime * ead->escape_lifetime);
        ead->losses_factor[i] = expm1(-ead->dt / t) / ead->dt;
    }
}

void free_escape_and_decay(escape_and_decay_t *ead)
{
    unsigned int i;
    population_t *p = ead->particles;

    for(i = 0; i < p->size; i++)
        ead->losses[i] = p->population[i] * ead->losses_factor[i];
}
