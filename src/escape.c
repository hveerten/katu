#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "escape.h"
#include "state.h"

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void init_escape(struct state_t *st, escape_t *esc, enum particle_type pt, double t)
{
    switch(pt)
    {
        case electron:      esc->particles = &st->electrons;      break;
        case photon:        esc->particles = &st->photons;        break;
        case proton:        esc->particles = &st->protons;        break;
        case neutron:       esc->particles = &st->neutrons;       break;

        case positive_pion: esc->particles = &st->positive_pions; break;
        case neutral_pion:  esc->particles = &st->neutral_pions;  break;
        case negative_pion: esc->particles = &st->negative_pions; break;

        case positive_left_muon:  esc->particles = &st->positive_left_muons;  break;
        case positive_right_muon: esc->particles = &st->positive_right_muons; break;
        case negative_left_muon:  esc->particles = &st->negative_left_muons;  break;
        case negative_right_muon: esc->particles = &st->negative_right_muons; break;

        case electron_neutrino:     esc->particles = &st->electron_neutrinos;       break;
        case electron_antineutrino: esc->particles = &st->electron_antineutrinos;   break;
        case muon_neutrino:         esc->particles = &st->muon_neutrinos;           break;
        case muon_antineutrino:     esc->particles = &st->muon_antineutrinos;       break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    MEM_PREPARE(esc->losses,        esc->particles->size);
    MEM_PREPARE(esc->escape_factor, esc->particles->size);

    esc->t = t;
    esc->dt = st->dt;

    unsigned int i;
    for(i = 0; i < esc->particles->size; i++)
        esc->escape_factor[i] = expm1(-esc->dt / t) / esc->dt;

    esc->escape_function = free_escape;
}

void update_escape(struct state_t *st, escape_t *esc, double t)
{
    esc->t = t;
    esc->dt = st->dt;

    unsigned int i;
    for(i = 0; i < esc->particles->size; i++)
        esc->escape_factor[i] = expm1(-esc->dt / t) / esc->dt;
}

void free_escape(escape_t *esc)
{
    unsigned int i;
    population_t *p = esc->particles;

    for(i = 0; i < p->size; i++)
        esc->losses[i] = p->population[i] * esc->escape_factor[i];
}
