#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <assert.h>

#include "constants.h"
#include "decay_and_escape.h"
#include "state.h"

static void free_decay_and_escape(decay_and_escape_t *dae);

#define MEM_PREPARE(X,S) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(double) * (S));\
        memset(X, '\0', sizeof(double) * (S));\
        \
    } while(0);

void init_decay_and_escape(struct state_t *st, decay_and_escape_t *dae, enum particle_type pt, double escape_lifetime)
{
    switch(pt)
    {
        case electron:      dae->particles = &st->electrons;      break;
        case photon:        dae->particles = &st->photons;        break;
        case proton:        dae->particles = &st->protons;        break;
        case neutron:       dae->particles = &st->neutrons;       break;

        case positive_pion: dae->particles = &st->positive_pions; break;
        case neutral_pion:  dae->particles = &st->neutral_pions;  break;
        case negative_pion: dae->particles = &st->negative_pions; break;

        case positive_left_muon:  dae->particles = &st->positive_left_muons;  break;
        case positive_right_muon: dae->particles = &st->positive_right_muons; break;
        case negative_left_muon:  dae->particles = &st->negative_left_muons;  break;
        case negative_right_muon: dae->particles = &st->negative_right_muons; break;

        case electron_neutrino:     dae->particles = &st->electron_neutrinos;       break;
        case electron_antineutrino: dae->particles = &st->electron_antineutrinos;   break;
        case muon_neutrino:         dae->particles = &st->muon_neutrinos;           break;
        case muon_antineutrino:     dae->particles = &st->muon_antineutrinos;       break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    switch(pt)
    {
        case neutron:
            dae->decay_lifetime = NEUTRON_LIFETIME;
            break;

        case neutral_pion:
            dae->decay_lifetime = NEUTRAL_PION_LIFETIME;
            break;

        case positive_pion:
        case negative_pion:
            dae->decay_lifetime = CHARGED_PION_LIFETIME;
            break;

        case positive_left_muon:
        case positive_right_muon:
        case negative_left_muon:
        case negative_right_muon:
            dae->decay_lifetime = MUON_LIFETIME;
            break;

        case electron:
        case photon:
        case proton:
        case electron_neutrino:
        case electron_antineutrino:
        case muon_neutrino:
        case muon_antineutrino:
            dae->decay_lifetime = INFINITY;
            break;

        // TODO: ADD FAIL MESSAGES
        default: assert(0); break;
    }

    MEM_PREPARE(dae->t,             dae->particles->size);
    MEM_PREPARE(dae->losses,        dae->particles->size);
    MEM_PREPARE(dae->losses_factor, dae->particles->size);

    dae->losses_function = free_decay_and_escape;

    update_decay_and_escape(st, dae, escape_lifetime);
}

void update_decay_and_escape(struct state_t *st, decay_and_escape_t *dae, double escape_lifetime)
{
    dae->escape_lifetime = escape_lifetime;
    dae->dt = st->dt;

    unsigned int i;
    for(i = 0; i < dae->particles->size; i++)
    {
        double t;

        if(isinf(dae->decay_lifetime))
            t = 1 / dae->escape_lifetime;
        else
            t = (dae->particles->energy[i] * dae->decay_lifetime + dae->escape_lifetime) / 
                (dae->particles->energy[i] * dae->decay_lifetime * dae->escape_lifetime);

        dae->t[i] = t;
        dae->losses_factor[i] = expm1(-dae->dt * t) / dae->dt;
    }
}

static void free_decay_and_escape(decay_and_escape_t *dae)
{
    unsigned int i;
    population_t *p = dae->particles;

    for(i = 0; i < p->size; i++)
        dae->losses[i] = p->population[i] * dae->losses_factor[i];
}
