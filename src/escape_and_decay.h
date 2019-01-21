#ifndef ESCAPE_AND_DECAY_H
#define ESCAPE_AND_DECAY_H
#pragma once

#include "population.h"

enum escape_and_decay_type
{
    free_esc_and_decay,
};

struct state_t;
typedef struct escape_and_decay_t
{
    double *losses;
    double *losses_factor;
    population_t *particles;
    double escape_lifetime;
    double decay_lifetime;
    double dt;      // Needed for implicit and log-log methods!

    void (*losses_function)(struct escape_and_decay_t *ead);
} escape_and_decay_t;

void init_escape_and_decay(struct state_t *st, escape_and_decay_t *ead, enum particle_type pt, double escape_lifetime);

void update_escape_and_decay(struct state_t *st, escape_and_decay_t *ead, double escape_lifetime);

void free_escape_and_decay(escape_and_decay_t *ead);

#endif /* end of include guard: ESCAPE_H */
