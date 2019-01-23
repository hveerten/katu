#ifndef DECAY_AND_ESCAPE_H
#define DECAY_AND_ESCAPE_H
#pragma once

#include "population.h"

struct state_t;
typedef struct decay_and_escape_t
{
    double *t;
    double *losses;
    double *losses_factor;
    population_t *particles;
    double decay_lifetime;
    double escape_lifetime;
    double dt;      // Needed for implicit and log-log methods!

    void (*losses_function)(struct decay_and_escape_t *dea);
} decay_and_escape_t;

void init_decay_and_escape(struct state_t *st, decay_and_escape_t *dea, enum particle_type pt, double escape_lifetime);

void update_decay_and_escape(struct state_t *st, decay_and_escape_t *dea, double escape_lifetime);

#endif /* end of include guard: DECAY_AND_ESCAPE_H */
