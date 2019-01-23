#ifndef DECAY_H
#define DECAY_H
#pragma once

#include "population.h"

struct state_t;
typedef struct decay_t
{
    double *losses;
    double *decay_factor;
    population_t *particles;
    double t;
    double dt;      // Needed for implicit and log-log methods!

    void (*decay_function)(struct decay_t *d);
} decay_t;

void init_decay(struct state_t *st, decay_t *d, enum particle_type pt);

void update_decay(struct state_t *st, decay_t *d);

#endif /* end of include guard: DECAY_H */
