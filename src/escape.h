#ifndef ESCAPE_H
#define ESCAPE_H
#pragma once

#include "population.h"

enum escape_type
{
    free_esc,
};

struct state_t;
typedef struct escape_t
{
    double *losses;
    double *escape_factor;
    population_t *particles;
    double t;
    double dt;      // Needed for implicit and log-log methods!

    void (*escape_function)(struct escape_t *esc);
} escape_t;

void init_escape(struct state_t *st, escape_t *esc, enum particle_type pt, double t);

void update_escape(struct state_t *st, escape_t *esc, double t);

void free_escape(escape_t *esc);

#endif /* end of include guard: ESCAPE_H */
