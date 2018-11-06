#ifndef ACCELERATION_H
#define ACCELERATION_H
#pragma once

#include "population.h"

enum acceleration_type
{
    linear,
    quadratic
};

struct state_t;
typedef struct acceleration_t
{
    double *temp1, *temp2;  // Needed for the derivatives!
    double *gains;
    population_t *particles;
    double t;

    double dt;      // Needed for implicit and log-log methods!

    void (*acceleration_function)(struct acceleration_t *acc);
} acceleration_t;

void init_acceleration(struct state_t *st, acceleration_t *acc, enum particle_type pt, double t);

void update_acceleration(struct state_t *st, acceleration_t *acc, double t);

void linear_acceleration(acceleration_t *acc);
void quadratic_acceleration(acceleration_t *acc);

#endif /* end of include guard: ACCELERATION_H */
