#ifndef STATE_STEP_H
#define STATE_STEP_H
#pragma once

#include <stdbool.h>
#include "state.h"

void step(state_t *st);
void step_tentative(state_t *st, bool try_new_step);

void step_heun(state_t *st);
void step_heun_tentative(state_t *st, bool try_new_step);

#endif /* end of include guard: STATE_STEP_H */
