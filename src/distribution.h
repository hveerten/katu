#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#pragma once

enum distribution_type
{
    maxwell_juttner,
    power_law,
    broken_power_law,
    power_law_with_exponential_cutoff,
    hybrid,
    connected_power_law
};

void generate_maxwell_juttner(double *population, double *energy, double theta, unsigned int size);
void generate_power_law(double *population, double *energy, double p, unsigned int size);
void generate_broken_power_law(double *population, double *energy, double gc, double p1, double p2, unsigned int size);
void generate_power_law_with_exponential_cutoff(double *population, double *energy, double p, double e, unsigned int size);
void generate_hybrid(double *population, double *energy, double theta, unsigned int size);
void generate_connected_power_law(double *population, double *energy, double gc, double p1, double p2, unsigned int size);

#endif /* end of include guard: DISTRIBUTION_H */
