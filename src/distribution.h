#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#pragma once

enum distribution_type
{
    maxwell_juttner,
    hybrid,
    power_law,
    power_law_with_exponential_cutoff,
    broken_power_law,
    connected_power_law,
    black_body
};

typedef struct _distribution_metadata_t
{
    enum distribution_type dt;
    double min, max;

    double t;       // maxwell juttner, hybrid and black body
    union {
        double gc;  // hybrid, broken power-law and connected power-law
        double e;   // power-law with exponential cutoff
    };
    union {
        double p;   // power-law and power-law with exponential cutoff
        double p1;  // broken power law and connected power-law
    };
    double p2;      // broken power law and connected power-law

} distribution_metadata_t;

void generate_distribution(double *population, double *energy, distribution_metadata_t *dm, unsigned int size);

double distribution_norm(distribution_metadata_t *dm);
double distribution_average(distribution_metadata_t *dm);

#endif /* end of include guard: DISTRIBUTION_H */
