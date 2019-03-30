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

typedef struct _distribution_metadata_t
{
    enum distribution_type dt;
    double min, max;

    union {
        double t;   // maxwell juttner && hybrid
        struct      // power-law && power-law with exp cuttoff
        {
            double p, e;
        };
        struct      // broken power-law && connected power-law
        {
            double gc, p1, p2;
        };
    };
    
} distribution_metadata_t;

void generate_distribution(double *population, double *energy, distribution_metadata_t *dm, unsigned int size);

void generate_maxwell_juttner(double *population, double *energy, double theta, unsigned int size);
void generate_power_law(double *population, double *energy, double p, unsigned int size);
void generate_broken_power_law(double *population, double *energy, double gc, double p1, double p2, unsigned int size);
void generate_power_law_with_exponential_cutoff(double *population, double *energy, double p, double e, unsigned int size);
void generate_hybrid(double *population, double *energy, double theta, unsigned int size);
void generate_connected_power_law(double *population, double *energy, double gc, double p1, double p2, unsigned int size);

double distribution_average(distribution_metadata_t *dm);

double maxwell_juttner_average(double theta);
double power_law_average(double energy_min, double energy_max, double p);
double broken_power_law_average(double energy_min, double energy_max,
        double energy_break, double p1, double p2);

#endif /* end of include guard: DISTRIBUTION_H */
