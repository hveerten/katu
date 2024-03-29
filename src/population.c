#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "constants.h"
#include "population.h"

void init_population(population_t *p, enum particle_type pt,
        double min_energy, double max_energy, unsigned int size)
{
/*
 *    switch(pt)
 *    {
 *        case photon:        p = &st->photons;        break;
 *        case electron:      p = &st->electrons;      break;
 *        case positron:      p = &st->positrons;      break;
 *        case proton:        p = &st->protons;        break;
 *        case neutron:       p = &st->neutrons;       break;
 *        case neutral_pion:  p = &st->neutral_pions;  break;
 *        case positive_pion: p = &st->positive_pions; break;
 *        case negative_pion: p = &st->negative_pions; break;
 *
 *        default: p = NULL; break;
 *    }
 */

    switch(pt)
    {
        case photon:        p->mass = ELECTRON_MASS; break;
        case electron:      p->mass = ELECTRON_MASS; break;
        case positron:      p->mass = ELECTRON_MASS; break;
        case proton:        p->mass = PROTON_MASS;   break;
        case neutron:       p->mass = PROTON_MASS;   break; // <----------
        case neutral_pion:  p->mass = PION_MASS;     break;
        case positive_pion: p->mass = PION_MASS;     break; // <----------
        case negative_pion: p->mass = PION_MASS;     break; // <----------

        case positive_left_muon:  p->mass = MUON_MASS; break;
        case positive_right_muon: p->mass = MUON_MASS; break;
        case negative_left_muon:  p->mass = MUON_MASS; break;
        case negative_right_muon: p->mass = MUON_MASS; break;

        // By default, neutrinos use electron mass
        case electron_neutrino:     p->mass = ELECTRON_MASS; break;
        case electron_antineutrino: p->mass = ELECTRON_MASS; break;
        case muon_neutrino:         p->mass = ELECTRON_MASS; break;
        case muon_antineutrino:     p->mass = ELECTRON_MASS; break;

        default: p->mass = ELECTRON_MASS; break;
    }

    p->size = size;
    posix_memalign((void **) &p->energy,               32, sizeof(double) * size);
    posix_memalign((void **) &p->log_energy,           32, sizeof(double) * size);
    posix_memalign((void **) &p->population,           32, sizeof(double) * size);
    posix_memalign((void **) &p->log_population,       32, sizeof(double) * size);
    posix_memalign((void **) &p->tentative_population, 32, sizeof(double) * size);

    double logde = (log(max_energy) - log(min_energy)) / (size - 1);

    for(unsigned i = 0; i < size; i++)
        p->log_energy[i] = log(min_energy) + logde * i;
    for(unsigned i = 0; i < size; i++)
        p->energy[i] = exp(p->log_energy[i]);
    for(unsigned i = 0; i < size; i++)
        p->log_population[i] = log(DBL_MIN);
    for(unsigned i = 0; i < size; i++)
        p->population[i] = DBL_MIN;
    for(unsigned i = 0; i < size; i++)
        p->tentative_population[i] = DBL_MIN;
}

void free_population(population_t *p)
{
    free(p->energy);
    free(p->log_energy);
    free(p->population);
    free(p->log_population);
    free(p->tentative_population);
}

void generate_population(population_t *p, distribution_metadata_t *dm)
{
    unsigned int i;

    generate_distribution(p->population, p->energy, dm, p->size);

    // Some distributions (MJ for example) can leave high energy particles
    // at 0, so to avoid infs and nans down the road, replace 0 by DBL_MIN
    for(i = 0; i < p->size; i++)
        p->population[i] = p->population[i] ? p->population[i] : DBL_MIN;

    for(i = 0; i < p->size; i++)
        p->log_population[i] = log(p->population[i]);
}

double calculate_population(population_t *p)
{
    unsigned int i;

    double pop = 0;
    double dlng = p->log_energy[1] - p->log_energy[0];

    pop += p->population[0]           * p->energy[0];
    pop += p->population[p->size - 1] * p->energy[p->size - 1];

    for(i = 1; i < p->size - 1; i++)
        pop += 2 * p->population[i] * p->energy[i];
    pop *= dlng / 2;

    return pop;
}

double calculate_energy(population_t *p)
{
    unsigned int i;

    double energy = 0;
    double dlng = p->log_energy[1] - p->log_energy[0];

    energy += p->population[0]           * p->energy[0]           * p->energy[0];
    energy += p->population[p->size - 1] * p->energy[p->size - 1] * p->energy[p->size - 1];

    for(i = 1; i < p->size - 1; i++)
        energy += 2 * p->population[i] * p->energy[i] * p->energy[i];
    energy *= dlng / 2;

    return energy;
}
