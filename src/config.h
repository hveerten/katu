#ifndef CONFIG_H
#define CONFIG_H
#pragma once

#include "distribution.h"

enum volume_type
{
    sphere,
    shell
};

typedef struct
{
    double density;
    double magnetic_field;
    double dt;
    double dt_max;
    double t_max;
    double eta;

    enum volume_type v;
    double R;
    double h;

    distribution_metadata_t electron_distribution;
    double electron_gamma_min;
    double electron_gamma_max;
    double electron_size;
    enum distribution_type electron_distribution_type;
    double electron_params[3];

    distribution_metadata_t proton_distribution;
    double proton_gamma_min;
    double proton_gamma_max;
    double proton_size;
    enum distribution_type proton_distribution_type;
    double proton_params[3];

    distribution_metadata_t photon_distribution;
    double photon_epsilon_min;
    double photon_epsilon_max;
    double photon_size;

    struct
    {
        double luminosity;
        double eta;

        distribution_metadata_t electron_distribution;
        enum distribution_type electron_distribution_type;
        double electron_params[3];

        distribution_metadata_t proton_distribution;
        enum distribution_type proton_distribution_type;
        double proton_params[3];

        double photon_luminosity;
        distribution_metadata_t photon_distribution;
        enum distribution_type photon_distribution_type;
        double photon_params[3];
    } ei;

} config_t;

void config_read_file(config_t *cfg, char *filename);

#endif /* end of include guard: CONFIG_H */
