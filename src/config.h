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
    double cfe_ratio;
    double tol;

    enum volume_type v;
    double R;
    double h;

    distribution_metadata_t electron_distribution;
    double electron_size;

    distribution_metadata_t proton_distribution;
    double proton_size;

    distribution_metadata_t photon_distribution;
    double photon_size;

    struct
    {
        double luminosity;
        double eta;

        distribution_metadata_t electron_distribution;
        distribution_metadata_t proton_distribution;

        double photon_luminosity;
        distribution_metadata_t photon_distribution;
    } ei;

} config_t;

void config_read_file(config_t *cfg, char *filename);

#endif /* end of include guard: CONFIG_H */
