#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "external_libs/toml.h"

#define TOML_READ_DOUBLE(X,Y,Z,A,D)                             \
    raw = toml_raw_in((X), Y);                                  \
    if(raw != 0)                                                \
    {                                                           \
        if(toml_rtod(raw, &(Z)))                             \
        {                                                       \
            fprintf(stderr,"Error reading the "A" ("Y")\n");    \
            exit(1);                                            \
        }                                                       \
    }                                                           \
    else                                                        \
    {                                                           \
        (Z) = D;                                                \
    };


static void config_read_distribution_type(
        toml_table_t *part_table, char *part_txt,
        distribution_metadata_t *dm)
{
    const char *raw;

    raw = toml_raw_in(part_table, "distribution_type");
    if(raw != 0)
    {
        char *dist_type;

        if(toml_rtos(raw, &dist_type))
        {
            fprintf(stderr,"Error reading the distribution type for %s\n", part_txt);
            exit(1);
        }

        if(strcmp(dist_type, "maxwell_juttner") == 0)                    dm->dt = maxwell_juttner;
        if(strcmp(dist_type, "black_body") == 0)                         dm->dt = black_body;
        if(strcmp(dist_type, "power_law") == 0)                          dm->dt = power_law;
        if(strcmp(dist_type, "broken_power_law") == 0)                   dm->dt = broken_power_law;
        if(strcmp(dist_type, "power_law_with_exponential_cutoff") == 0)  dm->dt = power_law_with_exponential_cutoff;
        if(strcmp(dist_type, "hybrid") == 0)                             dm->dt = hybrid;
        if(strcmp(dist_type, "connected_power_law") == 0)                dm->dt = connected_power_law;

        free(dist_type);
    }
    else
        dm->dt = broken_power_law;

    if(strcmp(part_txt, "protons") == 0)
    {
        TOML_READ_DOUBLE(part_table, "gamma_min", dm->min, "proton minimum gamma", 1e1);
        TOML_READ_DOUBLE(part_table, "gamma_max", dm->max, "proton maximum gamma", 1e5);
    }
    if(strcmp(part_txt, "electrons") == 0)
    {
        TOML_READ_DOUBLE(part_table, "gamma_min", dm->min, "electron minimum gamma", 1e1);
        TOML_READ_DOUBLE(part_table, "gamma_max", dm->max, "electron maximum gamma", 1e5);
    }
    if(strcmp(part_txt, "photons") == 0)
    {
        TOML_READ_DOUBLE(part_table, "epsilon_min", dm->min, "photon minimum epsilon", 1e-12);
        TOML_READ_DOUBLE(part_table, "epsilon_max", dm->max, "photon maximum epsilon", 1e5);
    }

    switch(dm->dt)
    {
        case maxwell_juttner:
            TOML_READ_DOUBLE(part_table, "temperature", dm->t, "temperature", 0.5)
            break;

        case black_body:
            TOML_READ_DOUBLE(part_table, "temperature", dm->t, "temperature", 0.5)
            break;

        case power_law:
            TOML_READ_DOUBLE(part_table, "slope", dm->p, "slope", 2.3)
            break;

        case broken_power_law:
            TOML_READ_DOUBLE(part_table, "break_point",  dm->gc, "break point", 1e4)
            TOML_READ_DOUBLE(part_table, "first_slope",  dm->p1, "first slope", 2.15)
            TOML_READ_DOUBLE(part_table, "second_slope", dm->p2, "second slope", 4.0)
            break;

        case power_law_with_exponential_cutoff:
            TOML_READ_DOUBLE(part_table, "slope",       dm->p, "slope", -2.3)
            TOML_READ_DOUBLE(part_table, "break_point", dm->e, "break point", 1e4)
            break;

        case hybrid:
            TOML_READ_DOUBLE(part_table, "temperature", dm->t, "temperature", 0.5)
            TOML_READ_DOUBLE(part_table, "slope",       dm->p, "slope", -2.3)
            break;

        case connected_power_law:
            TOML_READ_DOUBLE(part_table, "connection_point", dm->gc, "connection point", 1e4)
            TOML_READ_DOUBLE(part_table, "first_slope",      dm->p1, "first slope", -2.3)
            TOML_READ_DOUBLE(part_table, "second_slope",     dm->p2, "second slope", 4.0)
            break;

        default:
            assert(0);
    }
}

void config_read_file(config_t *cfg, char *filename)
{
    toml_table_t *conf;
    toml_table_t *general_table;
    toml_table_t *volume_table;
    toml_table_t *photon_table;
    toml_table_t *proton_table;
    toml_table_t *electron_table;

    toml_table_t *external_injection_table;
    toml_table_t *external_injection_photon_table;
    toml_table_t *external_injection_proton_table;
    toml_table_t *external_injection_electron_table;

    FILE *fin = fopen(filename, "r");
    const char *raw;
    char errbuf[256];

    if(fin == NULL)
    {
        fprintf(stderr,"Non existant config file: %s\n", filename);
        exit(1);
    }

    conf = toml_parse_file(fin, errbuf, sizeof(errbuf));
    fclose(fin);

    if(conf == 0)
    {
        fprintf(stderr,"Error parsing config file: %s\n", errbuf);
        exit(1);
    }

    general_table = toml_table_in(conf, "general");
    if(general_table != 0)
    {
        TOML_READ_DOUBLE(general_table, "density",        cfg->density,        "density",                   1e-22);
        TOML_READ_DOUBLE(general_table, "magnetic_field", cfg->magnetic_field, "magnetic field",            0.1);
        TOML_READ_DOUBLE(general_table, "dt",             cfg->dt,             "time step",                 0.1);
        TOML_READ_DOUBLE(general_table, "dt_max",         cfg->dt_max,         "maximum time step",         1e3);
        TOML_READ_DOUBLE(general_table, "t_max",          cfg->t_max,          "maximum time",              1e7);
        TOML_READ_DOUBLE(general_table, "t_acc",          cfg->t_acc,          "acceleration timescale",    1e100);
        TOML_READ_DOUBLE(general_table, "eta",            cfg->eta,            "proton to electron_ratio",  1.0);
        TOML_READ_DOUBLE(general_table, "cfe_ratio",      cfg->cfe_ratio,      "charged free escape ratio", 1.0);
        TOML_READ_DOUBLE(general_table, "tol",            cfg->tol,            "tolerance",                 1e-8);
    }

    volume_table = toml_table_in(conf, "volume");
    if(volume_table != 0)
    {
        TOML_READ_DOUBLE(volume_table, "R", cfg->R, "radius", 1e16);
        TOML_READ_DOUBLE(volume_table, "h", cfg->h, "height", 1e16);

        char *shape;
        raw = toml_raw_in(volume_table, "shape");

        if(raw != 0)
        {
            if(toml_rtos(raw, &shape))
            {
                fprintf(stderr,"Error reading the shape of the volume\n");
                exit(1);
            }

            if(strcmp(shape, "sphere") == 0) cfg->v = sphere;
            if(strcmp(shape, "shell") == 0)  cfg->v = shell;

            free(shape);
        }
        else
            cfg->v = sphere;
    }

    electron_table = toml_table_in(conf, "electrons");
    if(electron_table != 0)
    {
        TOML_READ_DOUBLE(electron_table, "size",      cfg->electron_size,      "number of points for electrons", 256);

        config_read_distribution_type(electron_table, "electrons", &cfg->electron_distribution);
    }

    proton_table = toml_table_in(conf, "protons");
    if(proton_table != 0)
    {
        TOML_READ_DOUBLE(proton_table, "size",      cfg->proton_size,      "number of points for protons", 256);

        config_read_distribution_type(proton_table, "protons", &cfg->proton_distribution);
    }

    photon_table = toml_table_in(conf, "photons");
    if(photon_table != 0)
    {
        TOML_READ_DOUBLE(photon_table, "size",        cfg->photon_size,        "number of points for photons", 160);

        config_read_distribution_type(photon_table, "photons", &cfg->photon_distribution);
    }

    external_injection_table  = toml_table_in(conf, "external_injection");
    if(external_injection_table != 0)
    {
        TOML_READ_DOUBLE(external_injection_table, "luminosity", cfg->ei.luminosity, "luminosity",               0.0);
        TOML_READ_DOUBLE(external_injection_table, "eta",        cfg->ei.eta,        "proton to electron_ratio", 1.0);

        external_injection_electron_table = toml_table_in(external_injection_table, "electrons");
        if(external_injection_electron_table != 0)
        {
            config_read_distribution_type(external_injection_electron_table, "electrons", &cfg->ei.electron_distribution);
        }

        external_injection_proton_table = toml_table_in(external_injection_table, "protons");
        if(external_injection_proton_table != 0)
        {
            config_read_distribution_type(external_injection_proton_table, "protons", &cfg->ei.proton_distribution);
        }

        external_injection_photon_table = toml_table_in(external_injection_table, "photons");
        if(external_injection_photon_table != 0)
        {
            TOML_READ_DOUBLE(external_injection_photon_table, "luminosity", cfg->ei.photon_luminosity, "luminosity", 0.0);

            config_read_distribution_type(external_injection_photon_table, "photons", &cfg->ei.photon_distribution);
        }
    }

    toml_free(conf);
}
