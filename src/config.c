#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "../external_libs/toml.h"

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
        enum distribution_type *part_dist_type, double *part_params)
{
    const char *raw;
    char *dist_type;

    raw = toml_raw_in(part_table, "distribution_type");
    if(raw != 0)
    {
        if(toml_rtos(raw, &dist_type))
        {
            fprintf(stderr,"Error reading the distribution type for %s\n", part_txt);
            exit(1);
        }

        if(strcmp(dist_type, "maxwell_juttner") == 0)                    *part_dist_type = maxwell_juttner;
        if(strcmp(dist_type, "power_law") == 0)                          *part_dist_type = power_law;
        if(strcmp(dist_type, "broken_power_law") == 0)                   *part_dist_type = broken_power_law;
        if(strcmp(dist_type, "power_law_with_exponential_cutoff") == 0)  *part_dist_type = power_law_with_exponential_cutoff;
        if(strcmp(dist_type, "hybrid") == 0)                             *part_dist_type = hybrid;
        if(strcmp(dist_type, "connected_power_law") == 0)                *part_dist_type = connected_power_law;
    }
    else
        *part_dist_type = broken_power_law;

    switch(*part_dist_type)
    {
        case maxwell_juttner:
            TOML_READ_DOUBLE(part_table, "temperature", part_params[0], "temperature", 0.5)
            break;

        case power_law:
            TOML_READ_DOUBLE(part_table, "slope", part_params[0], "slope", 2.3)
            break;

        case broken_power_law:
            TOML_READ_DOUBLE(part_table, "break_point",  part_params[0], "break point", 1e4)
            TOML_READ_DOUBLE(part_table, "first_slope",  part_params[1], "first slope", 2.15)
            TOML_READ_DOUBLE(part_table, "second_slope", part_params[2], "second slope", 4.0)
            break;

        case power_law_with_exponential_cutoff:
            TOML_READ_DOUBLE(part_table, "slope",       part_params[0], "slope", -2.3)
            TOML_READ_DOUBLE(part_table, "break_point", part_params[1], "break point", 1e4)
            break;

        case hybrid:
            TOML_READ_DOUBLE(part_table, "temperature", part_params[0], "temperature", 0.5)
            TOML_READ_DOUBLE(part_table, "slope",       part_params[1], "slope", -2.3)
            break;

        case connected_power_law:
            TOML_READ_DOUBLE(part_table, "connection_point", part_params[0], "connection point", 1e4)
            TOML_READ_DOUBLE(part_table, "first_slope",      part_params[1], "first slope", -2.3)
            TOML_READ_DOUBLE(part_table, "second_slope",     part_params[2], "second slope", 4.0)
            break;

        default:
            assert(0);
    }

    free(dist_type);
}

void config_read_file(config_t *cfg, char *filename)
{
    toml_table_t *conf;
    toml_table_t *general_table;
    toml_table_t *photon_table;
    toml_table_t *proton_table;
    toml_table_t *electron_table;

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
        TOML_READ_DOUBLE(general_table, "density",        cfg->density,        "density",                  1e-22);
        TOML_READ_DOUBLE(general_table, "magnetic_field", cfg->magnetic_field, "magnetic field",           0.1);
        TOML_READ_DOUBLE(general_table, "dt",             cfg->dt,             "time step",                0.1);
        TOML_READ_DOUBLE(general_table, "dt_max",         cfg->dt_max,         "maximum time step",        1e3);
        TOML_READ_DOUBLE(general_table, "t_max",          cfg->t_max,          "maximum time",             1e7);
        TOML_READ_DOUBLE(general_table, "R",              cfg->R,              "radius",                   1e16);
        TOML_READ_DOUBLE(general_table, "eta",            cfg->eta,            "proton to electron_ratio", 1.0);
    }

    electron_table = toml_table_in(conf, "electrons");
    if(electron_table != 0)
    {
        TOML_READ_DOUBLE(electron_table, "gamma_min", cfg->electron_gamma_min, "electron minimum gamma", 1e1);
        TOML_READ_DOUBLE(electron_table, "gamma_max", cfg->electron_gamma_max, "electron maximum gamma", 1e5);
        TOML_READ_DOUBLE(electron_table, "size",      cfg->electron_size,      "number of points for electrons", 256);

        config_read_distribution_type(electron_table, "electrons", &cfg->electron_distribution_type, cfg->electron_params);
    }

    proton_table = toml_table_in(conf, "protons");
    if(proton_table != 0)
    {
        TOML_READ_DOUBLE(proton_table, "gamma_min", cfg->proton_gamma_min, "proton minimum gamma", 1e1);
        TOML_READ_DOUBLE(proton_table, "gamma_max", cfg->proton_gamma_max, "proton maximum gamma", 1e5);
        TOML_READ_DOUBLE(proton_table, "size",      cfg->proton_size,      "number of points for protons", 256);

        config_read_distribution_type(proton_table, "protons", &cfg->proton_distribution_type, cfg->proton_params);
    }

    photon_table = toml_table_in(conf, "photons");
    if(photon_table != 0)
    {
        TOML_READ_DOUBLE(photon_table, "epsilon_min", cfg->photon_epsilon_min, "photon minimum epsilon", 1e-12);
        TOML_READ_DOUBLE(photon_table, "epsilon_max", cfg->photon_epsilon_max, "photon maximum epsilon", 1e5);
        TOML_READ_DOUBLE(photon_table, "size",        cfg->photon_size,        "number of points for photons", 160);
    }

    toml_free(conf);
}
