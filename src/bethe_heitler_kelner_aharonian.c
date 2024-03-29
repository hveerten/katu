#include "bethe_heitler_kelner_aharonian.h"
#include "constants.h"

#include "utils.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <stdio.h>

static double W1(double k, double gm, double d)
{
    double gp = k - gm;

    double pm2 = gm * gm - 1;
    double pp2 = gp * gp - 1;

    double T2 = pp2 + 2 * k * d;


    double aux0 = (d*d - 2*d*gm + 1) / pm2;
    double aux1 = (2 * gm*gm + 1) / (pm2 * d*d*d);

    double aux2 = 5 / d;
    double aux3 = 2 * (4 - gm*gp) / (pm2 * d);
    double aux4 = (pm2 - k*k) / (T2 * d);
    double aux5 = 2 * gp / pm2;

    return sqrt(pp2) * (4 * aux0 * aux1 + aux2 + aux3 + aux4 + aux5);
}

static double W2(double k, double gm, double d)
{
    double gp = k - gm;

    double pm2 = gm * gm - 1;
    double pp2 = gp * gp - 1;

    double pm = sqrt(pm2);
    double pp = sqrt(pp2);

    double Y = (gm > 1e6 && gp > 1e6) ?
                    log(4 * gm * gp * k / (2 *k*k + 1))
            :
                    log((gm * gp + pm * pp + 1) / k);
    Y *= 2 / pm2;


    double aux0 = (d*d - 2*d*gm + 1) / pm2;
    double aux1 = (3*k + pm2*gp) / (d*d*d);

    double aux2 = (2*gm*gm*(gm*gm + gp*gp) - 7*gm*gm - 3*gm*gp - gp*gp + 1) / d;
    double aux3 = k * (pm2 - gm*gp);

    return Y / pm * (-2*gm*aux0 * aux1 + aux2 + aux3);
}

static double W3(double k, double gm, double d)
{
    double gp = k - gm;

    double pm2 = gm * gm - 1;
    double pp2 = gp * gp - 1;

    double pp = sqrt(pp2);

    double T2 = pp2 + 2 * k * d;
    double T  = sqrt(T2);

    double dp = log1p(pp / (k*d) * (pp + T));


    double aux0 = 2 / d;
    double aux1 = -3*k;
    double aux2 = -k * (pm2 - k*k) / T2;

    return -dp / T * (aux0 + aux1 + aux2);
}

static double W4(double k, double gm, double d)
{
    double gp = k - gm;

    double bp = gp < 1e6 ?
                sqrt(1 - 1/(gp*gp))
            :
                1 - 1 / (2 * gp*gp);

    double aux1 = bp != 1 ?
                atanh(bp)
            :
                log(gp) + 1 / (8*gp*gp);


    return -4*aux1;
}

struct W_params { double k; double d; };
static double _W1(double x, void *params)
{
    double xx = exp(x);
    struct W_params p = *(struct W_params *)params;
    return xx * W1(p.k, xx, p.d);
}
static double _W2(double x, void *params)
{
    double xx = exp(x);
    struct W_params p = *(struct W_params *)params;
    return xx * W2(p.k, xx, p.d);
}
static double _W3(double x, void *params)
{
    double xx = exp(x);
    struct W_params p = *(struct W_params *)params;
    return xx * W3(p.k, xx, p.d);
}
static double _W4(double x, void *params)
{
    double xx = exp(x);
    struct W_params p = *(struct W_params *)params;
    return xx * W4(p.k, xx, p.d);
}


static double G(double k, double d)
{
    double e;
    struct W_params p;
    p.k = k;
    p.d = d;

    double lo = 0.5 * (d + 1/d);
    double hi = k - 1;

    if(lo == 1) lo = 1.0001;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);
    gsl_function F1 = (gsl_function) { .function = &_W1, .params = &p };
    gsl_function F2 = (gsl_function) { .function = &_W2, .params = &p };
    gsl_function F3 = (gsl_function) { .function = &_W3, .params = &p };
    gsl_function F4 = (gsl_function) { .function = &_W4, .params = &p };

    if(lo >= hi) return 0;

    if(lo < hi/2)
    {
        double G11 = 0, G12 = 0;
        double G21 = 0, G22 = 0;
        double G31 = 0, G32 = 0;
        double G41 = 0, G42 = 0;

        gsl_integration_qags(&F1, log(lo),   log(hi/2), 0, 1e-4, 256, w, &G11, &e);
        gsl_integration_qags(&F2, log(lo),   log(hi/2), 0, 1e-4, 256, w, &G21, &e);
        gsl_integration_qags(&F3, log(lo),   log(hi/2), 0, 1e-4, 256, w, &G31, &e);
        gsl_integration_qags(&F4, log(lo),   log(hi/2), 0, 1e-4, 256, w, &G41, &e);

        gsl_integration_qags(&F1, log(hi/2), log(hi),   0, 1e-4, 256, w, &G12, &e);
        gsl_integration_qags(&F2, log(hi/2), log(hi),   0, 1e-4, 256, w, &G22, &e);
        gsl_integration_qags(&F3, log(hi/2), log(hi),   0, 1e-4, 256, w, &G32, &e);
        gsl_integration_qags(&F4, log(hi/2), log(hi),   0, 1e-4, 256, w, &G42, &e);

        gsl_integration_workspace_free(w);

        double r = G11 + G12 + G21 + G22 + G31 + G32 + G41 + G42;

        if(!isnormal(r) || r < 0)
        {
            fprintf(stderr,"%lg\t%lg\n", k, d);
            fprintf(stderr,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    G11, G12, G21, G22, G31, G32, G41, G42, r);
        }

        return r;
    }
    else
    {
        double G1 = 0;
        double G2 = 0;
        double G3 = 0;
        double G4 = 0;

        gsl_integration_qags(&F1, log(lo), log(hi), 0, 1e-4, 256, w, &G1, &e);
        gsl_integration_qags(&F2, log(lo), log(hi), 0, 1e-4, 256, w, &G2, &e);
        gsl_integration_qags(&F3, log(lo), log(hi), 0, 1e-4, 256, w, &G3, &e);
        gsl_integration_qags(&F4, log(lo), log(hi), 0, 1e-4, 256, w, &G4, &e);

        gsl_integration_workspace_free(w);

        double r = G1 + G2 + G3 + G4;

        if(!isnormal(r) || r < 0)
        {
            fprintf(stderr,"%lg %lg\n", k, d);
            fprintf(stderr,"%lg\t%lg\t%lg\t%lg\t%lg\n",
                    G1, G2, G3, G4, r);
        }

        return G1 + G2 + G3 + G4;
    }
}

struct G_params { double d; };
static double _G(double x, void *params)
{
    double xx = exp(x);
    struct G_params p = *(struct G_params *)params;
    return G(xx, p.d) / (xx);
}

static double F(double gp, double ge, double e)
{
    double factor = FINE_STRUCTURE_CONSTANT * ELECTRON_RADIUS * ELECTRON_RADIUS / 2;
    double ret;
    double error;
    double d = ge / gp;

    double lo = 0.5 * (d + 1/d) + 1;
    double hi = 2 * gp * e;

    if(lo >= hi) return 0;

    struct G_params p = (struct G_params) {.d = d };

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(256);
    gsl_function G = (gsl_function) { .function = &_G, .params = &p };

    gsl_integration_qags(&G, log(lo), log(hi), 0, 1e-6, 256, w, &ret, &error);

    gsl_integration_workspace_free(w);

    return factor / d * ret;
}

void bethe_heitler_ka_process_lepton_gains(state_t *st)
{
    unsigned int i, j, k;

    double factor = LIGHT_SPEED;
    double dlngp = st->protons.log_energy[1] - st->protons.log_energy[0];
    double dlnx  = st->photons.log_energy[1] - st->photons.log_energy[0];

    double gains_inner[st->protons.size];

    for(i = 0; i < st->electrons.size; i++)
    {
        unsigned int index_base1 = i * st->protons.size;

        for(j = 0; j < st->protons.size; j++)
        {
            unsigned index_e_min = st->bethe_heitler_ka_LUT_e_min_index[index_base1 + j];
            unsigned index_e_max = st->photons.size - 1;

            gains_inner[j] = 0;

            if(index_e_min >= index_e_max) continue;

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            gains_inner[j] += st->photons.population[index_e_min] *
                              st->bethe_heitler_ka_LUT_reaction_rate[index_base2 + index_e_min] /
                              st->photons.energy[index_e_min];
            gains_inner[j] += st->photons.population[index_e_max] *
                              st->bethe_heitler_ka_LUT_reaction_rate[index_base2 + index_e_max] /
                              st->photons.energy[index_e_max];

            for(k = index_e_min + 1; k < index_e_max; k++)
            {
                double e = st->photons.energy[k];
                double n = st->photons.population[k];
                double r = st->bethe_heitler_ka_LUT_reaction_rate[index_base2 + k];

                gains_inner[j] += 2 * n * r / e;
            }

            gains_inner[j] *= dlnx / 2;
        }

        double gains = 0;

        for(j = 0; j < st->protons.size; j++)
        {
            double np = st->protons.population[j];
            double gp = st->protons.energy[j];

            gains += np / (gp*gp) * gains_inner[j];
        }

        gains *= dlngp / 2;

        st->bethe_heitler_lepton_gains[i] = factor * gains;
    }
}

void init_bethe_heitler_ka_LUT_lepton_gains(state_t *st)
{
    size_t size1 = st->electrons.size * st->protons.size;
    size_t size2 = st->electrons.size * st->protons.size * st->photons.size;

    MEM_PREPARE(st->bethe_heitler_ka_LUT_e_min_index,   size1, unsigned int);
    MEM_PREPARE(st->bethe_heitler_ka_LUT_reaction_rate, size2, double);
}

void calculate_bethe_heitler_ka_LUT_lepton_gains(state_t *st)
{
    unsigned int i, j, k;

    /* First, calculate the index of e_min */
    for(i = 0; i < st->electrons.size; i++)
    {
        double ge = st->electrons.energy[i];

        unsigned int index_base1 = i * st->protons.size;

        for(j = 0; j < st->protons.size; j++)
        {
            double gp = st->protons.energy[j];

            double delta = ge / gp;
            double e_min = (0.5 * (delta + 1 / delta) + 1) / (2 * gp);

            for(k = 0; k < st->photons.size && st->photons.energy[k] < e_min; k++)
            { }

            st->bethe_heitler_ka_LUT_e_min_index[index_base1 + j] = k;
        }
    }

    FILE *table;

    table = fopen("./tables/BH_KA_table.tab", "rb");

    if(table != NULL)
    {
        unsigned int e_size;
        unsigned int p_size;
        unsigned int g_size;

        double *e_energy;
        double *p_energy;
        double *g_energy;

        double *reaction_rate;

        fread(&e_size,  sizeof(unsigned int), 1,      table);
        fread(&p_size,  sizeof(unsigned int), 1,      table);
        fread(&g_size,  sizeof(unsigned int), 1,      table);

        e_energy = malloc(e_size * sizeof(double));
        p_energy = malloc(p_size * sizeof(double));
        g_energy = malloc(g_size * sizeof(double));
        reaction_rate = calloc(e_size * p_size * g_size, sizeof(double));

        fread(e_energy, sizeof(double),       e_size, table);
        fread(p_energy, sizeof(double),       p_size, table);
        fread(g_energy, sizeof(double),       g_size, table);

        fread(reaction_rate, sizeof(double), e_size * p_size * g_size, table);

        fclose(table);

        /* Let us check the original table, just in case */
        if(0)
        {
            for(unsigned int t = 0; t < 5; t++)
            {
                i = rand() / (double) RAND_MAX * e_size;
                j = rand() / (double) RAND_MAX * p_size;
                k = rand() / (double) RAND_MAX * g_size;

                double ge = exp(e_energy[i]);
                double gp = exp(p_energy[j]);
                double e  = exp(g_energy[k]);
                double r  = exp(reaction_rate[(i * p_size + j) * g_size + k]);

                double delta = ge / gp;
                double e_min = (0.5 * (delta + 1/delta) + 1) / (2 * gp);
                fprintf(stderr,"%03u %03u %03u (%lg, %lg, %lg-%lg) %lg <-> %lg\n",
                        i, j, k, ge, gp, e_min, e, r, F(gp, ge, e));
            }
        }

        for(i = 0; i < e_size * p_size * g_size; i++)
        {
            if(!isnormal(reaction_rate[i]))
                reaction_rate[i] = log(DBL_MIN);
        }

        /* NOTE!!!!
         *
         * We will do two interpolations using gsl.
         * First, a bicubical interpolation in a series of 'layers' until
         * we get a 'line' and then along this 'line' to obtain the values
         * we need.
         *
         * XXX: Note also that we store our values in a rather different
         * way as how GSL expects everything to be ordered so we need to be a
         * bit careful as how we init all interpolation functions for
         * everything to work.
         *
         * we: i*p_size*g_size + j*g_size + k;
         *
         * so i is the index for \gamma_e
         *    j is the index for \gamma_p
         *    k is the index for \varepsilon
         *
         * In GSL, things are sorted as
         *    j * xsize + i
         *
         * so, j=j and i=k
         */

        gsl_spline2d *spline[e_size];
        gsl_interp_accel *xacc[e_size];
        gsl_interp_accel *yacc[e_size];

        gsl_spline       *z_spline;
        gsl_interp_accel *zacc;


        for(i = 0; i < e_size; i++)
        {
            /*spline[i] = gsl_spline2d_alloc(gsl_interp2d_bicubic, g_size, p_size);*/
            spline[i] = gsl_spline2d_alloc(gsl_interp2d_bilinear, g_size, p_size);
            xacc[i]   = gsl_interp_accel_alloc();
            yacc[i]   = gsl_interp_accel_alloc();
        }

        double z_interpolated[e_size];
        z_spline = gsl_spline_alloc(gsl_interp_steffen, e_size);
        zacc = gsl_interp_accel_alloc();

        for(i = 0; i < e_size; i++)
        {
            gsl_spline2d_init(spline[i], g_energy, p_energy,
                    reaction_rate + (i * g_size * p_size),
                    g_size, p_size);
        }

        for(j = 0; j < st->protons.size; j++)
        {
            double lgp = st->protons.log_energy[j];

            for(k = 0; k < st->photons.size; k++)
            {
                double le = st->photons.log_energy[k];

                for(i = 0; i < e_size; i++)
                {
                    if(le  < g_energy[0] || le  > g_energy[g_size - 1] ||
                       lgp < p_energy[0] || lgp > p_energy[p_size - 1])
                        z_interpolated[i] = log(DBL_MIN);
                    else
                        z_interpolated[i] = gsl_spline2d_eval(spline[i], le, lgp, xacc[i], yacc[i]);
                }

                gsl_spline_init(z_spline, e_energy, z_interpolated, e_size);

                for(i = 0; i < st->electrons.size; i++)
                {
                    double lge = st->electrons.log_energy[i];
                    double r;

                    if(lge < e_energy[0] || lge > e_energy[e_size - 1])
                        r = log(DBL_MIN);
                    else
                        r = gsl_spline_eval(z_spline, lge, zacc);

                    // TODO: CHECK THIS WITH A BETTER TABLE
                    /*if(!isfinite(r))*/
                        /*r = log(DBL_MIN);*/

                    if(!isfinite(r))
                    {
                        fprintf(stderr,"Error interpolating\n");
                        fprintf(stderr,"%u %u %u %lg\n", i, j, k, r);
                        break;
                    }

                    st->bethe_heitler_ka_LUT_reaction_rate[
                        (i * st->protons.size + j) * st->photons.size + k] = exp(r);
                }
            }
        }

        gsl_spline_free(z_spline);
        gsl_interp_accel_free(zacc);

        for(i = 0; i < e_size; i++)
        {
            gsl_spline2d_free(spline[i]);
            gsl_interp_accel_free(xacc[i]);
            gsl_interp_accel_free(yacc[i]);
        }

        free(e_energy);
        free(p_energy);
        free(g_energy);
        free(reaction_rate);

        return;
    }

    for(i = 0; i < st->electrons.size; i++)
    {
        double ge = st->electrons.energy[i];

        unsigned int index_base1 = i * st->protons.size;

        for(j = 0; j < st->protons.size; j++)
        {
            double gp = st->protons.energy[j];

            double delta = ge / gp;
            double e_min = (0.5 * (delta + 1 / delta) + 1) / (2 * gp);

            unsigned int index_base2 = (index_base1 + j) * st->photons.size;

            unsigned int index_e_min = 0;

            for(k = 0; k < st->photons.size && st->photons.energy[k] < e_min; k++)
            { }

            index_e_min = k;

            for(/* empty */ ; k < st->photons.size; k++)
            {
                double e = st->photons.energy[k];

                double r = F(gp, ge, e);

                st->bethe_heitler_ka_LUT_reaction_rate[index_base2 + k] = r;

                if(!isnormal(r) || r < 0)
                {
                    fprintf(stderr,"%u %u %u %lg\n", i, j, k, r);
                    break;
                }

                /*fprintf(stderr,"Done %u %u %u %lg %lg %lg\n", i, j, k, ge, gp, e);*/
            }
            fprintf(stderr,"Done %u %u %lg %lg\n", i, j, ge, gp);
        }

        fprintf(stderr,"Done %u %lg\n", i, ge);
    }
}
