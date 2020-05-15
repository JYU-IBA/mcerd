/*
 *  get_espe
 *
 *   Copyright (C) 1996-2020  Kai Arstila, 2015-2020 Jaakko Julin
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <jibal.h>
#include "gauss.h"

#define NDATA   100000 /* Initial number of events for malloc. More space for events will be realloced as necessary. */
#define NSPE      20000
#define NDDIST      500
#ifdef TRACKID
#define NLINE 256
#else
#define NLINE       128
#endif

#define MCERD_SCALING_DEPTH 10*C_NM

#define BEAM_DOSE   1.0*C_UC
#define BEAM_ANGLE  (41.12*C_DEG)
#define BEAM_Z      17
#define BEAM_A        35
#define BEAM_M      (34.969*C_U)
#define BEAM_ENERGY (10.0*C_MEV)

#define TARGET_ANGLE   (20.0*C_DEG)
#define TARGET_DENSITY (4.98e22/C_CM3)

#define DET_SOLID       (0.6*0.226e-3)
#define DET_TIMERES     (0.170*C_NS/C_FWHM)
#define DET_ERES        (15.0*C_KEV)
#define DET_TOFLEN       0.623
#define DET_CIRCULAR 0
#define DET_RECT 1
#ifdef TRACKID
#define NEVENT_DATA  15
#else
#define NEVENT_DATA  11
#endif

#define TRUE  1
#define FALSE 0

typedef struct {
    char sca;
    char det;
    char sct;
    int scale;
    int real;
    int recoil;
    double energy;
    int z2;
    int a2;
    double m2;
    double depth;
    double weight;
    double time;
    double dx;
    double dy;
#ifdef TRACKID
    int64_t trackid;
    int status;
    int layer;
    double t1;
    double t2;
#endif
} Event;

typedef struct {
    int real;
    double ch;
    int ddist;
    int dwindow;
    int detsize;
    int avemass;
    int efromtime;
    int z1;
    int a1;
    int err;
    double m1;
    double timeres;
    double eres;
    double toflen;
    double m2;
    double dose;
    double theta;
    double tangle;
    double energy;
    double solid;
    double density;
    double depmin;
    double depmax;
    double scale;
    double dx;
    double dy;
    double emin;
    double emax;
    int output_events;
    int detgeo;
    int nbins;
} Params;

typedef struct {
    int n;
    double x[NDDIST];
    double y[NDDIST];
} Ddist;

typedef struct {
    Ddist ddist;
    double avemass;
    double masssum;
    double scaleweightsum; /* Sum of weights of all scaling ions */
    double scaleweight; /* Above, divided by nscale */
    double realweightsum; /* Sum of weights of all real ions before depth distribution weights */
    double realweight; /* Above, divided by nreal */
    double drealweightsum; /* Sub of weights of all read ions after depth distribution weights */
    double drealweight;
    int nevents;
    int nerrors;
    int nscale;
    int nreal;
    int nreal_participating;
    int ndout;
} Data;

typedef struct {
    double cs;
    double total;
    double final;
    double m2; /* m2 used for cross section calculation */
    int z2; /* z2 used for cross sectoin calculation */
} Scale;

void readparams(int, char **, jibal *, Params *, Data *);
void read_ddist(char *, Ddist *);
double get_ddist_value(double, Ddist *);
void usage(void);
int convert_and_check_validity(Event *, Data *, Params *);

double ipow2(double);
double ipow(double, int);

double Srbs(double, double, double, double, double, double);
double Serd(double, double, double, double, double, double);
double Srbs_mc(double, double, double, double);
int cm_angle(double, double, double *);
double mc2lab_scatc(double, double, double);
int print_events(Event *, int, double);
void read_element(jibal_isotope *isotopes, const char *sym, int *z, double *m, int *a);


void read_element(jibal_isotope *isotopes, const char *sym, int *z, double *m, int *a) {
    const jibal_isotope *i = jibal_isotope_find(isotopes, sym, 0, 0);
    if (!i) {
        fprintf(stderr, "ERROR: can not find isotope %s\n", sym);
        return;
    }
    *z = i->Z;
    *m = i->mass;
    *a = i->A;
}

int print_events(Event *event, int nevents, double scale) {
    int i;
    for (i = 0; i < nevents; i++) {
        if (!(event[i].scale)) {
#ifdef TRACKID
            fprintf(stdout, "%c %c %c " /* 1, 2, 3 */
                            "%12"PRIu64" %i %i " /* 4, 5, 6 */
                            "%8.4f %3i %6.2f " /* 7, 8, 9 */
                            "%10.4f %14.7e " /* 10, 11 */
                            "%7.3f %7.3f " /* 12, 13 */
                            "%6.2f %6.2f\n" /* 14, 15 */
                event[i].sca, event[i].det, event[i].sct, /* 1, 2, 3 */
                event[i].trackid, event[i].status, event[i].layer, /* 4, 5, 6 */
                event[i].energy/C_MEV, event[i].z2, event[i].m2/C_U, /* 7, 8, 9 */
                event[i].depth/C_NM, event[i].weight*scale,  /* 10, 11 */
                event[i].t1, event[i].t2, /* 12, 13 */
                event[i].dx/C_MM, event[i].dy/C_MM /* 14, 15 */
                );
#else
            fprintf(stdout, "%c %c %c " /* 1, 2, 3 */
                            "%8.4f %3i %6.2f " /* 4, 5, 6 */
                            "%10.4f %14.7e %7.3f " /* 7, 8, 9 */
                            "%6.2f %6.2f\n", /* 10, 11 */
                    event[i].sca, event[i].det, event[i].sct, /* 1, 2, 3 */
                    event[i].energy / C_MEV, event[i].z2, event[i].m2 / C_U, /* 4, 5, 6 */
                    event[i].depth / C_NM, event[i].weight * scale, event[i].time,  /* 7, 8, 9 */
                    event[i].dx / C_MM, event[i].dy / C_MM /* 10, 11 */
            );
#endif
        }
    }
    return 1;
}

void output_spectrum(Params *params, Event *event, int recoil, const Data *data, double scale) {
    int i, j, maxj = 0;
    double emin, emax, time, ene, intensity = 0.0, m2;
    double *hist_counts = calloc((size_t) params->nbins, sizeof(double));
    double *hist_weight = calloc((size_t) params->nbins, sizeof(double));
    emin = 1000 * C_MEV;
    emax = -1.0 * C_MEV;

    for (i = 0; i < data->nevents; i++) {
        if (params->avemass) {
            if (params->m2 > 0.0)
                m2 = params->m2;
            else
                m2 = data->avemass;
        } else {
            m2 = event[i].m2;
        }
        if (!recoil)
            m2 = params->m1;
        if (params->efromtime) {
            time = event[i].time + params->timeres * gaussian();
            event[i].energy = 0.5 * m2 * ipow2(params->toflen / time);
        } else {
            event[i].energy += params->eres * gaussian();
        }

        if (event[i].energy < emin)
            emin = event[i].energy;
        if (event[i].energy > emax)
            emax = event[i].energy;
    }

    if (params->emax > 0.0) {
        emin = params->emin;
        emax = params->emax;
    }

    emin = (floor(emin / params->ch) - 1) * params->ch;
    emax = (floor(emax / params->ch) + 1) * params->ch;

    fprintf(stderr, "Minimum energy:                      %10.3f MeV\n", emin / C_MEV);
    fprintf(stderr, "Maximum energy:                      %10.3f MeV\n", emax / C_MEV);

    for (i = 0; i < data->nevents; i++) {
        if (!(event[i].scale)) {
            j = (int) (0.5 + (event[i].energy - emin) / params->ch);
            if (j > 0 && j < params->nbins - 1) {
                hist_weight[j] += event[i].weight;
                (hist_counts[j])++;
                if (j > maxj)
                    maxj = j;
            }
        }
    }
    maxj++; /* Empty last bin in output */
    if (params->emax > 0.0) {
        maxj = (int) (0.5 + (emax - emin) / params->ch) + 1;
    }

    for (i = 0; i <= maxj; i++) {
        ene = (params->ch * i + emin);
        intensity += scale * hist_weight[i];
        if (params->err) {
            printf("%10.3f %10.3f %10.3f %6.0f\n", ene / C_MEV, scale * hist_weight[i],
                   scale * hist_weight[i] / sqrt(hist_counts[i] + 1), hist_counts[i]);
        } else {
            printf("%10.3f %10.3f\n", ene / C_MEV, scale * hist_weight[i]);
        }
    }
    fprintf(stderr, "Total spectrum intensity:        %12.1f\n", intensity);
    free(hist_counts);
    free(hist_weight);
}

Scale calculate_scale(const Params *params, const Data *data, int z2, int recoil) {
    Scale s;
    s.z2 = z2;
    if (params->m2 > 0.0) {
        s.m2 = params->m2;
    } else {
        s.m2 = data->avemass;
    }
    if (recoil) {
        s.cs = Serd(params->z1, params->m1, s.z2, s.m2, params->theta, params->energy);
    } else {
        s.cs = Srbs(params->z1, params->m1, s.z2, s.m2, params->theta, params->energy);
    }
    s.total = s.cs * params->solid * (params->density/MCERD_SCALING_DEPTH) * params->dose * (1.0 / sin(params->tangle));
    if (params->scale > 0.0) {
        s.final = params->scale;
    } else {
        s.final = s.total / data->scaleweightsum;
    }
    return s;
}

void calculate_data_averages(Data *data) {
    data->avemass = data->masssum / data->nevents;
    data->scaleweight = data->scaleweightsum / data->nscale;
    data->realweight = data->realweightsum / data->nreal;
    data->drealweight = data->drealweightsum / data->nreal;
}

int main(int argc, char *argv[]) {
    Event *event;
    Params params = {0};
    Data data = {0};
    Scale scale;
    jibal jibal;
    char line[NLINE];
    int i;
    int valid;
    int z2, recoil, nerr;
    int nevents_allocated = NDATA;
    event = (Event *) malloc(sizeof(Event) * nevents_allocated);
    for (i = 0; i < argc; i++) {
        fprintf(stderr, "%s%s", argv[i], i < argc - 1 ? " " : "\n");
    }
    z2 = 0;
    recoil = TRUE;

    jibal = jibal_init(NULL);
    if (jibal.error) {
        fprintf(stderr, "Initializing JIBAL failed with error code: %i (%s)\n", jibal.error,
                jibal_error_string(jibal.error));
        return -1;
    }

    readparams(argc, argv, &jibal, &params, &data);
    fprintf(stderr, "Beam energy:                           %7.2f MeV \n", params.energy / C_MEV);
    fprintf(stderr, "Detector angle:                        %7.2f deg \n", params.theta / C_DEG);
    fprintf(stderr, "Target angle:                          %7.2f deg \n", params.tangle / C_DEG);
    fprintf(stderr, "Detector solid angle:                  %7.2f msr\n", params.solid / C_MSR);
    fprintf(stderr, "Surface density for scaling:              %6.3e 1/cm3\n", params.density * C_CM3);
    fprintf(stderr, "MCERD scaling ion max depth:              %5.1e nm\n", MCERD_SCALING_DEPTH/C_NM);
    fprintf(stderr, "Beam dose:                              %10.2e (%.3f particle-uC)\n", params.dose,
            params.dose / C_UC);
    if (params.scale > 0.0) {
        fprintf(stderr, "Preset scaling factor:                %11.8f\n", params.scale);
    }

    nerr = 0;

    int lineno = 0;
    while (fgets(line, NLINE, stdin) != NULL) {
        lineno++;
        if (*line == '#')
            continue;
        Event *e = &event[data.nevents];
#ifdef TRACKID
#define
        int ncols = sscanf(line,"%c %c %c %"PRIi64" %i %i %lf %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &(e->sca),&(e->det),&(e->sct), /* 1, 2, 3 */
                   &(e->trackid), &(e->status), &(e->layer), /* 4, 5, 6 */
                   &(e->energy),&(e->z2),&(e->m2), /* 7, 8, 9 */
                   &(e->depth),&(e->weight), /* 10, 11 */
                   &(e->t1), &(e->t2),  /* 12, 13 */
                   &(e->dx),&(e->dy), /* 14, 15 */
        );
        event->time=event->t2-event->t1;
#else
        int ncols = sscanf(line, "%c %c %c %lf %i %lf %lf %lf %lf %lf %lf",
                           &(e->sca), &(e->det), &(e->sct),
                           &(e->energy), &(e->z2), &(e->m2),
                           &(e->depth), &(e->weight), &(e->time),
                           &(e->dx), &(e->dy));
#endif
        if (ncols != NEVENT_DATA) {
            fprintf(stderr, "Wrong number of columns (%i != %i) on line %i: %s", ncols, NEVENT_DATA, lineno, line);
            return -1;
        }

        valid = convert_and_check_validity(e, &data, &params);
        if (valid && data.nevents == 1) {
            recoil = e->recoil;
            z2 = event->z2;
        }

        if (data.nevents >= nevents_allocated) {
            nevents_allocated *= 2;
            event = (Event *) realloc(event, sizeof(Event) * nevents_allocated);
            if (!event) {
                fprintf(stderr, "Allocation failure nevents_allocated=%i", nevents_allocated);
            }
        }
#if 0
        if(data.nevents && (data.nevents % 10000 == 0)) {
            calculate_data_averages(&data);
            scale = calculate_scale(&params, &data, z2, recoil);
        }
#endif
    }

    if (data.nevents == 0) {
        fprintf(stderr, "Could not read any valid events.\n");
        exit(1);
    }

    calculate_data_averages(&data);

    fprintf(stderr, "Total events:   %7i\n", data.nevents);
    fprintf(stderr, "Real events:    %7i\n", data.nreal);
    fprintf(stderr, "Scaling events: %7i\n", data.nscale);
#ifdef DEBUG
    fprintf(stderr,"Zero:           %7i\n", data.nevents-data.nreal-data.nscale);
#endif
    fprintf(stderr, "Errors:         %7i\n", nerr);
    fprintf(stderr, "Discarded:      %7i\n\n", data.ndout);


    scale = calculate_scale(&params, &data, z2, recoil);

    fprintf(stderr, "Z1, M1, A1:   %3i %6.2f u %3i\n", params.z1, params.m1 / C_U, params.a1);
    fprintf(stderr, "Z2, M2:       %3i %6.2f u\n\n", scale.z2, scale.m2 / C_U);

    fprintf(stderr, "%s cross section for scaling ions:    %7.2f b/sr\n", recoil ? "ERD" : "RBS", scale.cs / C_BARN);
    fprintf(stderr, "Average mass of recoils:                 %5.2f u\n", data.avemass / C_U);
    fprintf(stderr, "Expected number of scaling ions: %12.1f\n", scale.total);
    fprintf(stderr, "Weight of scaling ions:          %12.1f\n", data.scaleweightsum);
    fprintf(stderr, "Average weight of scaling ions:     %10.2f\n", data.scaleweight);
    fprintf(stderr, "Average weight of real ions:        %10.2f\n", data.realweight);
    fprintf(stderr, "Avg weight with depth distribution: %10.2f\n", data.drealweight);
    fprintf(stderr, "Final scaling factor:                    %11.8f\n", scale.final);


    if (params.output_events) {
        print_events(event, data.nevents, scale.final);
    } else {
        output_spectrum(&params, event, recoil, &data, scale.final);
    }

    jibal_free(&jibal);
    free(event);
    exit(0);
}

double ipow2(double x) {
    return (x * x);
}

double ipow(double x, int n) {
    double value = 1.0;
    int i;

    for (i = 0; i < n; i++) {
        value *= x;
    }

    return (value);
}

void readparams(int argc, char *argv[], jibal *jibal, Params *params, Data *data) {
    int i = 1, step, err = FALSE;

    params->real = FALSE;
    params->ch = (0.1 * C_MEV);
    params->dwindow = FALSE;
    params->ddist = FALSE;

    params->err = FALSE;

    params->efromtime = TRUE;
    params->avemass = FALSE;
    params->m2 = 0.0;

    params->solid = DET_SOLID;
    params->eres = DET_ERES;
    params->timeres = DET_TIMERES;
    params->toflen = DET_TOFLEN;

    params->z1 = BEAM_Z;
    params->a1 = BEAM_A;
    params->m1 = BEAM_M;
    params->dose = BEAM_DOSE;
    params->theta = BEAM_ANGLE;
    params->energy = BEAM_ENERGY;
    params->tangle = TARGET_ANGLE;
    params->density = TARGET_DENSITY;

    params->emin = 1000 * C_MEV;
    params->emax = -1.0 * C_MEV;

    params->scale = -1.0;
    params->dx = -1.0;
    params->dy = -1.0;

    params->output_events = FALSE;
    params->detgeo = DET_CIRCULAR;

    params->nbins = NSPE; /* TODO: make configurable */

#ifdef DEBUG
    for(i=0; i < argc; i++) {
        fprintf(stderr, "argv[%i]=\"%s\"%s", i, argv[i], i!=argc-1?", ":"\n");
    }
#endif


    for (i = 1; i < argc && !err; i += step) {
        step = 0;
        if (!strcmp(argv[i], "-real")) {
            params->real = TRUE;
            step = 1;
        }
        if (!strcmp(argv[i], "-ch")) {
            if ((i + 1) < argc) {
                params->ch = atof(argv[i + 1]);
                if (!(params->ch > 0.0 && params->ch <= 10.0))
                    err = TRUE;
                if (!err)
                    params->ch *= C_MEV;
                step = 2;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-depth")) {
            if ((i + 2) < argc) {
                params->dwindow = TRUE;
                params->depmin = atof(argv[i + 1]);
                params->depmax = atof(argv[i + 2]);
                if (!(params->depmin >= 0.0 && params->depmax >= 0.0))
                    err = TRUE;
                params->depmin *= C_NM;
                params->depmax *= C_NM;
                step = 3;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-dist")) {
            if ((i + 1) < argc) {
                read_ddist(argv[i + 1], &(data->ddist));
                params->ddist = TRUE;
                step = 2;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-m2")) {
            params->avemass = TRUE;
            if ((i + 1) < argc) {
                params->m2 = atof(argv[i + 1]);
                if (!(params->m2 > 0.0 && params->m2 < 300.0))
                    err = TRUE;
                params->m2 *= C_U;
                step = 2;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-avemass")) {
            params->avemass = TRUE;
            step = 1;
        }
        if (!strcmp(argv[i], "-scale")) {
            if ((i + 1) < argc) {
                params->scale = atof(argv[i + 1]);
                if (params->scale <= 0.0)
                    err = TRUE;
                step = 2;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-err")) {
            params->err = TRUE;
            step = 1;
        }
        if (!strcmp(argv[i], "-detsize")) {
            if ((i + 2) < argc) {
                params->dx = atof(argv[i + 1]);
                params->dy = atof(argv[i + 2]);
                if (!(params->dx > 0.0 && params->dx <= 10000.0))
                    err = TRUE;
                if (!(params->dy > 0.0 && params->dy <= 10000.0))
                    err = TRUE;
                params->dx *= C_MM;
                params->dy *= C_MM;
                step = 3;
                if (!err)
                    params->detsize = TRUE;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-erange")) {
            if ((i + 2) < argc) {
                params->emin = atof(argv[i + 1]);
                params->emax = atof(argv[i + 2]);
                if (!(params->emin >= 0.0 && params->emin <= 1000.0))
                    err = TRUE;
                if (!(params->emax > 0.0 && params->emax <= 1000.0))
                    err = TRUE;
                if (params->emin >= params->emax)
                    err = TRUE;
                params->emin *= C_MEV;
                params->emax *= C_MEV;
                step = 3;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-timeres")) {
            if ((i + 1) < argc) {
                params->efromtime = TRUE;
                params->timeres = atof(argv[i + 1]);
                if (!(params->timeres >= 0.0 && params->timeres <= 10000.0))
                    err = TRUE;
                params->timeres *= (C_PS / C_FWHM);
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-eres")) {
            if ((i + 1) < argc) {
                params->efromtime = FALSE;
                params->eres = atof(argv[i + 1]);
                if (!(params->eres > 0.0 && params->eres <= 100000.0))
                    err = TRUE;
                params->eres *= (C_KEV / C_FWHM);
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-toflen")) {
            if ((i + 1) < argc) {
                params->toflen = atof(argv[i + 1]);
                if (!(params->toflen > 0.0 && params->toflen <= 50.0))
                    err = TRUE;
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-beam")) {
            if ((i + 1) < argc) {
                read_element(jibal->isotopes, argv[i + 1], &(params->z1), &(params->m1), &(params->a1));
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-energy")) {
            if ((i + 1) < argc) {
                params->energy = atof(argv[i + 1]);
                if (!(params->energy > 0.0 && params->energy <= 1000.0))
                    err = TRUE;
                params->energy *= C_MEV;
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-dose")) {
            if ((i + 1) < argc) {
                params->dose = atof(argv[i + 1]);
                if (!(params->dose > 0.0 && params->dose <= 1.0e25))
                    err = TRUE;
                if (params->dose < 1000.0)
                    params->dose *= (1.0e-6 / C_EV);
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-theta")) {
            if ((i + 1) < argc) {
                params->theta = atof(argv[i + 1]);
                if (!(params->theta > 0.0 && params->theta <= 180.0))
                    err = TRUE;
                params->theta *= C_DEG;
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-tangle")) {
            if ((i + 1) < argc) {
                params->tangle = atof(argv[i + 1]);
                if (!(params->tangle > 0.0 && params->tangle <= 90.0)) {
                    fprintf(stderr, "Target angle (tangle) must be > 0.0 and <= 90.0 degrees! You gave: %lf\n",
                            params->tangle);
                    err = TRUE;
                }
                params->tangle *= C_DEG;
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-solid")) {
            if ((i + 1) < argc) {
                params->solid = atof(argv[i + 1]);
                params->solid *= C_MSR;
                if (!(params->solid > 0.0 && params->solid <= 4.0 * C_PI))
                    err = TRUE;
                step = 2;
            } else
                err = TRUE;
        }

        if (!strcmp(argv[i], "-density")) {
            if ((i + 1) < argc) {
                params->density = atof(argv[i + 1]);
                if (params->density <= 0.0)
                    err = TRUE;
                params->density /= C_CM3;
                step = 2;
            } else
                err = TRUE;
        }
        if (!strcmp(argv[i], "-convert")) {
            params->output_events = TRUE;
            step = 1;
        }

        if (!strcmp(argv[i], "-detgeo")) {
            if ((i + 1) < argc) {
                if (!strcmp(argv[i + 1], "rect")) {
                    params->detgeo = DET_RECT;
                }
                step = 2;
            } else
                err = TRUE;
        }

        if (err || step == 0) {
            fprintf(stderr, "Error with argument %i: %s\n", i, i < argc ? argv[i] : "i out of bounds");
            usage();
        }
    }
}

void read_ddist(char *fname, Ddist *ddist) {
    FILE *fp;
    int i = 0;

    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", fname);
        exit(1);
    }

    while (i < NDDIST && fscanf(fp, "%lf %lf", &(ddist->x[i]), &(ddist->y[i])) == 2) {
        ddist->x[i] *= C_NM;
        i++;
    }

    if (i == 0) {
        fprintf(stderr, "Could not read any depth data from file %s\n", fname);
        exit(2);
    }
    if (i >= NDDIST) {
        fprintf(stderr, "Too many data points (> %i) in file %s\n", NDDIST, fname);
        exit(3);
    }

    ddist->n = i;

    fclose(fp);
}

double get_ddist_value(double depth, Ddist *d) {
    double value;
    int n, i = 0;

    n = d->n;

    if (depth <= d->x[0] || depth >= d->x[n - 1])
        return (0.0);

    while (depth > d->x[i])
        i++;

    value = d->y[i - 1] + (depth - d->x[i - 1]) * (d->y[i] - d->y[i - 1]) /
                          (d->x[i] - d->x[i - 1]);

    return (value);

}

void usage(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "get_espe - Calculate an energy spectrum from simulated ERD data\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-real    only real events are handled\n");
    fprintf(stderr, "\t-ch      channel width in the output (MeV)\n");
    fprintf(stderr, "\t-depth   depth range of the events (nm, two values)\n");
    fprintf(stderr, "\t-dist    file name for depth distribution\n");
    fprintf(stderr, "\t-m2      average mass of the secondary particles (u)\n");
    fprintf(stderr, "\t-avemass use average mass for calculating energy from TOF\n");
    fprintf(stderr, "\t-scale   scale the total intensity to value\n");
    fprintf(stderr, "\t-err     give statistics in the third column\n");
    fprintf(stderr, "\t-detsize limit in the size of the detector foil (mm)\n");
    fprintf(stderr, "\t-erange  energy range in the output spectrum (MeV)\n");

    fprintf(stderr, "\t-timeres time resolution of the TOF-detector (ps, FWHM)\n");
    fprintf(stderr, "\t-eres    energy resolution (keV, FWHM) of the SSD, (energy signal used!)\n");
    fprintf(stderr, "\t-toflen  time-of-flight length (m)\n");

    fprintf(stderr, "\t-beam    mass number and the chemical symbol of the primary ion\n");
    fprintf(stderr, "\t-dose    dose of the beam (particle-uC)\n");
    fprintf(stderr, "\t-energy  beam energy (MeV)\n");
    fprintf(stderr, "\t-theta   scattering angle (deg)\n");
    fprintf(stderr, "\t-tangle  angle between target surface and beam (deg)\n");
    fprintf(stderr, "\t-solid   solid angle of the detector (msr)\n");
    fprintf(stderr, "\t-density surface atomic density of the first 10 nm layer (at./cm3), e.g. 4.98e22 for Si.\n");
    exit(10);
}

int convert_and_check_validity(Event *event, Data *data, Params *params) {
    double dweight;
    int nerr = 0, valid = TRUE;

    if (event->sca == 'S') {
        event->scale = TRUE;
    } else if (event->sca == 'R') {
        event->scale = FALSE;
    } else {
        fprintf(stderr, "Not a scaling ion (S) nor real one (R)\n");
        valid = FALSE;
        nerr++;
    }

    if (event->det == 'R') {
        event->real = TRUE;
    } else if (event->det == 'V') {
        event->real = FALSE;
    } else {
        fprintf(stderr, "Not a real ion (R) nor virtual one (V)\n");
        valid = FALSE;
        nerr++;
    }

    if (event->sct == 'R') {
        event->recoil = TRUE;
    } else if (event->sct == 'S') {
        event->recoil = FALSE;
    } else {
        fprintf(stderr, "Not a recoil atom (R) nor a scattered one (S)\n");
        valid = FALSE;
        nerr++;
    }

    if (event->energy < 0.0 || event->energy > 1000.0) {
        fprintf(stderr, "Particle energy must be 0 - 1000 MeV\n");
        valid = FALSE;
        nerr++;
    } else {
        event->energy *= C_MEV;
    }

    if (event->m2 < 1.0 || event->m2 > 300.0) {
        fprintf(stderr, "M2 must be 1 - 300\n");
        valid = FALSE;
        nerr++;
    } else {
        event->a2 = (int) (event->m2 + 0.5);
        event->m2 *= C_U;
    }

    if (event->z2 < 1 || event->z2 > 100) {
        fprintf(stderr, "Z2 must be 1 - 100\n");
        valid = FALSE;
        nerr++;
    }

    if (event->depth < 0.0) {
        fprintf(stderr, "Depth must be positive\n");
        valid = FALSE;
        nerr++;
    } else {
        event->depth *= C_NM;
    }

    if (event->weight < 0.0) {
        fprintf(stderr, "Weight must be positive and I have some issues with this: %g\n", event->weight);
        valid = FALSE;
        nerr++;
    }
    if (event->weight == 0.0) {
        fprintf(stderr, "WARNING: event with EXACTLY (%g) weight!\n", event->weight);
    }

    if (event->time <= 0.0 && event->scale) {
        fprintf(stderr, "TOF time must be positive (for scaling ions). It is %g.\n", event->time);
        valid = FALSE;
        nerr++;
    } else {
        event->time *= C_NS;
    }

    event->dx *= C_MM;
    event->dy *= C_MM;

    if (params->real) {
        if (!(event->real)) {
            valid = FALSE;
            data->ndout++;
        }
    }

    if (params->detsize) {
        switch (params->detgeo) {
            case DET_CIRCULAR:
                if (ipow2(event->dx / params->dx) + ipow2(event->dy / params->dy) > 1.0) {
                    valid = FALSE;
                    data->ndout++;
                }
            case DET_RECT:
                if (fabs(event->dx / params->dx) > 1.0 || fabs(event->dy / params->dy)) {
                    valid = FALSE;
                    data->ndout++;
                }
                break;
            default:
                break;

        }
    }

    if (params->dwindow && !event->scale) {
        if (event->depth < params->depmin || event->depth > params->depmax) {
            fprintf(stderr, "Depth outside range from %g to %g.\n", params->depmin, params->depmax);
            valid = FALSE;
            data->ndout++;
        }
    }

    if (params->ddist && !event->scale) {
        dweight = get_ddist_value(event->depth, &(data->ddist));
        if (dweight > 0) {
            data->realweightsum += event->weight;
            data->nreal++;
            event->weight *= dweight;
            data->drealweightsum += event->weight;
        } else {
            fprintf(stderr, "Depth distribution thinks here (depth=%g) should be nothing.\n", event->depth);
            event->weight = 0.0;
            valid = FALSE;
            data->ndout++;
        }
    }

    if (valid && event->scale) {
        data->scaleweightsum += event->weight;
        data->nscale++;
    }

    if (valid) {
        data->masssum += event->m2;
        data->nevents++;
    }

    data->nerrors += nerr;

    return (valid);

}

double Srbs(double z1, double m1, double z2, double m2, double t, double E) {
    double value, Ecm, tcm[2];
    int n;

    Ecm = m2 * E / (m1 + m2);
    n = cm_angle(t, m1 / m2, tcm);

    value = mc2lab_scatc(Srbs_mc(z1, z2, tcm[0], Ecm), tcm[0], t);

    return (value);
}

double Srbs_mc(double z1, double z2, double t, double E) {
    double value;

    value = ipow2((z1 * z2 * C_E * C_E) / (4.0 * C_PI * C_EPSILON0)) * ipow2(1.0 / (4.0 * E)) *
            ipow(1.0 / sin(t / 2.0), 4);

    return (value);
}

double Serd(double z1, double m1, double z2, double m2, double t, double E) {
    double value;

    value = ipow2(z1 * z2 * C_E * C_E / (8 * C_PI * C_EPSILON0 * E)) * ipow2(1.0 + m1 / m2) /
            ipow(cos(t), 3);

    return (value);
}

int cm_angle(double lab_angle, double r, double a[]) {
    double stmp;

    if (r > 1.0 && sin(lab_angle) > 1.0 / r) {
        return (0);
    }

    stmp = asin(r * sin(lab_angle));

    if (r > 1.0) {
        a[0] = lab_angle + stmp;
        a[1] = C_PI + lab_angle - stmp;
        return (2);
    } else {
        a[0] = lab_angle + stmp;
        return (1);
    }

}

double mc2lab_scatc(double mcs, double tcm, double t) {
    double value;

    value = (mcs * ipow2(sin(tcm))) / (ipow2(sin(t)) * cos(tcm - t));

    return (value);
}
