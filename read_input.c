#include "general.h"
#include "prototypes.h"
#include "read_input.h"
#include "rough_surface.h"

void read_input(Global *global, Ion **ion, Target *target, Detector *detector) {
    FILE *fp, *fout;
    Point p1, p2, p3;
    Fvalue *fval;
    char buf[MAXLEN], *c, *word;
    double number, unit;
    double rec_dist_unit, M;
    int i, j, rinput[NINPUT], n;

    for (i = 0; i < NINPUT; i++)
        rinput[i] = FALSE;

    fout = stderr;

    if (global->master.argc > 1) {
        fprintf(stderr, "Reading configuration from %s\n", global->master.argv[1]);
        fp = fopen(global->master.argv[1], "r");
        if (fp == NULL)
            fatal_error("Could not open command file\n");
    } else {
        fatal_error("No command file given");
        exit(11);
    }

    int lineno = 0;
    while (fgets(buf, MAXLEN, fp) != NULL) {
        lineno++;
        buf[strcspn(buf, "\r\n")] = 0; /* Strips newlines */
        i = 0;
        while (i < NINPUT && (c = strstr(buf, inlines[i])) == NULL)
            i++;
        fprintf(stderr, "%3i %s\n", lineno, buf);
        if (i < NINPUT) {
            c += strlen(inlines[i]);
            c = trim_space(c);
            fprintf(stderr, "Line %i, got inline #%i: %s %s\n", lineno, i, inlines[i], c);
        } else {
            fprintf(stderr, "Unknown configuration line %i, ignoring it: %s", lineno, buf);
        }

        switch (i) {
            case I_TYPE:
                word = get_word(c, &n);  /* we should free word-memory also! */
                if (strcmp(word, "ERD") == 0) {
                    global->simtype = SIM_ERD;
                    global->nions = 2; /* Incident and recoil */
                } else if (strcmp(word, "RBS") == 0) {
                    global->simtype = SIM_RBS;
                    global->nions = 3;  /* Incident, target and scattered incident */
                } else {
                    fatal_error("no such type for simulation");
                }
                *ion = malloc(sizeof(Ion) * global->nions);
                if (!*ion) {
                    fatal_error("Can not allocate ions.");
                }
                memset(*ion, 0, sizeof(Ion) * global->nions);
                break;
            case I_ION:
                fprintf(stderr, "Trying to find ion %s\n", c);
                get_ion(global->jibal->elements, c, &((*ion)[PRIMARY].Z), &((*ion)[PRIMARY].A), &((*ion)[PRIMARY].I));
                (*ion)[PRIMARY].type = PRIMARY;
                fprintf(fout, "Beam ion: Z=%.0f, M=%.3f\n", (*ion)[PRIMARY].Z, (*ion)[PRIMARY].A / C_U);
                break;
            case I_ENERGY:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_MEV);
                global->E0 = number * unit;
                global->ionemax = 1.01 * global->E0;
                break;
            case I_TARGET:
                fprintf(stderr, "Reading target file.\n");
                read_target_file(c, global, target);
                break;
            case I_DETECTOR:
                fprintf(stderr, "Reading detector file.\n");
                read_detector_file(c, global, detector, target);
                break;
            case I_RECOIL:
                get_ion(global->jibal->elements, c, &((*ion)[SECONDARY].Z), &((*ion)[SECONDARY].A),
                        &((*ion)[SECONDARY].I));
                (*ion)[SECONDARY].type = SECONDARY;
                fprintf(fout, "Recoil atom: Z=%.0f, M=%.3f\n", (*ion)[SECONDARY].Z, (*ion)[SECONDARY].A / C_U);
                fprintf(fout, "Recoil atom isotopes:\n");
                for (i = 0; i < (*ion)[SECONDARY].I.n; i++)
                    fprintf(fout, "%10.3f %10.3f%%\n", (*ion)[SECONDARY].I.A[i] / C_U,
                            (*ion)[SECONDARY].I.c[i] * 100.0);
                if ((*ion)[SECONDARY].I.n > 0)
                    fprintf(fout, "Most abundant isotope: %10.3f\n", (*ion)[SECONDARY].I.Am / C_U);
                break;
            case I_RECDIST:
                fval = read_file(c, 2, &n);
                if (n < 2)
                    fatal_error("Too few points in the recoiling material distribution\n");
                if (fval->a[0] < 0.0)
                    fatal_error("Recoil distribution can not start from a negative depth\n");
                rec_dist_unit = C_NM;
                target->recmaxd = fval->a[n - 1] * rec_dist_unit;
                target->effrecd = 0.0;
                for (i = 0; i < n; i++) {
                    target->recdist[i].x = fval->a[i] * rec_dist_unit;
                    target->recdist[i].y = fval->b[i] * rec_dist_unit;
                    if (i > 0 && (target->recdist[i - 1].y > 0.0 ||
                                  target->recdist[i].y > 0.0)) {
                        target->effrecd += (target->recdist[i].x - target->recdist[i - 1].x);
                    }
                }
                target->nrecdist = n;
                free(fval);
                break;
            case I_TANGLE:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_DEG);
                global->beamangle = PI / 2.0 - number * unit;
                break;
            case I_SPOTSIZE:
                c = get_number(c, &number);
                global->bspot.x = 0.5 * number;
                c = get_number(c, &number);
                global->bspot.y = 0.5 * number;
                get_unit_value(c, &unit, C_KEV);
                global->bspot.x *= unit;
                global->bspot.y *= unit;
                break;
            case I_MINANGLE:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_DEG);
                global->minangle = number * unit;
                break;
            case I_MINENERGY:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_KEV);
                global->emin = number * unit;
                break;
            case I_NIONS:
                c = get_number(c, &number);
                global->nsimu = number;
                break;
            case I_NPREIONS:
                c = get_number(c, &number);
                global->npresimu = number;
                break;
            case I_RECAVE:
                c = get_number(c, &number);
                global->nrecave = number;
                break;
            case I_SEED:
                c = get_number(c, &number);
                global->seed = number;
                rnd(0.0, 1.0, RND_SEED, -1 * (long) global->seed);
                break;
            case I_RECWIDTH:
                word = get_word(c, &n);
                if (strcmp(word, "wide") == 0 || strcmp(word, "WIDE") == 0) {
                    fprintf(stderr, "Got wide.\n");
                    global->recwidth = REC_WIDE;
                } else if (strcmp(word, "narrow") == 0 || strcmp(word, "NARROW") == 0) {
                    fprintf(stderr, "Got narrow.\n");
                    global->recwidth = REC_NARROW;
                } else {
                    fprintf(stderr, "Recoil angle width should be either wide or narrow.");
                    exit(666);
                }
                break;
            case I_PREDATA:
                fval = read_file(c, 2, &n);
                for (i = 0; i < n; i++) {
                    target->recpar[i].x = fval->a[i] * C_DEG / C_NM;
                    target->recpar[i].y = fval->b[i] * C_DEG;
                }
                free(fval);
                global->predata = TRUE;
                break;
            case I_MINSCAT:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_DEG);
                global->costhetamin = cos(number * unit);
                break;
            case I_BDIV:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_DEG);
                global->beamdiv = cos(number * unit);
                break;
            case I_BPROF:
                word = get_word(c, &n);
                if (strcmp(word, "flat") == 0 || strcmp(word, "FLAT") == 0)
                    global->beamprof = BEAM_FLAT;
                else
                    global->beamprof = BEAM_GAUSS;
                break;
            case I_SURFACE:
                word = get_word(c, &n);
                read_afm(word, &(target->surface));
                global->rough = TRUE;
                break;
            case I_SURFSIZE:
                c = get_number(c, &number);
                get_unit_value(c, &unit, C_DEG);
                target->surface.size = number * unit;
                break;
            case I_NSCALE:
                c = get_number(c, &number);
                global->nscale = number;
                break;
            case I_TRACKP:
                word = get_word(c, &n);
                if (strcmp(word, "true") == 0 || strcmp(word, "TRUE") == 0) {
                    global->output_trackpoints = TRUE;
                }
                break;
            case I_MISSES:
                word = get_word(c, &n);
                if (strcmp(word, "true") == 0 || strcmp(word, "TRUE") == 0) {
                    global->output_misses = TRUE;
                }
                break;
            case I_CASCADES:
                word = get_word(c, &n);
                if (strcmp(word, "true") == 0 || strcmp(word, "TRUE") == 0) {
                    global->cascades = TRUE;
                }
                break;
            case I_NCASCADES:
                c = get_number(c, &number);
                global->ncascades = min(global->nions, max(number, 1000));
                break;
            case I_ADVANCED:
                word = get_word(c, &n);
                if (strcmp(word, "true") == 0 || strcmp(word, "TRUE") == 0) {
                    global->advanced_output = TRUE;
                }
                break;
            case I_NOMC:
                word = get_word(c, &n);
                if (strcmp(word, "true") == 0 || strcmp(word, "TRUE") == 0) {
                    global->nomc = TRUE;
                }
                break;
            default:
                break;
        }
    }


    p1.x = p1.y = p1.z = 0.0;
    p2 = p1;
    p3 = p1;

    p2.y = 1.0;
    p3.x = sin(PI / 2.0 - global->beamangle);
    p3.z = cos(PI / 2.0 - global->beamangle);

    target->plane = get_plane_params(&p1, &p2, &p3);

    if (global->recwidth == REC_WIDE && global->npresimu > 0) {
        global->npresimu = 0;
        fprintf(stderr, "Presimulation not needed with wide recoil angle scheme\n");
    }

    if (global->recwidth == REC_WIDE)
        global->predata = FALSE;

    if (global->predata) {
        global->npresimu = 0;
    } else if (global->npresimu == 0 && global->recwidth == REC_NARROW) {
        fatal_error("No precalculated recoiling data in the narrow recoiling scheme\n");
    }

    if (global->npresimu > 0) {
        global->cpresimu = 0;
        global->simstage = PRESIMULATION;

        global->presimu = (Presimu *)
                malloc(global->npresimu * global->nrecave * 2 * sizeof(Presimu));

        if (global->presimu == NULL)
            fatal_error("Could not allocate enough memory for the presimulation\n");
    } else {
        global->simstage = REALSIMULATION;
    }

    target->ntarget = target->nlayers - detector->nfoils;

    detector->tdet[0] += target->ntarget;
    detector->tdet[1] += target->ntarget;
    detector->edet[0] += target->ntarget;

    for (i = target->nlayers - 1; i > target->ntarget; i--) {
        target->layer[i].dlow = 0.0;
        target->layer[i].dhigh -= target->layer[i - 1].dhigh;
    }

    global->nsimu += global->npresimu;

    M = 4.0 * (*ion)[PRIMARY].A * (*ion)[SECONDARY].A /
        ipow2((*ion)[PRIMARY].A + (*ion)[SECONDARY].A);

    if (global->simtype == SIM_ERD) {
        global->costhetamax = sqrt(global->emin / (global->E0 * M));
        global->costhetamin = 1.0;
    } else if (global->simtype == SIM_RBS) {
        if ((*ion)[PRIMARY].A <= (*ion)[SECONDARY].A) {
            global->costhetamax = -1.0;
        } else {
            global->costhetamax = sqrt(1.0 - ipow2((*ion)[SECONDARY].A / (*ion)[PRIMARY].A));
        }
        /* costhetamax is calculated for average recoil mass! */
    }
    fclose(fp);

}

Fvalue *read_file(char *fname, int cols, int *n) {
    FILE *fp;
    Fvalue *f;
    char buf[LINE];
    int i = 0;

    f = (Fvalue *) malloc(sizeof(Fvalue));

    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", fname);
        exit(10);
    }

    while (fgets(buf, LINE, fp) != NULL) {
        switch (cols) {
            case 1:
                sscanf(buf, "%lf", &(f->a[i]));
                break;
            case 2:
                sscanf(buf, "%lf %lf", &(f->a[i]), &(f->b[i]));
                break;
            case 3:
                sscanf(buf, "%lf %lf %lf", &(f->a[i]), &(f->b[i]), &(f->c[i]));
                break;
        }
        i++;
    }
    fclose(fp);

    *n = i;

    return (f);

}

void get_atom(const jibal_element *elements, char *c, double *z) {
    char *symbol;
    int n;
    symbol = get_word(c, &n);
    const jibal_element *element = jibal_element_find(elements, symbol);
    if (!element) {
        fprintf(stderr, "Element %s not found from database\n", symbol);
        exit(12);
    }
    *z = (double) element->Z;
}

void get_ion(const jibal_element *elements, char *c, double *z, double *m, Isotopes *I) {
    char *symbol;
    int mass = -1, found = FALSE, n;
    double number, Amax = -1.0;

    if (isdigit(*c)) {
        c = get_number(c, &number);
        mass = (int) (number + 0.5);
    }
    symbol = get_word(c, &n);
    n = 0;
    jibal_element *element = jibal_element_copy(jibal_element_find(elements, symbol), JIBAL_NAT_ISOTOPES);
    if (!element) {
        fprintf(stderr, "Element %s not found from database\n", symbol);
        exit(12);
    }
    *z = element->Z;
    int i;
    I->c_sum = 0.0;
    for (i = 0; i < element->n_isotopes; i++) {
        const jibal_isotope *isotope = element->isotopes[i];
        double conc = element->concs[i];
        if (mass > 0.0) { /* Try to find a single isotope */
            if (mass == isotope->A) {
                *m = isotope->mass;
                I->c_sum = 1.0;
                found = TRUE;
                break;
            }
        } else {
            found = TRUE;
            I->A[n] = isotope->mass;
            I->c[n] = conc;
            I->c_sum += conc;
            if (isotope->abundance > Amax) {
                Amax = isotope->abundance;
                I->Am = isotope->mass;
                *m = I->Am;
            }
            n++;
        }

    }
    I->n = n;
    if (n && I->c_sum != 1.0) {
        fprintf(stderr, "WARNING: Concentrations of element %s don't sum exactly to 100%%, but instead to %lf*100\n",
                symbol, I->c_sum * 100.0);
    }
    free(element);
    if (!found) {
        fprintf(stderr, "Isotope(s) for %s not found from database\n", c);
        exit(12);
    }
}

char *get_string(char *c) {

/* 
   String can contain white spaces, thus there can be only one string per line
*/
    char *p;
    size_t i = 0;

    while (*(c + i) != '\0')
        i++;

    i--;

    while (i > 0 && isspace(*(c + i)))
        i--;

    i++;

    if (*(c + i) == '\0')
        return (NULL);

    p = (char *) malloc(sizeof(char) * (i + 1));

    p = strncpy(p, c, i);

    p[i] = '\0';

    return (p);

}

char *get_word(char *c, int *n) {
    char *p;
    size_t i = 0;

    while (*(c + i) != '\0' && !isspace(*(c + i)))
        i++;

    i--;

    while (i > 0 && isspace(*(c + i)))
        i--;

    i++;

    if (i == 0)
        return (NULL);

    p = (char *) malloc(sizeof(char) * (i + 1));

    p = strncpy(p, c, i);

    p[i] = '\0';

    *n = i;
    return (p);

}

char *get_number(char *c, double *number) {
    int t;
    char *s;
    size_t i = 0;

    while (mc_isnumber(*(c + i)))
        i++;
    if ((*(c + i) == 'e' || *(c + i) == 'E') && i != 0 && *(c + i) != '\0'
        && *(c + i + 1) != '\0' &&
        (mc_isnumber(*(c + i + 1)) || *(c + i + 1) == '-' || *(c + i + 1) == '+')) {
        t = i;
        i++;
        while (mc_isnumber(*(c + i)))
            i++;
        if (t == (i + 1))
            i = t;
    }

    if (i == 0)
        return (NULL);
    if (*(c + i) == '\0') {
        *number = atof(c);
    } else {
        s = (char *) malloc((i + 1) * sizeof(char));
        strncpy(s, c, i);
        s[i] = '\0';
        *number = atof(s);
/*
      free(s);
*/
    }

    c = trim_space(c + i);
    return (c);
}

char *get_unit_value(char *c, double *unit, double def) {
    int j = -1, v;
    char *u;
    size_t i = 0;
/*   
   while(*(c+i) != '\0' && !mc_isnumber(*(c+i)) && !isspace(*(c+i)))
      i++;
*/
    while (*(c + i) != '\0' && !isspace(*(c + i)))
        i++;

    if (i == 0) {
        *unit = def;
        return (c);
    }

    u = (char *) malloc(sizeof(char) * (i + 1));
    strncpy(u, c, i);
    u[i] = '\0';

    do {
        j++;
        v = strcmp(u, units[j].unit);
    } while (units[j].value > 0.0 && v != 0);

    if (units[j].value < 0.0)
        fatal_error("Could not find value for given unit");

    *unit = units[j].value;

    c = trim_space(c + i);

    return (c);
}

char *trim_space(char *c) {
    while (*c != '\0' && isspace(*c))
        c++;

    return (c);
}

int mc_isnumber(char c) {
    if (isdigit(c) || c == '-' || c == '.' || c == '+')
        return (TRUE);
    else
        return (FALSE);

}
