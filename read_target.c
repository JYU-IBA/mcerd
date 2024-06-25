#include "general.h"
#include "prototypes.h"

#include "read_input.h"

int line_not_empty(char *);

void read_target_file(char *fname, Global *global, Target *target) {
    FILE *fp;
    char buf[LINE], *c, *stofile, *ltype, *state;
    double number, unit, Z, M, con[MAXELEMENTS], atom[MAXELEMENTS], sumcon = 0.0;
    double d = 0.0, density, minN = 1.0e100, sumM, sumN;
    int i, j, n = 0, nlayer, natoms, norigatom, noriglayer;

    fp = fopen(fname, "r");

    fprintf(stderr, "Opening target file \"%s\"\n", fname);
    if (fp == NULL) {
        fatal_error("Could not open the target description file\n");
    }

    natoms = norigatom = target->natoms;
    nlayer = noriglayer = target->nlayers;

    while (fgets(buf, LINE, fp) != NULL && line_not_empty(buf)) {
        c = get_number(buf, &M);
        M *= C_U;
        get_atom(global->jibal->elements, c, &Z);
        target->ele[natoms].Z = Z;
        target->ele[natoms].A = M;
        printf("Atom: %3i %5.4f %f\n", natoms, Z, M / C_U);
        natoms++;
    }
    target->natoms = natoms;

    while (fgets(buf, LINE, fp) != NULL) {
        printf("\nlayer: %i\n", nlayer);
        c = buf;
        target->layer[nlayer].type = TARGET_FILM;
        target->layer[nlayer].gas = FALSE;
        if (isalpha(*c)) {
            ltype = get_word(c, &n);
            c += n;
            c = trim_space(c);
            printf("layer type: %s\n", ltype);
            if (strstr(ltype, "bulk") != NULL) {
                target->layer[nlayer].type = TARGET_BULK;
                target->layer[nlayer].dlow = d;
                target->layer[nlayer].dhigh = 1.0;  /* 1.0 meters */
            }
        }
        if (target->layer[nlayer].type == TARGET_FILM) {
            c = get_number(c, &number);
            get_unit_value(c, &unit, C_NM);
            printf("thickness %.3f nm\n", number * unit / C_NM);
            target->layer[nlayer].dlow = d;
            target->layer[nlayer].dhigh = d + number * unit;
            d = target->layer[nlayer].dhigh;
        }

        fgets(buf, LINE, fp);
        c = buf;
        stofile = get_word(c, &n);
        printf("stofile: %s\n", stofile);
        if (!strcmp(stofile, "ZBL")) {
            *(target->layer[nlayer].stofile_prefix) = '\0';
        } else {
            strncpy(target->layer[nlayer].stofile_prefix, stofile, MAXSTOFILEPREFIXLEN);
        }
        fgets(buf, LINE, fp); /* Second "ZBL" most likely, ignored. TODO: do something about this */
        fgets(buf, LINE, fp);
        c = get_number(c, &number);
        if (!c) {
            fprintf(stderr, "ERROR! density: %s is not an acceptable value\n", buf);
        }
        get_unit_value(c, &unit, C_G_CM3);
        density = number * unit;
        printf("density: %.3f %s\n", density / C_G_CM3, SYM_G_CM3);

#ifdef GAS_OR_SOLID
        fgets(buf, LINE, fp);
        c = buf;
        state=get_word(c, &n);
        if(!strcmp(state, "gas")) {
            target->layer[nlayer].gas = TRUE;
        } else if(!strcmp(state, "solid")) {
            target->layer[nlayer].gas = FALSE;
        } else {
            fprintf(stderr, "Error in target (or foils). Specify whether the layer is gaseous \"gas\" or solid \"solid\""
                   " on a line immediately after specifying the density. Thank you. You gave: %s\n", buf);
            exit(-1);
        }
#endif
        printf(target->layer[nlayer].gas ? "gaseous\n" : "solid\n");

        i = 0;
        sumM = 0.0;
        sumcon = 0.0;
        while (fgets(buf, LINE, fp) != NULL && line_not_empty(buf)) {
            c = buf;
            c = get_number(c, atom + i);
            c = get_number(c, con + i);
            j = (int) (atom[i] + norigatom + 0.5);
            target->layer[nlayer].atom[i] = j;
            sumM += target->ele[j].A * con[i];
            sumcon += con[i];
            printf("atom: %i %.1f, con: %.3f %.3f\n", j, atom[i], con[i],
                   target->ele[j].A / C_U);
            i++;
        }
        sumM /= sumcon;
        n = i;
        target->layer[nlayer].natoms = n;
        sumN = density / sumM;
        target->layer[nlayer].Ntot = sumN;
        printf("sumM: %.3f u\nsumN %.4e %s\n", sumM / C_U, sumN, SYM_1_M3);
        for (i = 0; i < n; i++) {
            natoms = target->layer[nlayer].atom[i];
            con[i] /= sumcon;
            target->layer[nlayer].N[i] = con[i] * sumN;
            printf("%.2f %.0f %.4e %s\n", target->ele[natoms].A / C_U,
                   target->ele[natoms].Z, target->layer[nlayer].N[i], SYM_1_M3);
        }
        nlayer++;
    }
    target->nlayers = nlayer;

    for (i = 0; i < nlayer; i++) {
        if (target->layer[i].Ntot < minN) {
            target->minN = target->layer[i].Ntot;
            minN = target->minN;
        }
    }

    printf("\nnumber of layers: %i\nminN: %.4e", target->nlayers, target->minN);
    printf("\nnumber of atoms: %i\n", target->natoms);
    fclose(fp);

}

int line_not_empty(char *buf) {
    int empty = TRUE, i;

    for (i = 0; (i < LINE) && (buf[i] != '\0') && (empty == TRUE); i++)
        if (!isspace(buf[i]))
            empty = FALSE;

    return (!empty);
}
