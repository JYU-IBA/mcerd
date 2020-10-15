#include "general.h"
#include "prototypes.h"
#include <jibal_gsto.h>
#include <jibal_stop.h>

#define NSTO 500 /* Keep this below or at most equal to MAXSTO */

void calc_stopping_and_straggling(Global *global, Ion *ion, Target *target, int nlayer) {
    Target_layer *layer;
    double minv, maxv, vstep, v;
    int i, j, p, z1, z2, nion;

    nion = ion->scatindex;

    layer = &(target->layer[nlayer]);


    minv = 0.0; /* Don't change */
    maxv = 1.1 * sqrt(2.0 * global->ionemax / (ion->A));
    vstep = (maxv - minv) /
            NSTO; /* I would use (NSTO-1) in the divider, but since this value is used in "stodiv", which is used to compute the index in an evenly spaced stopping table this convention will have to be respected for now. */
    z1 = (int) (ion->Z + 0.5);

    for (i = 0; i < MAXSTO; i++)
        layer->sto[nion].sto[i] = 0.0;

    layer->sto[nion].stodiv = 1.0 / vstep;
    layer->sto[nion].n_sto = NSTO;

    for (i = 0; i < layer->natoms; i++) {
        p = layer->atom[i];
        z2 = target->ele[p].Z;
        if (!jibal_gsto_auto_assign(global->jibal->gsto, z1, z2)) {
            fprintf(stderr, "Error in assigning stopping Z1=%i, Z2=%i.\n", z1, z2);
            exit(14);
        }
    }
    if (!jibal_gsto_load_all(global->jibal->gsto)) {
        fprintf(stderr, "Error in loading stopping.\n");
        exit(14);
    }

    for (i = 0; i < layer->natoms; i++) {
        p = layer->atom[i];
        z2 = target->ele[p].Z;
        for (j = 0; j < layer->sto[nion].n_sto; j++) {
            v = j * vstep;
            layer->sto[nion].vel[j] = v;
            double em = energy_per_mass(v);
            double stop = jibal_gsto_get_em(global->jibal->gsto, GSTO_STO_ELE, z1, z2, em) * layer->N[i];
            double stragg = jibal_gsto_get_em(global->jibal->gsto, GSTO_STO_STRAGG, z1, z2, em) * layer->N[i];
            if ((j && stop == 0.0) || !isfinite(stop)) {
                fprintf(stderr, "Warning: stopping at v=%e m/s (j=%i) for atom %i (i=%i), Z1=%i, Z2=%i is: %g\n", v,
                        j, p, i, z1, z2, stop);
            }
            if ((j && stragg == 0.0) || !isfinite(stragg)) {
                fprintf(stderr, "Warning: straggling at v=%e m/s (j=%i) for atom %i (i=%i), Z1=%i, Z2=%i is: %g\n", v,
                        j, p, i, z1, z2, stragg);
            }
            layer->sto[nion].sto[j] += stop;
            layer->sto[nion].stragg[j] += stragg;
        }
    }
}
