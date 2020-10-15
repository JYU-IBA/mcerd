#include "general.h"
#include "prototypes.h"

void init_params(Global *global, Target *target, int argc, char *argv[]) {
    int i, j;

    global->master.argv = argv;
    global->master.argc = argc;

    global->ionemax = -1.0;
    global->emin = 10 * C_KEV;
    global->E0 = -1.0;
    global->nsimu = 10;
    global->predata = FALSE;

    global->beamprof = BEAM_NONE;
    global->beamdiv = 0.0;

    global->nions = 2; /* Will be adjusted later */
    global->ncascades = 10; /* Default, can be overridden */

    global->nmclarge = 0;
    global->nmc = 0;

    global->costhetamin = cos(1.0 * C_DEG);

    target->natoms = 0;
    target->nlayers = 0;
    target->natoms = 0;

    for (i = 0; i <= SECONDARY; i++)
        for (j = 0; j < NIONSTATUS; j++)
            global->finstat[i][j] = 0;

    global->master.argv = argv;

    global->rough = FALSE;
    global->output_misses = FALSE;
    global->output_trackpoints = FALSE;
    global->cascades = FALSE;
    global->advanced_output = FALSE;
    global->nomc = FALSE; /* REMINDER: default should be FALSE */

}

void init_recoiling_angle(Global *global, Ion *ion, Target *target, Detector *detector) {
    double angle[50][2], depth[50][2], low, high, rlow, rhigh, tlow, thigh;
    double angave, lensum, a, b, x0, x1;
    int i = 1, j = 0, k = 0, finished = FALSE, recd_nonzero, layer[50], delay = 0, n;

    rlow = target->recdist[i - 1].x;
    rhigh = target->recdist[i].x;
    tlow = target->layer[j].dlow;
    thigh = target->layer[j].dhigh;

    while (!(tlow <= rlow && thigh >= rlow)) {
        j++;
        delay = 1;
        tlow = target->layer[j].dlow;
        thigh = target->layer[j].dhigh;
    }

    high = tlow;

    do {
        low = high;
        if (target->recdist[i - 1].y > 0 || target->recdist[i].y > 0)
            recd_nonzero = TRUE;
        else
            recd_nonzero = FALSE;
        if (rhigh < thigh) {
            high = rhigh;
            rlow = rhigh;
            i++;
            if (i < target->nrecdist)
                rhigh = target->recdist[i].x;
            else
                finished = TRUE;
        } else {
            high = thigh;
            tlow = thigh;
            j++;
            delay = 1;
            if (j < target->ntarget)
                thigh = target->layer[j].dhigh;
            else
                finished = TRUE;
        }
        if (!finished && recd_nonzero) {
            n = j - delay;
            depth[k][0] = low;
            depth[k][1] = high;
            layer[k] = n;
            angle[k][0] = low * target->recpar[n].x + target->recpar[n].y;
            angle[k][1] = high * target->recpar[n].x + target->recpar[n].y;
            k++;
        }
        delay = 0;
    } while (!finished);

    if (j >= target->ntarget && recd_nonzero) {
        n = j - delay;
        depth[k][0] = low;
        depth[k][1] = high;
        layer[k] = n;
        angle[k][0] = low * target->recpar[n].x + target->recpar[n].y;
        angle[k][1] = high * target->recpar[n].x + target->recpar[n].y;
        k++;
        delay = 0;
    }

    angave = 0.0;
    lensum = 0.0;
    for (i = 0; i < k; i++) {
        for (j = 0; j < 2; j++)
            printf("%3i %10.3f %10.3f\n", layer[i], depth[i][j] / C_NM, angle[i][j] / C_DEG);
        printf("\n");

        n = layer[i];
        a = target->recpar[n].x;
        b = target->recpar[n].y;
        x0 = depth[i][0];
        x1 = depth[i][1];
        angave += ((1. / 3.) * a * a * (x0 * x0 + x0 * x1 + x1 * x1) + a * b * (x0 + x1) + b * b) *
                  (x1 - x0);
        lensum += x1 - x0;
    }

    target->angave = angave / lensum;

    printf("angave %10.7f\n", target->angave / C_DEG);

}

void init_io(Global *global, Ion *ion, Target *target) {
    Target_layer *layer;
    char fname[NFILE], label[NFILE];
    double Emin, vel, x;
    int i, j;

    strcpy(global->basename, global->master.argv[1]);
    sprintf(label, ".%u", global->seed);
    strcat(global->basename, label);

    sprintf(fname, "%s.%s", global->basename, "out");

    global->master.fpout = fopen(fname, "w");

    if (global->master.fpout == NULL)
        fatal_error("Could not open a file for general output\n");

    sprintf(fname, "%s.%s", global->basename, "dat");

    global->master.fpdat = fopen(fname, "w");

    if (global->master.fpdat == NULL)
        fatal_error("Could not open file for formatted output\n");

    sprintf(fname, "%s.%s", global->basename, "erd");

    global->master.fperd = fopen(fname, "w");
    if (global->master.fperd == NULL)
        fatal_error("Could not open ERD output file\n");

    sprintf(fname, "%s.%s", global->basename, "range");

    global->master.fprange = fopen(fname, "w");

    if (global->master.fprange == NULL)
        fatal_error("Could not open range output file\n");

    if (global->output_trackpoints) {
        sprintf(fname, "%s.%s", global->basename, "track");
        global->master.fptrack = fopen(fname, "w");

        if (global->master.fptrack == NULL)
            fatal_error("Could not open track output file\n");
    }

#if 0

    Emin = global->emin;

    j = 0;
    for(i=0;i<target->ntarget;i++){
       layer = &(target->layer[i]);
       for(x=layer->dlow;x < layer->dhigh && j<NRECDIST;x+=C_NM){
          vel = sqrt(2.0*Emin/ion[SECONDARY].A);
          Emin += inter_sto(&(layer->sto[SECONDARY]),vel,STOPPING)*C_NM;
          target->efin[j] = Emin;
 /*         printf("%8.3f %8.3f\n",x/C_NM,Emin/C_MEV); */
          j++;
       }
    }

#endif

#ifdef DEBUG
    sprintf(fname,"%s.%s",global->basename,"debug");

    global->master.fpdebug = fopen(fname,"w");

    if(global->master.fperd == NULL)
       fatal_error("Could not open file for debugging output\n");
#endif
}
