#include "general.h"
#include "prototypes.h"

#define NPRESIMU 500

int presimu_comp_angle(Presimu *, Presimu *);
int presimu_comp_depth(Presimu *, Presimu *);
void fit_linear(double *, double *, int, double *, double *);

void finish_presimulation(Global *global, Detector *detector, Ion *recoil) {
/*
 *    In the presimulation stage we save the recoil angle relative to the
 *    detector when the recoil comes out of the target. Also the recoil
 *    depth and layer are saved for later use.
 */

    double theta_lab, fii_lab, theta, fii;

    /* Recoil direction in the laboratory coordinate system */

    rotate(recoil->lab.theta, recoil->lab.fii, recoil->theta, recoil->fii,
           &theta_lab, &fii_lab);

    /* Recoil direction in the detector coordinate system */

    rotate(detector->angle, C_PI, theta_lab, fii_lab, &theta, &fii);

    global->presimu[global->cpresimu].depth = recoil->hist.tar_recoil.p.z;
    global->presimu[global->cpresimu].angle = theta;
    global->presimu[global->cpresimu].layer = recoil->hist.layer;
    global->cpresimu++;

}

void analyze_presimulation(Global *global, Target *target, Detector *detector) {
/*
 *    We determine here the solid angle as function of recoiling depth
 *    and layer, which contains all but PRESIMU_LEVEL portion of the
 *    recoils (typically 99%). The result is given as a linear fit
 *    for half-angle vs. depth for each target layer.
 */

    FILE *fp;
    char fname[NFILE];
    double sumdepth, depth[MAXLAYERS][NPRESIMU], angle[MAXLAYERS][NPRESIMU];
    double a = 0.0, b = 0.0;
    int i, j, n, npre = 0, nlevel, layer, nlayer[MAXLAYERS];

    for (i = 0; i < MAXLAYERS; i++)
        nlayer[i] = 0;

    nlevel = 10.0 / (PRESIMU_LEVEL);
    nlevel = max(nlevel, global->cpresimu / (NPRESIMU + 1));

    qsort((void *) (global->presimu), (size_t) (global->cpresimu),
          sizeof(Presimu), (int (*)(const void *, const void *))
                  presimu_comp_depth);

    while (npre < (global->cpresimu - nlevel / 2.0)) {
        layer = global->presimu[npre].layer;
        n = 0;
        sumdepth = 0;
        while (n < nlevel && (n + npre) < global->cpresimu && global->presimu[npre + n].layer == layer) {
            sumdepth += global->presimu[npre + n].depth;
            n++;
        }
        if (n > nlevel / 2.0) {
            qsort((void *) (global->presimu + npre), (size_t) (n),
                  sizeof(Presimu), (int (*)(const void *, const void *))
                          presimu_comp_angle);
            i = n * PRESIMU_LEVEL;
            depth[layer][nlayer[layer]] = sumdepth / n;
            angle[layer][nlayer[layer]] = global->presimu[npre + i].angle;
            nlayer[layer]++;
        }
        npre += n;
    }

    for (i = 0; i < target->ntarget; i++) {
        for (j = 0; j < nlayer[i]; j++)
            fprintf(global->master.fpout, "%3i %10.4f %10.4f\n", i, depth[i][j] / C_NM,
                    angle[i][j] / C_DEG);

        fprintf(global->master.fpout, "\n");
        fprintf(global->master.fpout, "%3i %10.5f %10.5f + 1.1 * %6.3f * %6.2f\n",
                i, a * C_NM / C_DEG, b / C_DEG, detector->thetamax / C_DEG,
                max(1.0, max(detector->vsize[0], detector->vsize[1])));
        if (nlayer[i] > 2) {
            fit_linear(depth[i], angle[i], nlayer[i], &a, &b);
            target->recpar[i].x = a;
            target->recpar[i].y = b +
                                  1.1 * detector->thetamax * max(1.0, max(detector->vsize[0], detector->vsize[1]));
        } else {
            if (i > 0) {
                target->recpar[i].x = target->recpar[i - 1].x;
                target->recpar[i].y = target->recpar[i - 1].y;
            } else {
                target->recpar[i].x = 0.0;
                target->recpar[i].y = 5.0 * C_DEG +
                                      1.1 * detector->thetamax * max(1.0, max(detector->vsize[0], detector->vsize[1]));
            }
        }
    }

    sprintf(fname, "%s.%s", global->basename, "pre");
    fp = fopen(fname, "w");
    if (fp == NULL)
        fatal_error("Could not open file for presimulation output\n");

    fprintf(global->master.fpout, "\n");
    for (i = 0; i < target->ntarget; i++) {
        fprintf(global->master.fpout, "%3i %10.5f %10.5f\n", i,
                target->recpar[i].x * C_NM / C_DEG, target->recpar[i].y / C_DEG);
        fprintf(fp, "%14.5e %14.5e\n", target->recpar[i].x * C_NM / C_DEG,
                target->recpar[i].y / C_DEG);
    }
    fclose(fp);

    fprintf(stderr, "Presimulation finished\n");

    free(global->presimu);

    global->simstage = REALSIMULATION;
}

int presimu_comp_angle(Presimu *a, Presimu *b) {
/*
 *    We want to sort in reverse order here.
 */

    if (a->angle < b->angle)
        return (1);
    else if (a->angle > b->angle)
        return (-1);
    else
        return (0);
}

int presimu_comp_depth(Presimu *a, Presimu *b) {
    if (a->depth < b->depth)
        return (-1);
    else if (a->depth > b->depth)
        return (1);
    else
        return (0);
}

void fit_linear(double *x, double *y, int n, double *a, double *b) {
    /*
     *    This is a weighted least-squares fit of a straight line
     *    according to P.R. Bevington's Data Reduction and Error
     *    Analysis for the Physical Sciences.
     *
     *    Here we use constantly weight 1.0 since were have same
     *    statistics for every point.
     */

    double Delta, xsum, ysum, xysum, x2sum, wsum, w = 1.0;
    int i;

    if (n < 3)
        fatal_error("Too few points in the LS-fit in the pre-simulation analysis\n");

    xsum = ysum = xysum = x2sum = wsum = 0.0;

    for (i = 0; i < n; i++) {
        xsum += x[i] * w;
        ysum += y[i] * w;
        xysum += x[i] * y[i] * w;
        x2sum += x[i] * x[i] * w;
        wsum += w;
    }

    Delta = 1.0 / (wsum * x2sum - xsum * xsum);

    *a = Delta * (wsum * xysum - xsum * ysum);

    *b = Delta * (x2sum * ysum - xsum * xysum);

}
