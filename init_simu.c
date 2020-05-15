#include "general.h"
#include "prototypes.h"

void scattering_table(Global *global, Ion *ion, Target *target, Scattering *scat, Potential *pot, int natom) {
    double e, y, emin, emax, ymin, ymax, estep, ystep;
    double targetZ, targetA;
    int i, j;

    double angle;

    targetZ = target->ele[natom].Z;
    targetA = target->ele[natom].A;

    scat->a = 0.8854 * C_BOHR_RADIUS / (pow(ion->Z, 0.23) + pow(targetZ, 0.23));
    scat->E2eps = 4.0 * C_PI * C_EPSILON0 * scat->a * targetA /
                  ((ion->A + targetA) * ion->Z * targetZ * ipow2(C_E));
/*
   fprintf(global->master.fpout,"Z, A: %.0f %.5g\n",targetZ,targetA/C_U);
*/
    emin = log(global->emin * scat->E2eps);
    emax = log(global->ionemax * scat->E2eps);

    ymin = log(1.0 / (2.0 * exp(emax) * tan(0.5 * C_PI * MAXANGLE / 180.0)));
    ymax = log(1.0 / (2.0 * scat->a * pow(target->minN, 1. / 3.))); /* Think this! */

    estep = (emax - emin) / (EPSNUM - 2);
    ystep = (ymax - ymin) / (YNUM - 2);

    scat->logemin = emin;
    scat->logymin = ymin;

    scat->logediv = 1.0 / estep;
    scat->logydiv = 1.0 / ystep;

    e = emin;

    fprintf(global->master.fpout, "a, E2eps %e %e\n", scat->a, scat->E2eps);
    fprintf(global->master.fpout, "emin, emax: %g %g\n", emin, emax);
    fprintf(global->master.fpout, "ymin,ymax: %g %g\n", ymin, ymax);
    fprintf(global->master.fpout, "estep,ystep: %g %g\n", estep, ystep);
    fflush(global->master.fpout);
/*
   printf("Z,A: %10.3f %10.3f %10.3f %10.3f\n",ion->Z,ion->A/C_U,targetZ,targetA/C_U);
*/
    for (i = 1; i < EPSNUM; i++) {
        scat->angle[i][0] = exp(e);
        ion->opt.e = exp(e);
        y = ymin;
        for (j = 1; j < YNUM; j++) {
            ion->opt.y = exp(y);
            scat->angle[i][j] = scattering_angle(pot, ion);
/*
         printf("%10.5f %10.5f %10.5f\n",e,y,scat->angle[i][j]);
*/
            y += ystep;
        }
        e += estep;
    }

    y = ymin;
    for (j = 1; j < YNUM; j++) {
        scat->angle[0][j] = exp(y);
        y += ystep;
    }
/*
   printf("Z,A: %10.3f %10.3f %10.3f %10.3f\n",ion->Z,ion->A/C_U,targetZ,targetA/C_U);

   for(i=1;i<EPSNUM;i++){
      for(j=1;j<YNUM;j++){
         y = ymin + ystep*(j-1);
         y = (exp(y)*scat->a)/C_ANGSTROM;
         e = estep*(i-1) + emin;
         e = (exp(e)/scat->E2eps)/C_MEV;
         angle = scat->angle[i][j]/C_DEG;
         printf("%14.5e %14.5e %10.5f\n",y,e,angle);
      }
   }
*/
}
