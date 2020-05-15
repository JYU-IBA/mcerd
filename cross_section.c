#include "general.h"
#include "prototypes.h"

void calc_cross_sections(Global *global, Scattering *scat, Potential *pot) {
    int i;
    double e;

    scat->cross.emin = log(0.99 * global->emin * scat->E2eps);
    scat->cross.emax = log(1.01 * global->ionemax * scat->E2eps);

    scat->cross.estep = (scat->cross.emax - scat->cross.emin) / (EPSIMP - 1);

    e = scat->cross.emin;

    scat->cross.b = (double *) malloc(EPSIMP * sizeof(double));

    for (i = 0; i < EPSIMP; i++) {
        scat->cross.b[i] = calc_cross(global->minangle, exp(e), scat, pot);
/*
      printf("%10.4f %12.7f\n",exp(e)/(scat->E2eps*C_MEV),
             scat->cross.b[i]*scat->a/C_ANGSTROM);
*/
        e += scat->cross.estep;
    }


}

#define DISTEPS 5.0e-3
#define MAXSTEPS 50

double calc_cross(double angle, double e, Scattering *scat, Potential *pot) {
    Ion ion;
    double diff, ylogold, ylognew, alognew, alogold, y;
    double anew, aold, ynew, yold;
/*
 *  double alow,ahigh,maxdist,ylow,yhigh,ang;
 */
    int step = 0;

/*      
 *  We should be very very careful with maximum impact parameters, because
 *  setting a density proportional limit causes problems when comparing
 *  impact parameters with different target atoms.
 */

/*
   maxdist = 1.0/(2.0*scat->a*pow(atom->N,1./3.));

   maxdist = 5.0*C_ANGSTROM/scat->a;  
   ion.opt.e = e;
   ion.opt.y = maxdist;

   printf("Entering calc_cross-routine\n");
   printf("ion.e: %.7g, ion.y: %.7g\n",ion.opt.e,ion.opt.y);

   if(scattering_angle(pot,&ion) > angle){
      fatal_error("calc_cross: Impact parameter exceeds 5.0 A\n");
   }

   ion.opt.e = e;
   ylow = 1.0;
   yhigh = 3.0*C_ANGSTROM/scat->a;


   printf("ylow, yhigh: %.5g %.5g\n",ylow,yhigh);

   ion.opt.y = ylow;
   alow = scattering_angle(pot,&ion);
   ion.opt.y = yhigh;
   ahigh = scattering_angle(pot,&ion);

      printf("%5i %12.7g %12.7g %12.7f %12.7f\n",step,ylow,yhigh,alow/C_DEG,
        ahigh/C_DEG);
*/
    yold = 10.0;
    ylogold = log(yold);
    ynew = 2.0;
    ylognew = log(ynew);
    ion.opt.e = e;
    ion.opt.y = yold;
    aold = scattering_angle(pot, &ion);
    ion.opt.y = ynew;
    anew = scattering_angle(pot, &ion);
    alogold = log(aold);
    alognew = log(anew);
/*
   for(ynew=0.5;ynew<15;ynew+=0.1){
      ion.opt.y = ynew;
      anew = scattering_angle(pot,&ion);
      printf("%14.5e%14.5e\n",ynew,anew/C_DEG);
   }
   exit(0);
*/
    do {
/*
      printf("%5i %12.7g %12.7g %12.7f %12.7f\n",step+2,ynew,yold,anew/C_DEG,
        aold/C_DEG);
*/
        step++;
        y = exp(ylognew + (ylogold - ylognew) * (log(angle) - alognew) /
                          (alogold - alognew));
        if (y < 0.0)
            y = ynew / 2.0;
        ion.opt.y = y;

        yold = ynew;
        ylogold = ylognew;
        ynew = y;
        ylognew = log(y);
        aold = anew;
        alogold = alognew;
        anew = scattering_angle(pot, &ion);
        alognew = log(anew);

        diff = fabs(anew - angle) / angle;
    } while (diff > DISTEPS && step < MAXSTEPS);
/*
      printf("%5i %12.7g %12.7g %12.7f %12.7f\n",step+2,ynew,yold,anew/C_DEG,
        aold/C_DEG);


   printf("%15.10f\n",0.5*(anew+aold)/C_DEG);

   ion.opt.y = 0.5*(ylow + yhigh);

   printf("\n");
*/

    if (diff > DISTEPS)
        fprintf(stderr, "get_impact: convergence stopped \n");

/*
   printf("E,b: %.2f MeV, %.4f A\n",e/(scat->E2eps*C_MEV),ynew*scat->a/C_ANGSTROM);
*/
    return (ynew);
}

double get_cross(Ion *ion, Scattering *scat) {
/*
 *     Interpolate the cross section ie. the maximum impact parameter
 *     for current ion energy.
 *
 */

    double b, e;
    int i;

    e = log(ion->E * scat->E2eps);

/* 
 *  Check that the maximum impact parameter does not exceed the half of the
 *  atom distance. (Actually we do not do it here).
 */

    i = (int) ((e - scat->cross.emin) / scat->cross.estep);

    b = scat->cross.b[i] + (scat->cross.b[i + 1] - scat->cross.b[i]) *
                           (e - (i * scat->cross.estep + scat->cross.emin)) / scat->cross.estep;

    b *= scat->a;
    b = C_PI * ipow2(b);


    if (!(b > 0 && b < 1e-15)) {
        fprintf(stderr, "Cross section seems awfully low: \n");
        fprintf(stderr, "%e\n", b);
        fprintf(stderr, "%i %g %g\n", i, ion->E, scat->E2eps);
        fprintf(stderr, "%g %g %g\n", e, scat->cross.emin, scat->cross.estep);
        fprintf(stderr, "%g %g %g\n", scat->cross.b[i], scat->cross.b[i + 1], scat->a);
    }

    return (b);

}
