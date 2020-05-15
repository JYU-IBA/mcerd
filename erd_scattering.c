#include <assert.h>
#include "general.h"
#include "prototypes.h"

double get_isotope(Isotopes *);
double cross_section(double, double, double, double, double, double,
                     double [][NSANGLE], int, int);
double Serd(double, double, double, double, double, double);
double Srbs_mc(double, double, double, double);
double mc2lab_scatc(double, double, double);
int cm_angle(double, double, double *a);
double inter_cross_table(double, double, double [][NSANGLE]);

int erd_scattering(Global *global, Ion *ion, Ion *recoil, Target *target,
                   Detector *detector) {
    Ion *tatom;
    Vector sc_lab, sc_target, sc_ion;
    double sq_tmp, K;
    int i, table;

    if (ion->scale && ion->p.z > SCALE_DEPTH)
        return (FALSE);

    if (global->simtype ==
        SIM_RBS) /* FIXME: no pointer arithmetic on "ion", it is not a table! Although it works, most of the time... */
        tatom = ion + TARGET_ATOM;

    recoil->w = ion->w;

    if (ion->wtmp > 0.0)
        recoil->w *= ion->wtmp;

    if (global->recwidth == REC_NARROW) {

        /* Recoiling direction in the laboratory coordinate system */

        sc_lab = get_recoiling_dir(global, ion, recoil, target, detector);

        /* Recoiling direction in the target coordinate system */

        rotate(global->beamangle, 0.0, sc_lab.theta, sc_lab.fii,
               &(sc_target.theta), &(sc_target.fii));

        /* Recoiling direction in the ion coordinate system */

        rotate(ion->theta, ion->fii - C_PI, sc_target.theta, sc_target.fii,
               &(sc_ion.theta), &(sc_ion.fii));

    } else {

        /* Recoiling direction in the ion coordinate system */

        sc_ion = get_recoiling_dir(global, ion, recoil, target, detector);

        rotate(ion->theta, ion->fii, sc_ion.theta, sc_ion.fii,
               &(sc_target.theta), &(sc_target.fii));

        rotate(global->beamangle, C_PI, sc_target.theta, sc_target.fii,
               &(sc_lab.theta), &(sc_lab.fii));
    }

    if (recoil->w > 1e7)
        fprintf(stderr, "%.3e %.3e %.3e\n", recoil->w, ion->w, ion->wtmp);

    recoil->scale = ion->scale;

    if (global->simtype == SIM_ERD) {
        if (sc_ion.theta < C_PI / 2.0) {
            if (recoil->I.n > 0)
                recoil->A = get_isotope(&(recoil->I));
            recoil->E = ion->E * ipow2(cos(sc_ion.theta)) * 4.0 * (ion->A * recoil->A) /
                        ipow2(ion->A + recoil->A);
            recoil->theta = sc_target.theta;
            recoil->fii = sc_target.fii;
            recoil->opt.cos_theta = cos(recoil->theta);
            recoil->opt.sin_theta = sin(recoil->theta);
            recoil->p = ion->p;
            recoil->tlayer = ion->tlayer;
            recoil->nsct = 0;
            if (ion->scale) {
                table = FALSE; /* Rutherford cross sections for scaling ions! */
            } else {
                table = target->table;
            }
            recoil->w *= cross_section(ion->Z, ion->A, recoil->Z, recoil->A, ion->E,
                                       sc_ion.theta, target->cross, table, SIM_ERD);
            recoil->time = 0.0;
            recoil->type = SECONDARY;
            recoil->status = NOT_FINISHED;
            recoil->lab.p = ion->lab.p;
            recoil->lab.theta = ion->lab.theta;
            recoil->lab.fii = ion->lab.fii;
            recoil->virtual = FALSE;
/*         for(i=0;i<(target->nlayers - target->ntarget);i++) {
 	        recoil->Ed[i] = 0.0;
            recoil->hit[i].x = 0.0;
            recoil->hit[i].y = 0.0;
         }*/
#ifdef REC_ANG
            recoil->rec_ion_theta = sc_ion.theta/C_DEG;
            recoil->rec_lab_theta = sc_lab.theta/C_DEG;
            if(sc_ion.fii > PI)
               sc_ion.fii -= (2.0*PI);
            if(sc_lab.fii > PI)
               sc_lab.fii -= (2.0*PI);
            recoil->rec_ion_fii = sc_ion.fii/C_DEG;
            recoil->rec_lab_fii = sc_lab.fii/C_DEG;
#endif
            /*
             *   Save the physical state of the primary and the secondary ion
             *   at the scattering point.
             */

            save_ion_history(global, ion, recoil, detector, sc_lab, sc_ion,
                             recoil->Z, recoil->A);

#ifdef DEBUG
            print_ion_position(global,ion,"R",ANYSIMULATION);
#endif
            return (TRUE);
        } else {
            recoil->status = ion->status = FIN_NO_RECOIL;
            return (FALSE);
        }
    } else if (global->simtype == SIM_RBS) {
        if (tatom->I.n > 0)
            tatom->A = get_isotope(&(tatom->I));
        sq_tmp = ipow2(tatom->A) - ipow2(ion->A * sin(sc_ion.theta));

        if (sq_tmp > 0.0) {
            sq_tmp = sqrt(sq_tmp);
            K = ipow2((sq_tmp + ion->A * cos(sc_ion.theta)) / (ion->A + tatom->A));
            recoil->E = ion->E * K;
            if (ion->scale) {
                table = FALSE; /* Rutherford cross sections for scaling ions! */
            } else {
                table = target->table;
            }
            recoil->w *= cross_section(ion->Z, ion->A, tatom->Z, tatom->A, ion->E,
                                       sc_ion.theta, target->cross, table, SIM_RBS);

        } else { /* rare case because of the average masses */
            sq_tmp = 0.0;
            K = ipow2((sq_tmp + ion->A * cos(sc_ion.theta)) / (ion->A + tatom->A));
            recoil->E = ion->E * K;
            recoil->w = 0.0;
        }
        recoil->theta = sc_target.theta;
        recoil->fii = sc_target.fii;
        recoil->opt.cos_theta = cos(recoil->theta);
        recoil->opt.sin_theta = sin(recoil->theta);
        recoil->p = ion->p;
        recoil->tlayer = ion->tlayer;
        recoil->nsct = 0;
        recoil->time = 0.0;
        recoil->type = SECONDARY;
        recoil->status = NOT_FINISHED;
        recoil->lab.p = ion->lab.p;
        recoil->lab.theta = ion->lab.theta;
        recoil->lab.fii = ion->lab.fii;
        recoil->virtual = FALSE;
        for (i = 0; i < (target->nlayers - target->ntarget); i++)
            recoil->Ed[i] = 0.0;
#ifdef REC_ANG
        recoil->rec_ion_theta = sc_ion.theta/C_DEG;
        recoil->rec_lab_theta = sc_lab.theta/C_DEG;
        if(sc_ion.fii > PI)
           sc_ion.fii -= (2.0*PI);
        if(sc_lab.fii > PI)
           sc_lab.fii -= (2.0*PI);
        recoil->rec_ion_fii = sc_ion.fii/C_DEG;
        recoil->rec_lab_fii = sc_lab.fii/C_DEG;
#endif

        save_ion_history(global, ion, recoil, detector, sc_lab, sc_ion,
                         tatom->Z, tatom->A);

#ifdef DEBUG
        print_ion_position(global,ion,"R",ANYSIMULATION);
#endif
        return (TRUE);
    } else {
        fatal_error("Unknown simulation type\n");
    }
    return (FALSE);
}

double get_isotope(Isotopes *I) {
    double r, c = 0.0, mass;
    int i = 0;

    r = rnd(0.0, I->c_sum, RND_OPEN, RND_CONT);
    assert(I->c_sum > 0.9999 && I->c_sum < 1.0001);

    while (i < I->n && r > c) {
        c += I->c[i];
        i++;
    }

    if (i == I->n && r > c) {
        fprintf(stderr,
                "Error calculating random isotope mass. Got isotope %i out of %i, which doesn't make sense.\nOur winning lottery ticket number is: %g, which is greater than %g.\n",
                i + 1, I->n, r, c);
        double sum = 0.0;
        for (i = 0; i < I->n; i++) {
            fprintf(stderr, "Isotope %i: mass %g u, concentration: %g%%\n", i + 1, I->A[i] / C_U, I->c[i] * 100.0);
            sum += I->c[i];
        }
        fprintf(stderr, "Total concentration %g%%.", sum * 100.0);
        exit(10);
    }
    mass = I->A[i - 1];

    return (mass);
}

void save_ion_history(Global *global, Ion *ion, Ion *recoil, Detector *detector,
                      Vector sc_lab, Vector sc_ion, double Z, double A) {
    recoil->hist.tar_recoil.p = recoil->p;
    recoil->hist.tar_recoil.theta = recoil->theta;
    recoil->hist.tar_recoil.fii = recoil->fii;

    recoil->hist.lab_recoil = sc_lab;
    recoil->hist.lab_recoil.p = ion->lab.p;

    recoil->hist.ion_recoil = sc_ion;

    recoil->hist.tar_primary.p = ion->p;
    recoil->hist.tar_primary.theta = ion->theta;
    recoil->hist.tar_primary.fii = ion->fii;

    recoil->hist.lab_primary = ion->lab;

    recoil->hist.ion_E = ion->E;
    recoil->hist.recoil_E = recoil->E;

    recoil->hist.layer = recoil->tlayer;
    recoil->hist.nsct = ion->nsct;
    recoil->hist.time = ion->time;

    recoil->hist.w = recoil->w;

    recoil->hist.Z = Z;  /* Atomic number of the target particle */
    recoil->hist.A = A;  /* Mass of the target particle */
}

Vector get_recoiling_dir(Global *global, Ion *ion, Ion *recoil, Target *target,
                         Detector *det) {
/*
 *  We calculate here the recoiling angle in the laboratory coordinates
 *  (narrow scheme) or in the ion coordinates (wide scheme) 
 */

    Vector d, d_ion, d_target;
    double t, costhetamax, thetamax, theta, fii;
    double w, min_ion_theta, min_theta, min_fii;
    int n, i = 0, j, nimportance;

    if (global->simstage == PRESIMULATION) {
        theta = 0.0;
        fii = 0.0;
        rotate(det->angle, 0.0, theta, fii, &(d.theta), &(d.fii));
    } else if (global->recwidth == REC_NARROW) {
        n = ion->tlayer;
        thetamax = target->recpar[n].x * ion->p.z + target->recpar[n].y;
/*
      thetamax = max(0.0,thetamax);
      thetamax += det->thetamax;
*/
        costhetamax = cos(thetamax);

        /*
         *   This may be dangerous if thetamax is very small
         *   because costhetamax is then very close to 1.0.
         *   This is why we decrease costhetamax if loop
         *   is gone through for more than 5 times.
         */


        do {
            i++;
            t = rnd(costhetamax, 1.0, RND_CLOSED, RND_CONT);
            theta = acos(t);
            fii = rnd(0, 2.0 * C_PI, RND_OPEN, RND_CONT);
            rotate(det->angle, 0.0, theta, fii, &(d.theta), &(d.fii));
            rotate(global->beamangle, 0.0, d.theta, d.fii,
                   &(d_target.theta), &(d_target.fii));
            rotate(ion->theta, ion->fii - C_PI, d_target.theta, d_target.fii,
                   &(d_ion.theta), &(d_ion.fii));
            if (i % 5 == 0) {
                thetamax = min(thetamax * 2.0, C_PI / 2.0);
                costhetamax = cos(thetamax);
                printf("recoil: %6.2f %6.2f %6.2f %6.2f\n", d_ion.theta / C_DEG,
                       d_ion.fii / C_DEG, thetamax / C_DEG, acos(global->costhetamin) / C_DEG);
                printf("%i calculations for recoil angle, increasing thetamax\n", i);
            }
            if (i > 100) {
                printf("%i calculations for recoil angle, exiting...\n", i);
                exit(11);
            }
        } while (cos(d_ion.theta) > global->costhetamin);
        /*
        if(i>1)
           printf("Scattering angle recalculated\n");
        */
/*
      printf("T %10.5f %10.5f\n",d_ion.theta/C_DEG,ipow2(thetamax/det->thetamax));
      printf("P %10.5f %10.5f\n",theta/C_DEG,ipow2(thetamax/det->thetamax));
*/
        rotate(det->angle, 0.0, theta, fii, &(d.theta), &(d.fii));
        recoil->w *= ipow2(thetamax / det->thetamax);

    } else { /* Wide simulation scheme, direction in ion coordinates */
        t = rnd(global->costhetamax, global->costhetamin, RND_CLOSED, RND_CONT);
        d.theta = acos(t);
        d.fii = rnd(0, 2.0 * C_PI, RND_OPEN, RND_CONT);
    }

    return (d);
}

double cross_section(double z1, double m1, double z2, double m2, double E,
                     double theta, double cross[][NSANGLE], int table,
                     int ctype) {
    /*
     *   This function calculates the cross section for recoiling or scattering,
     *   either Rutherford values or interpolates from a given table and returns
     *   the value in units of b/sr.
     */

    double value, tcm[2], Ecm;
    int n;

    if (ctype == SIM_ERD) {
        if (theta >= (C_PI / 2.0))
            value = -1.0;
        else
            value = Serd(z1, m1, z2, m2, theta, E);
    } else { /* SIM_RBS */
        n = cm_angle(theta, m1 / m2, tcm);
        Ecm = m2 * E / (m1 + m2);
        if (n == 0)
            value = -1.0;
        else
            value = mc2lab_scatc(Srbs_mc(z1, z2, tcm[0], Ecm), tcm[0], theta);
    }

    if (table) {
        value *= inter_cross_table(E, theta, cross);
    }

    return (value / C_BARN);

}

double Serd(double z1, double m1, double z2, double m2, double t, double E) {
    double value;

    value = ipow2(z1 * z2 * C_E * C_E / (8.0 * C_PI * C_EPSILON0 * E)) * ipow2(1.0 + m1 / m2) /
            ipow(cos(t), 3);

    return (value);
}

double Srbs_mc(double z1, double z2, double t, double E) {
    double value;

    value = ipow2((z1 * z2 * C_E * C_E) / (4.0 * C_PI * C_EPSILON0)) * ipow2(1.0 / (4.0 * E)) *
            ipow(1.0 / sin(t / 2.0), 4);

    return (value);
}

double mc2lab_scatc(double mcs, double tcm, double t) {
    double value;

    value = (mcs * ipow2(sin(tcm))) / (ipow2(sin(t)) * cos(tcm - t));

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

double inter_cross_table(double E, double theta, double cross[][NSANGLE]) {
    int i, j;

    i = j = 1;
/*   
   if(E < cross[1][0])
      i = 1;
   else {
   } 
*/

    return (1.0);
}
