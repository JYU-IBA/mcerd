#include "general.h"
#include "prototypes.h"


int is_in_energy_detector(Global *global, Ion *ion, Target *target, Detector *detector, int beoynd_detector_ok) {
    if (detector->edet[0] <= target->ntarget) /* Probably not set or otherwise invalid. */
        return FALSE;
    if (ion->tlayer == detector->edet[0]) { /* We are actually in the (first layer of) energy detector */
        return TRUE;
    }
    if (ion->tlayer >= detector->edet[0] && beoynd_detector_ok) {
        return TRUE; /* If beyond_detector_ok is TRUE, we will return true also beyond the detector */
    }
    return FALSE;

}


void move_to_erd_detector(Global *global, Ion *ion, Target *target, Detector *detector) {
/*   
 *    Here the recoil moves from the target to a detector foil or 
 *    from one detector foil to the next.
 */

    Line ion_line;
    Det_foil *foil;
    Point pout, cross, fcross;
    Vector dir;
    double out_theta, out_fii, theta, fii, dist;
    int is_cross, n, i;

    if (ion->tlayer < target->ntarget) {    /* Just recoiled from the target */
        ion->tlayer = target->ntarget;
        for (i = 0; i < detector->nfoils; i++) {
            ion->Ed[i] = 0.0;
            ion->dt[i] = 0.0;
            ion->hit[i].x = 0.0;
            ion->hit[i].y = 0.0;
        }
    }

    if (ion->tlayer >= target->nlayers) {
        ion->status = FIN_DET;
#ifdef DEBUG
        fprintf(global->master.fpdebug, "T It's the final (detector) layer!\n");
#endif
        return;
    }

    n = ion->tlayer - target->ntarget; /* n = foil number */
#ifdef DEBUG
    fprintf(global->master.fpdebug, "T We're in detector layer %i, tlayer=%i.\n", n, ion->tlayer);
#endif
    ion->Ed[n] = ion->E;

    if (n < 0 || n >= detector->nfoils)
        fatal_error("Ion in a strange detector foil\n");

    foil = &(detector->foil[n]);

    pout = coord_transform(ion->lab.p, ion->lab.theta, ion->lab.fii, ion->p, FORW);

    rotate(ion->lab.theta, ion->lab.fii, ion->theta, ion->fii, &out_theta, &out_fii);

    dir.theta = out_theta;
    dir.fii = out_fii;

    ion_line = get_line_params(&dir, &pout);

    is_cross = get_line_plane_cross(&ion_line, &(foil->plane), &cross);
#ifdef DEBUG
    fprintf(global->master.fpdebug, "T is_cross\n");
#endif
    if (is_cross) {
        ion->lab.p = cross;
        ion->p.x = 0.0;
        ion->p.y = 0.0;
        ion->p.z = 0.0;
        theta = detector->angle;
        fii = C_PI;
        fcross = coord_transform(foil->center, theta, fii, cross, BACK);
        ion->hit[n] = fcross;
    }
#ifdef  DEBUG
    print_ion_position(global,ion,"F",ANYSIMULATION);
#endif
/*
   printf("POS %10.3f %10.3f\n",ion->hist.tar_recoil.p.z/C_NM,
          get_distance(&cross,&(detector->vfoil.center))/C_MM);
*/
    if (is_cross && is_on_right_side(&pout, &dir, &cross) &&
        is_in_foil(&fcross, foil)) {
#ifdef DEBUG
        fprintf(global->master.fpdebug, "T is_cross and is_on_right_size and is in foil\n");
#endif
        dist = get_distance(&pout, &cross);
        ion->time += dist / sqrt(2.0 * ion->E / ion->A);
        ion->dt[n] = ion->time;
        ion->status = NOT_FINISHED;
        ion->lab.theta = detector->angle;
        ion->lab.fii = 0.0;
        rotate(detector->angle, C_PI, out_theta, out_fii, &(ion->theta), &(ion->fii));
        ion->opt.cos_theta = cos(ion->theta);
        ion->opt.sin_theta = sin(ion->theta);
    } else {
        if (n == 0 && is_cross && global->virtualdet &&
            is_on_right_side(&pout, &dir, &cross) &&
            is_in_foil(&fcross, &(detector->vfoil))) {
            ion->status = NOT_FINISHED;
#ifdef DEBUG
            print_ion_position(global,ion,"W",ANYSIMULATION);
#endif
            hit_virtual_detector(global, ion, target, detector, &cross, &pout);
            ion->dt[n] = ion->time;
            ion->virtual = TRUE;
#ifdef DEBUG
            fprintf(global->master.fpdebug, "T VIRTUAL!\n");
#endif
        } else {
            ion->lab.theta = detector->angle;
            ion->lab.fii = 0.0;
            if (n > 0) {
                ion->status = FIN_OUT_DET;
#ifdef DEBUG
                fprintf(global->master.fpdebug, "T Out of detector, n=%i, tlayer=%i, T1=%g T2=%g\n", n, ion->tlayer, ion->dt[0]/C_NS, ion->dt[1]/C_NS);
#endif
            } else {
                ion->status = FIN_MISS_DET;
            }
        }
    }

#ifdef DEBUG
    if(ion->virtual)
       print_ion_position(global,ion,"V",ANYSIMULATION);
#endif

}

Line get_line_params(Vector *d, Point *p) {
/*
 *    This calculates the line parameters for a line pointing to the
 *    direction d and going through point p.
 *
 */


    Line k;
    double x, y, z;

    z = cos(d->theta);
    y = sin(d->theta) * sin(d->fii);
    x = sin(d->theta) * cos(d->fii);

    if (z == 0.0)
        if (x == 0.0) {
            k.type = L_YAXIS;
            k.a = p->x;
            k.b = p->z;
        } else {
            k.type = L_XYPLANE;
            k.a = y / x;
            k.b = p->y - k.a * p->x;
            k.c = p->z;
        }
    else {
        k.type = L_GENERAL;
        k.a = x / z;
        k.b = p->x - k.a * p->z;
        k.c = y / z;
        k.d = p->y - k.c * p->z;
    }

    return (k);

}

int get_line_plane_cross(Line *line, Plane *plane, Point *cross) {

    Plane p;
    Line q;
    Point k;
    int c = NO_CROSS;

    p = *plane;
    q = *line;

    switch (p.type) {
        case GENERAL_PLANE:
            switch (q.type) {
                case L_GENERAL:
                    if ((q.a - p.b - p.a * q.c) != 0.0) {
                        c = CROSS;
                        k.z = (p.a * q.b + p.b * q.d + p.c) / (1.0 - p.a * q.a - p.b * q.c);
                        k.x = q.a * k.z + q.b;
                        k.y = q.c * k.z + q.d;
                    } else
                        c = NO_CROSS;
                    break;
                case L_XYPLANE:
                    c = CROSS;
                    k.x = -q.b / q.a - (-p.a * q.b - q.a * (-p.c - p.b * q.c)) /
                                       (q.a * (p.a - q.a));
                    k.y = -(-p.a * q.b - q.a * (-p.c - p.b * q.c)) / (p.a - q.a);
                    k.z = q.c;
                    break;
                case L_YAXIS:
                    k.x = q.a;
                    k.y = p.a * q.a + p.b * q.b + p.c;
                    k.z = q.b;
                    c = CROSS;
                    break;
            }
            break;
        case Z_PLANE:
            switch (q.type) {
                case L_GENERAL:
                    c = CROSS;
                    k.x = q.c * p.a + q.d;
                    k.y = q.a * p.a + q.b;
                    k.z = p.a;
                    break;
                case L_XYPLANE:
                    c = NO_CROSS;
                    break;
                case L_YAXIS:
                    c = NO_CROSS;
                    break;
            }
            break;
        case X_PLANE:
            switch (q.type) {
                case L_GENERAL:
                    c = CROSS;
                    k.x = p.a;
                    k.z = (p.a - q.d) / (q.c + 1.0e-20);
                    k.y = k.z * q.a + q.b;
                    break;
                case L_XYPLANE:
                    c = CROSS;
                    k.x = p.a;
                    k.y = p.a * q.a + q.b;
                    k.z = q.c;
                    break;
                case L_YAXIS:
                    c = NO_CROSS;
                    break;
            }
            break;
    }

    *cross = k;

    return (c);

}

double get_distance(Point *p1, Point *p2) {
    double r;

    r = sqrt(ipow2(p1->x - p2->x) + ipow2(p1->y - p2->y) + ipow2(p1->z - p2->z));

    return (r);
}

int is_on_right_side(Point *p1, Vector *d, Point *p2) {
    Point dp;
    double r1, r2;

    r1 = get_distance(p1, p2);

    dp.x = p1->x + r1 * sin(d->theta) * cos(d->fii);
    dp.y = p1->y + r1 * sin(d->theta) * sin(d->fii);
    dp.z = p1->z + r1 * cos(d->theta);

    r2 = get_distance(&dp, p2);

    if (r2 < r1)
        return (TRUE);
    else
        return (FALSE);

}

int is_in_foil(Point *p, Det_foil *foil) {
    /* Here we assume that point p is on the plane of circle c */

    switch (foil->type) {
        case FOIL_CIRC:
            if (foil->virtual) {
                if (sqrt(ipow2(p->x / foil->size[0]) + ipow2(p->y / foil->size[1])) <= 1.0)
                    return (TRUE);
                else
                    return (FALSE);
            } else {
                if (sqrt(ipow2(p->x) + ipow2(p->y)) <= foil->size[0])
                    return (TRUE);
                else
                    return (FALSE);
            }
            break;
        case FOIL_RECT:
            if (fabs(p->x) <= foil->size[0] && fabs(p->y) <= foil->size[1])
                return (TRUE);
            else
                return (FALSE);
            break;
        default:
            break;
    }
    return (FALSE);
}

double get_angle_between_points(Point *p, Point *p1, Point *p2) {
    Point dp1, dp2;
    double r1, r2, angle;

    dp1.x = p1->x - p->x;
    dp1.y = p1->y - p->y;
    dp1.z = p1->z - p->z;

    dp2.x = p2->x - p->x;
    dp2.y = p2->y - p->y;
    dp2.z = p2->z - p->z;

    r1 = sqrt(ipow2(dp1.x) + ipow2(dp1.y) + ipow2(dp1.z));
    r2 = sqrt(ipow2(dp2.x) + ipow2(dp2.y) + ipow2(dp2.z));

    angle = acos((dp1.x * dp2.x + dp1.y * dp2.y + dp1.z * dp2.z) / (r1 * r2));

    return (angle);
}

#if 0
int is_in_rect(Point *p,Rect *r)
{
   Point pab,pbc;    /* These are projections of point p to the */
                     /* lines between rectangle corners a-b and b-c */
   double d_ab,d_bc; /* Distances between a-b and a-c */
   double d,d1,d2,theta;

   d1 = get_distance(&(r->p[0]),p);
   d2 = get_distance(&(r->p[1]),p);   

   theta = get_angle_between_points(&(r->p[0]),&(r->p[1]),p);

   d_ab = get_distance(&(r->p[0]),&(r->p[1]));
   
   d = d1/d_ab;
   d *= cos(theta);

   pab.x = r->p[0].x + (r->p[1].x - r->p[0].x)*d;
   pab.y = r->p[0].y + (r->p[1].y - r->p[0].y)*d;
   pab.z = r->p[0].z + (r->p[1].z - r->p[0].z)*d;

   theta = get_angle_between_points(&(r->p[1]),&(r->p[2]),p);
   
   d_bc = get_distance(&(r->p[1]),&(r->p[2]));
   
   d = d2/d_bc;
   d *= cos(theta);

   pbc.x = r->p[1].x + (r->p[2].x - r->p[1].x)*d;
   pbc.y = r->p[1].y + (r->p[2].y - r->p[1].y)*d;
   pbc.z = r->p[1].z + (r->p[2].z - r->p[1].z)*d;

   if(get_distance(&pab,&(r->p[0])) < d_ab && 
      get_distance(&pab,&(r->p[1])) < d_ab &&
      get_distance(&pbc,&(r->p[1])) < d_bc && 
      get_distance(&pbc,&(r->p[2])) < d_bc)
      return(TRUE);
   else
      return(FALSE);   

}
#endif
