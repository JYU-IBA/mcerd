#include "general.h"
#include "prototypes.h"

#define MAXVRATIO 0.80  /* Maximum relative change of ion energy
                           in energy corrections for energy loss and 
                           kinematics */

double get_eloss_corr(Ion *, Target *, double, double, double);
double get_eloss(double, Ion *, Target *, double);

double rbs_kin(double, double, double, int *);
double rbs_cross(double, double, double, int *);

void hit_virtual_detector(Global *global, Ion *ion, Target *target,
                          Detector *det, Point *p_virt_lab, Point *p_out_tar) {
/*
 *    Here we calculate a random projection point from virtual detector to
 *    the real detector, and correct the ion energy according to the
 *    changed kinematics.
 */

    Vector v_real_foil, v_virt_lab, v_real_lab, v_diff_lab;
    Vector v_recvirt_pri, v_recreal_pri, v_rec_tar, v_rec_lab, v_virt_tar, v_real_tar;
    Point p_rec_tar;
    double theta, fii, dx, dy, r, dE1, dE2, dw, dist, eloss, m1, m2;
    int err = FALSE, err1, err2, n;

/*
 *    For kinematic correction we calculate first the angular change we have 
 *    to make in order to move from the virtual hit point to the real hit
 *    point and then apply this to the original recoiling direction to 
 *    calculate the new initial recoiling direction. Since the original 
 *    recoiling direction can be almost anything (especially in the wide 
 *    mode) it is possible that we end up to a unphysical new recoiling 
 *    direction. Same angles than for the kinematic correction are used for
 *    the cross section correction.
 *
 *    For stopping correction the energy loss and direct path length in the 
 *    target material is calculated for the initial ion path hitting to the 
 *    virtual detector and the energy loss is corrected for the change in 
 *    the path length to the real hit point.
 * 
 */

/*
 *    Variable naming convention: first part distinguishes between vectors
 *    (v) and points (p), second part (after p or v) gives the direction 
 *    from the recoiling point [rec = recoil direction, virt = direction to the 
 *    virtual hit point, real = direction to the real (new) hit point].
 *    Last part gives the frame of reference [foil = foil coordinate system,
 *    pri = primary ion (at the moment of recoiling) coordinate system,
 *    tar = target coordinate system, lab = laboratory coordinate system].
 *
 *    *p_virt_lab:   hit point in the virtual foil in the lab. coordinates
 *    v_real_foil:   new hit position in the real foil coordinates
 *    v_recvirt_pri: recoil vector in primary ion coord. at the recoil moment
 *    v_rec_tar:     recoil vector in target ion coordinates at the recoil moment
 *    v_rec_lab:     recoil vector in lab. coordinates at the recoil moment
 *    v_virt_lab:    vector from the recoil point to the virtual hit point
 *    v_real_lab:    vector from the recoil point to the new real hit point
 *    v_diff_lab:    angle between old and new hit point from recoil point
 *    v_recreal_pri: recoil vector to the new real hit point in pri. ion coord.
 *
 */

    m1 = m2 = dE1 = dw = 0.0;

    switch (det->vfoil.type) {
        case FOIL_CIRC:
            r = det->foil[0].size[0] * sqrt(rnd(0.0, 1.0, RND_CLOSED, RND_CONT));
            fii = rnd(0.0, 2.0 * C_PI, RND_CLOSED, RND_CONT);
            v_real_foil.p.x = r * cos(fii);
            v_real_foil.p.y = r * sin(fii);
            break;
        case FOIL_RECT:
            dx = det->foil[0].size[0];
            v_real_foil.p.x = rnd(-dx, dx, RND_CLOSED, RND_CONT);
            dy = det->foil[0].size[1];
            v_real_foil.p.y = rnd(-dy, dy, RND_CLOSED, RND_CONT);
            break;
    }
    v_real_foil.p.z = 0.0;
    theta = det->angle;
    fii = 0.0;

    v_recvirt_pri = ion->hist.ion_recoil;
    v_rec_tar = ion->hist.tar_recoil;
    v_rec_lab = ion->hist.lab_recoil;

/*  
 *  Because of the beam spot effect v_virt_lab and v_real_lab must be given
 *  in the coordinates relative to the recoil point. Directions are still in
 *  the laboratory coordinates.
 */

    v_virt_lab.p = *p_virt_lab;

    v_virt_lab.p.x -= v_rec_lab.p.x;
    v_virt_lab.p.y -= v_rec_lab.p.y;
    v_virt_lab.p.z -= v_rec_lab.p.z;

    r = sqrt(ipow2(v_virt_lab.p.x) + ipow2(v_virt_lab.p.y) + ipow2(v_virt_lab.p.z));
    v_virt_lab.theta = acos(v_virt_lab.p.z / r);
    v_virt_lab.fii = atan2(v_virt_lab.p.y, v_virt_lab.p.x);

    v_real_lab.p = coord_transform(det->vfoil.center, theta, fii, v_real_foil.p, FORW);

    v_real_lab.p.x -= v_rec_lab.p.x;
    v_real_lab.p.y -= v_rec_lab.p.y;
    v_real_lab.p.z -= v_rec_lab.p.z;

    r = sqrt(ipow2(v_real_lab.p.x) + ipow2(v_real_lab.p.y) + ipow2(v_real_lab.p.z));

    v_real_lab.theta = acos(v_real_lab.p.z / r);
    v_real_lab.fii = atan2(v_real_lab.p.y, v_real_lab.p.x);

    rotate(v_virt_lab.theta, C_PI + v_virt_lab.fii, v_real_lab.theta, v_real_lab.fii,
           &(v_diff_lab.theta), &(v_diff_lab.fii));

    rotate(v_recvirt_pri.theta, v_recvirt_pri.fii, v_diff_lab.theta, v_diff_lab.fii,
           &(v_recreal_pri.theta), &(v_recreal_pri.fii));

    if (global->simtype == SIM_RBS) {
        m1 = ion[PRIMARY].A;
        m2 = ion[SECONDARY].A;   /* mass of the scattering target atom */

        /* Ion actually points here to the SECONDARY ie. scattered ion, thus
           ion[SECONDARY] means the target atom. */
    }

#ifdef VIR_DEBUG
    printf("KIN: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           v_recvirt_pri.theta/C_DEG,v_recvirt_pri.fii/C_DEG,
           v_diff_lab.theta/C_DEG,v_diff_lab.fii/C_DEG,
           v_recreal_pri.theta/C_DEG,v_recreal_pri.fii/C_DEG);
#endif

/*
 *    dE1: relative change in ion energy in the kinematics correction
 *    dE2: energy loss correction
 */

    if (global->simtype == SIM_ERD)
        dE1 = (ipow2(cos(v_recreal_pri.theta) / cos(v_recvirt_pri.theta)) - 1.0);
    else if (global->simtype == SIM_RBS) {
        dE1 = rbs_kin(m1, m2, v_recreal_pri.theta, &err1) /
              rbs_kin(m1, m2, v_recvirt_pri.theta, &err2) - 1.0;
        if (err1 || err2)
            err = TRUE;
    }

/*
 *    p_rec_tar:   recoiling point in target coordinates
 *    v_virt_tar: old (virtual) recoiling direction in target coordinates
 *    v_real_tar: new (real) recoiling direction in target coordinates
 *
 *    v_rec_tar.p is in the coordinates of the specific ion
 *    p_rec_tar is in the coordinates of the target coordinate system
 */

    p_rec_tar = coord_transform(v_rec_lab.p, global->beamangle, C_PI, v_rec_tar.p, FORW);

    v_virt_tar.p.x = p_out_tar->x - p_rec_tar.x;
    v_virt_tar.p.y = p_out_tar->y - p_rec_tar.y;
    v_virt_tar.p.z = p_out_tar->z - p_rec_tar.z;

    r = sqrt(ipow2(v_virt_tar.p.x) + ipow2(v_virt_tar.p.y) + ipow2(v_virt_tar.p.z));

    v_virt_tar.theta = acos(v_virt_tar.p.z / r);
    v_virt_tar.fii = atan2(v_virt_tar.p.y, v_virt_tar.p.x);

/*
 *    Here v_virt_tar direction is in laboratory coordinates. We now calculate
 *    the directions in target coordinates for both virtual and real.
 */

    rotate(v_virt_tar.theta, v_virt_tar.fii, v_diff_lab.theta, v_diff_lab.fii,
           &(v_real_tar.theta), &(v_real_tar.fii));

    rotate(global->beamangle, 0.0, v_virt_tar.theta, v_virt_tar.fii,
           &(v_virt_tar.theta), &(v_virt_tar.fii));

    rotate(global->beamangle, 0.0, v_real_tar.theta, v_real_tar.fii,
           &(v_real_tar.theta), &(v_real_tar.fii));

#ifdef VIR_DEBUG
    printf("TAR: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           v_virt_tar.theta/C_DEG,v_virt_tar.fii/C_DEG,
           v_virt_lab.theta/C_DEG,v_virt_lab.fii/C_DEG,
           v_real_tar.theta/C_DEG,v_real_tar.fii/C_DEG,
           v_recreal_pri.theta/C_DEG,v_recreal_pri.fii/C_DEG);
#endif

    eloss = ion->hist.recoil_E - ion->E;

    dE2 = get_eloss_corr(ion, target, dE1, v_virt_tar.theta, v_real_tar.theta);

/*
#ifdef VIR_DEBUG
   printf("DIST %10.3f %10.3f %10.3f %10.3f\n",d1/C_NM,d2/C_NM,eloss/C_MEV,dE2);
#endif
*/
    dE1 *= ion->hist.recoil_E;
    dE2 *= -eloss;

    if (fabs(dE1 + dE2) / ion->E > MAXVRATIO) {
        ion->status = FIN_MISS_DET;
        fprintf(stderr, "Energy would change too much in virtual detector ");
        fprintf(stderr, "%10.3f MeV\n", (dE1 + dE2) / C_MEV);
    }

/*
 *    dw is the correction factor cross section when we move from the virtual
 *    detector to the real detector
 */

    if (global->simtype == SIM_ERD)
        dw = ipow(cos(v_recvirt_pri.theta) / cos(v_recreal_pri.theta), 3);
    else if (global->simtype == SIM_RBS) {
        dw = rbs_cross(m1, m2, v_recreal_pri.theta, &err1) /
             rbs_cross(m1, m2, v_recvirt_pri.theta, &err2);
        if (err1 || err2)
            err = TRUE;
    }

    n = ion->tlayer - target->ntarget; /* n = foil number */

#ifdef VIR_DEBUG
    printf("VIR1 %7.3f %8.3f %7.3f %8.3f %7.3f %8.3f %7.2f %7.2f %6.3f\n",
           v_virt_lab.theta/C_DEG,v_virt_lab.fii/C_DEG,
           v_real_lab.theta/C_DEG,v_real_lab.fii/C_DEG,
           ion->hist.recoil_E/C_MEV,
           ion->hit[n].x/C_MM,ion->hit[n].y/C_MM,
           dE1/C_MEV,dE2/C_MEV);
    printf("VIR2 %7.3f %8.3f %7.3f %8.3f %7.3f %8.3f %7.2f %7.2f %6.3f\n",
           v_recvirt_pri.theta/C_DEG,v_recvirt_pri.fii/C_DEG,
           v_recreal_pri.theta/C_DEG,v_recreal_pri.fii/C_DEG,
           ion->hist.recoil_E/C_MEV,
           ion->hit[n].x/C_MM,ion->hit[n].y/C_MM,
           dE1/C_MEV,dE2/C_MEV);
#endif

    ion->w *= dw;
    ion->E += dE1 + dE2;
    ion->Ed[n] = ion->E;

    /* We may end up with a negative energy */

    if (ion->E > global->ionemax || ion->E < global->emin) {
        ion->status = FIN_MISS_DET;
        fprintf(stderr, "Energy would change too much in virtual detector\n");
    }

    if (err) {
        ion->status = FIN_MISS_DET;
        fprintf(stderr, "Unphysical RBS-scattering in virtual detector\n");
    }

    if (ion->status == NOT_FINISHED) {

        v_real_lab.p.x += v_rec_lab.p.x;
        v_real_lab.p.y += v_rec_lab.p.y;
        v_real_lab.p.z += v_rec_lab.p.z;

        dist = get_distance(p_out_tar, &(v_real_lab.p));
        ion->time += dist / sqrt(2.0 * ion->E / ion->A);

        ion->lab.p = v_real_lab.p;
        ion->lab.theta = det->angle;
        ion->lab.fii = 0.0;

        ion->p.x = 0.0;
        ion->p.y = 0.0;
        ion->p.z = 0.0;

        rotate(det->angle, C_PI, v_real_lab.theta, v_real_lab.fii,
               &(ion->theta), &(ion->fii));

        ion->opt.cos_theta = cos(ion->theta);
        ion->opt.sin_theta = sin(ion->theta);
    }

}

double get_eloss_corr(Ion *ion, Target *target, double dE1, double theta1,
                      double theta2) {

    double eloss, Ed1, Ed2, E;

    E = ion->hist.recoil_E;

    eloss = ion->hist.recoil_E - ion->E;

/*
 *    We calculate here the effective (according to the direct trajectory) 
 *    energy loss Ed1 to the virtual detector point using initial recoil 
 *    energy E, and the effective energy loss Ed2 to the real detector point
 *    using the kinematically corrected recoil energy (E + dE1*E) in order
 *    to take into account the systematic changes in the energy loss.
 */

    Ed1 = get_eloss(E, ion, target, cos(theta1));

    Ed2 = get_eloss(E + dE1 * E, ion, target, cos(theta2));

    if (Ed1 > E || Ed2 > (E + dE1 * E))
        return (-2.0);


    if (Ed1 > 0.0)
        Ed2 = Ed2 / Ed1 - 1.0;
    else
        Ed2 = 0.0;

#ifdef VIR_DEBUG
    printf("ELO %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
           (theta1-PI/2.0)/C_DEG,(theta2-PI/2.0)/C_DEG,
           E/C_MEV,(E + dE1*E)/C_MEV,
           eloss/C_MEV,Ed1/C_MEV,Ed2/C_MEV);
#endif

    return (Ed2);

}

double get_eloss(double E, Ion *ion, Target *target, double cos_theta) {
    Target_layer *layer;
    Vector rec;
    double d, v, sto, v2, sto2, nextz, Eorig;
    int nlayer, cont = TRUE, next_layer, surf = FALSE;

    nlayer = ion->hist.layer;
    layer = &(target->layer[nlayer]);

    rec.p.z = ion->hist.tar_recoil.p.z;

    Eorig = E;

    while (cont) {
        d = 20.0 * C_NM;
        nextz = rec.p.z + d * cos_theta;
        v = sqrt(2.0 * E / ion->A);
        next_layer = FALSE;
        if (nextz < layer->dlow) {
            nextz = layer->dlow + 0.1 * C_ANGSTROM * cos_theta;
            d = (nextz - rec.p.z) / cos_theta;
            next_layer = TRUE;
        }

        sto = inter_sto(&(layer->sto[ion->scatindex]), v, STOPPING);

        if ((E - d * sto) < 0.0)
            cont = FALSE;
        else {
            v2 = sqrt(2.0 * (E - d * sto) / ion->A);
            sto2 = inter_sto(&(layer->sto[ion->scatindex]), v2, STOPPING);
            E -= 0.5 * (sto + sto2) * d;
        }

        if (E < 0.0)
            cont = FALSE;

        if (next_layer) {
            nlayer--;
            if (nlayer < 0) {
                cont = FALSE;
                surf = TRUE;
            } else
                layer = &(target->layer[nlayer]);
        }

        rec.p.z = nextz;

#ifdef VIR_DEBUG
        printf("STO %10.3f %10.3f %4i\n",E/C_MEV,rec.p.z/C_NM,nlayer);
#endif
    }

    if (surf)
        return (Eorig - E);
    else
        return (1.1 * Eorig);
}

double rbs_cross(double m1, double m2, double theta, int *err) {
    double sq_tmp, value;

    *err = FALSE;

    if (sin(theta) >= m2 / m1) {
        *err = TRUE;
        return (1.0);
    }

    if (fmod(theta, C_PI) == 0) {
        *err = TRUE;
        return (1.0);
    }

    sq_tmp = sqrt(ipow2(m2) - ipow2(m1 * sin(theta)));

    value = ipow2(sq_tmp + m2 * cos(theta)) / (m2 * pow(sin(theta), 4) * sq_tmp);

    return (value);
}

double rbs_kin(double m1, double m2, double theta, int *err) {
    double sq_tmp, value;

    *err = FALSE;

    if (sin(theta) >= m2 / m1) {
        *err = TRUE;
        return (1.0);
    }

    sq_tmp = sqrt(ipow2(m2) - ipow2(m1 * sin(theta)));

    value = ipow2((sq_tmp + m1 * cos(theta)) / (m1 + m2));

    return (value);
}
