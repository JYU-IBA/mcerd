#include "general.h"
#include "prototypes.h"
#include "rough_surface.h"
#include <assert.h>

int recdist_crossing(Global *, Ion *, Target *, double, double *);
int recdist_nonzero(Global *, Ion *, Target *, double);

void create_ion(Global *global, Ion *ion, Target *target) {
/*
 *     Create the ion ie. calculate and give all the parameters for the
 *     new ion (primary or secondary).
 *
 *     About coordinate systems: 
 *       - laboratory coordinates: z-axis is along the beam line
 *       - target: z-axis is along the target normal
 *       - ion: z-axis is along the direction of ion movement
 *       - detector: z-axis is along the detector direction
 */

    Point pt;
    double theta, fii, t, btheta, bfii;

    switch (global->simtype) {
        case SIM_ERD:
            /*
             *   Here we also give the random hit point to the target surface.
             */
            pt.x = rnd(-global->bspot.x, global->bspot.x, RND_CLOSED, RND_CONT);
            pt.y = rnd(-global->bspot.y, global->bspot.y, RND_CLOSED, RND_CONT);
            if (global->beamangle > 0)
                pt.z = pt.x / tan(PI / 2.0 - global->beamangle);
            else
                pt.z = 0.0;
            break;
        case SIM_RBS:
            /*
             *   Here we also give the random hit point to the target surface.
             */
            pt.x = rnd(-global->bspot.x, global->bspot.x, RND_CLOSED, RND_CONT);
            pt.y = rnd(-global->bspot.y, global->bspot.y, RND_CLOSED, RND_CONT);
            if (global->beamangle > 0)
                pt.z = pt.x / tan(PI / 2.0 - global->beamangle);
            else
                pt.z = 0.0;
            break;
        default:
            fatal_error("no such simulation type\n");
            break;
    }

    ion->lab.p = pt;

    if (global->simstage == REALSIMULATION &&
        global->cion % (global->nscale + 1) == 0) {
        ion->scale = TRUE;
        ion->lab.p.x = 0.0;
        ion->lab.p.y = 0.0;
        ion->lab.p.z = 0.0;
        ion->effrecd = SCALE_DEPTH;
        ion->w = global->nscale * SCALE_DEPTH / target->effrecd;
    } else {
        ion->scale = FALSE;
        ion->effrecd = target->effrecd;
        ion->w = 1.0;      /* Original ion weight is 1 */
    }

/*
   switch(global->beamprof){
      case BEAM_GAUSS:
         bfii = rnd(0,2.0*PI,RND_OPEN,RND_CONT);         
         break;
      case BEAM_FLAT:
         t = rnd(cos(global->beamdiv),1.0,RND_CLOSED,RND_CONT);
         btheta = acos(t);
         bfii = rnd(0,2.0*PI,RND_OPEN,RND_CONT);
         break;
      default:
         btheta = 0.0;
         bfii = 0.0;
         break;
   }

   rotate(global->beamangle,0.0,btheta,bfii,&theta,&fii);
*/
    theta = global->beamangle;
    fii = 0.0;

    ion->tlayer = 0;
    ion->time = 0.0;

    ion->p.x = 0.0;
#ifdef DEBUG
    ion->p.x = rnd(0.0,2000.0,RND_OPEN,RND_CONT)*C_NM;
#endif
    ion->p.y = 0.0;
    ion->p.z = 0.001 * C_ANGSTROM;
    ion->E = global->E0;
    ion->nsct = 0;
    ion->theta = theta;
    ion->fii = fii;
    ion->opt.cos_theta = cos(theta);
    ion->opt.sin_theta = sin(theta);
    ion->type = PRIMARY;
    ion->status = NOT_FINISHED;

    ion->lab.theta = global->beamangle;
    ion->lab.fii = PI;

#ifdef DEBUG
    print_ion_position(global,ion,"L",ANYSIMULATION);
       print_ion_position(global,ion,"M",ANYSIMULATION);
#endif
}

void move_target(Target *target) {
    double theta;

    target->surface.origin.x = target->surface.size * rnd(0.0, 1.0, RND_CLOSED, RND_CONT);
    target->surface.origin.y = target->surface.size * rnd(0.0, 1.0, RND_CLOSED, RND_CONT);
    theta = rnd(0.0, 2 * PI, RND_OPEN, RND_CONT);

#ifdef DEBUG
    target->surface.origin.x = 0.0;
    target->surface.origin.y = 0.0;
    theta = 0.0;
#endif

    target->surface.cos_r = cos(theta);
    target->surface.sin_r = sin(theta);
    target->surface.step = target->surface.size / (target->surface.nsize - 1);

}

void next_scattering(Global *global, Ion *ion, Target *target,
                     Scattering **scat, SNext *snext) {
/*
 *     Calculate the distance to the next scattering point and determine 
 *     the scattering atom. We loop through the atoms in the current layer.
 */

    Target_layer *layer;
    double cross = 0.0, b[MAXATOMS], rcross;
    int i, p;
    if(global->nomc) {
        snext->d = 100.0 * C_NM;
        return;
    }
    layer = &(target->layer[ion->tlayer]);

    for (i = 0; i < layer->natoms; i++) {
        p = layer->atom[i];
        b[i] = layer->N[i] * get_cross(ion, &(scat[ion->scatindex][p]));
        cross += b[i];
    }

    rcross = rnd(0.0, cross, RND_OPEN, RND_CONT);

    i = 0;
    while (i < layer->natoms && rcross >= 0.0) {
        rcross -= b[i];
        i++;
    }
    if (i < 1 || i > layer->natoms)
        fatal_error("Scattering atom calculated wrong\n");

    i--;

    snext->natom = layer->atom[i];

    ion->opt.y =
            sqrt(-rcross / (PI * layer->N[i])) / scat[ion->scatindex][snext->natom].a; /* NOTE! ion->type is wrong */
/*
   printf("%5.2f %5.2f %5.2f %7.5f %10.4f %15.7f\n",
         ion->Z,ion->A/C_U,target->ele[snext->natom].Z,
         target->ele[snext->natom].A/C_U,ion->E/C_MEV,
         sqrt((b[i]/layer->N[i])/PI)/C_ANGSTROM);
*/
    snext->d = -log(rnd(0.0, 1.0, RND_RIGHT, RND_CONT)) / (cross);
#ifdef NO_SUBNM_STEPS_WARNING_UNPHYSICAL
    if(snext->d < 1.0 * C_NM)
        snext->d = 1.0*C_NM;
    }
#endif
}

int move_ion(Global *global, Ion *ion, Target *target, SNext *snext) {
/*
 *     Move ion to the next point, which is the closest of following
 *      i) next MC-scattering point, distance already calculated
 *     ii) next crossing of the target layer, distance calculated
 *         using surface_crossing -subroutine
 *    iii) next ERD scattering, in case of a primary ion, random distance 
           to the scattering point is calculated here according
 *         to the simulation mode (wide or narrow), which is accepted only
 *         if recoil material distribu<tion value is nonzero in that point 
 *     iv) point where the relative electronic stopping power has
 *         changed more than MAXELOSS
 *
 *     To make things fast we do not evaluate all these values at the
 *     beginning but with a stepwise procedure
 *      (since we are in the innermost loop we have to worry about the speed)
 */

    Target_layer *layer;
    Point nextp, p1cross, p2cross;
    double stopping, nextz, vel, vel1, vel2, sto1, sto2, dE, drec, straggling, eloss;
    double recang, d, d1, d2, dreclayer;
    int cross_layer = FALSE, cross_recdist = FALSE, sto_dec = FALSE, sc = MC_SCATTERING;
    int n, scross, s1cross, s2cross, nz;

    d = snext->d;
    layer = &(target->layer[ion->tlayer]);

    nextp.x = d * ion->opt.sin_theta * cos(ion->fii) + ion->p.x;
    nextp.y = d * ion->opt.sin_theta * sin(ion->fii) + ion->p.y;
    nextp.z = d * ion->opt.cos_theta + ion->p.z;

    nextz = nextp.z;

/*
 *  In case we have a rough sample surface we have different handling for
 *  target and detector layer (later we may have roughness also in the 
 *  detector). With rough surfaces we can't say which surface we are going to
 *  cross first!
 */

    if (global->rough && ion->tlayer < target->ntarget) {
        scross = s1cross = s2cross = FALSE;
        if (ion->tlayer == 0) { /* We are in the virtual surface layer */
            if (nextz < layer->dlow) {
                cross_layer = -1;
                d = fabs((layer->dlow - ion->p.z) / ion->opt.cos_theta);
                scross = TRUE;
            }
        } else {
            target->surface.move = layer->dlow;
            s1cross = surface_crossing(&(target->surface), &(ion->p), &(nextp), &p1cross);
        }
        target->surface.move = layer->dhigh;
        s2cross = surface_crossing(&(target->surface), &(ion->p), &(nextp), &p2cross);
        if (s1cross) {
            scross = s1cross;
            cross_layer = -1;
            d1 = get_distance(&(ion->p), &p1cross);
            d = d1;
        }
        if (s2cross) {
            scross = s2cross;
            cross_layer = 1;
            d2 = get_distance(&(ion->p), &p2cross);
            d = d2;
        }
        if (s1cross && s2cross) {
            if (d1 < d2) {
                cross_layer = -1;
                d = d1;
            } else {
                cross_layer = 1;
                d = d2;
            }
        }
        if (scross)
            sc = NO_SCATTERING;
    } else {
        if (nextz < layer->dlow || nextz > layer->dhigh) {
            sc = NO_SCATTERING;
            if (nextz < layer->dlow) {
                d = fabs((layer->dlow - ion->p.z) / ion->opt.cos_theta);
                cross_layer = -1;
            } else {
                d = fabs((ion->p.z - layer->dhigh) / ion->opt.cos_theta);
                cross_layer = 1;
            }
        }
    }

    if (ion->type == PRIMARY && ion->tlayer > 0 && ion->tlayer < target->ntarget) {
        cross_recdist = recdist_crossing(global, ion, target, d, &dreclayer);
        if (cross_recdist) {
            cross_layer = FALSE;
            d = dreclayer;
            sc = NO_SCATTERING;
        }
        drec = -log(rnd(0.0, 1.0, RND_RIGHT, RND_CONT));
        if (global->recwidth == REC_WIDE || global->simstage == PRESIMULATION) {
            drec *= (ion->effrecd / (cos(global->beamangle) * global->nrecave));
            ion->wtmp = -1.0;
        } else {
            n = ion->tlayer;
            recang = ipow2(target->recpar[n].x * ion->p.z + target->recpar[n].y);
            ion->wtmp = target->angave / recang;
            drec *= (ion->effrecd / (cos(global->beamangle) * global->nrecave)) /
                    (recang / target->angave);
        }
        if (drec < d) {
            nz = recdist_nonzero(global, ion, target, drec);
            if (nz) { /* We have a recoiling event */
                d = drec;
                sc = ERD_SCATTERING;
                cross_layer = FALSE;
                cross_recdist = FALSE;
            }
        }
    }

    vel1 = sqrt(2.0 * ion->E / ion->A);
    sto1 = inter_sto(&(layer->sto[ion->scatindex]), vel1, STOPPING);

    dE = sto1 * d;

    vel2 = sqrt(2.0 * (max(ion->E - dE, 0.0)) / ion->A);
    sto2 = inter_sto(&(layer->sto[ion->scatindex]), vel2, STOPPING);

    while (fabs(sto1 - sto2) / sto1 > MAXELOSS || dE >= ion->E) {
        sc = NO_SCATTERING;
        cross_layer = FALSE;
        sto_dec = TRUE;
        d /= 2.0;
        dE = sto1 * d;
        vel2 = sqrt(2.0 * (max(ion->E - dE, 0.0)) / ion->A);
        sto2 = inter_sto(&(layer->sto[ion->scatindex]), vel2, STOPPING);
    }

    stopping = 0.5 * (sto1 + sto2);
    vel = 0.5 * (vel1 + vel2);


    if (sc == NO_SCATTERING)
        d += 0.01 * C_ANGSTROM; /* To be sure that the layer surface is crossed */

    if (cross_layer) {
        (ion->tlayer) += cross_layer;
    }

/*
 *  Now we know where to proceed, let's calculate straggling too and proceed
 *  with the ion.
 *
 */

    straggling = sqrt(inter_sto(&(layer->sto[ion->type]), vel, STRAGGLING) * d);

    /* we always want a positive energy loss, ie. we don't increase ion energy */


    straggling *= gaussian(RND_CONT);

#ifdef NO_STRAGGLING
    straggling = 0.0;
#endif

    eloss = d * stopping;

    if (fabs(straggling) < eloss) {  /* Here we actually forget the straggling */
        eloss += straggling;       /* if it is larger than the energy loss */
        /*
     printf("Straggling within limits %.4f keV for %.4f nm distance appliend\n",straggling/C_KEV,d/C_NM);
        */
    } else {
        if (!(ion->scale))
            /*
            printf("Maximum straggling %.4f keV for %.4f nm distance appliend\n",
               (eloss*straggling/fabs(straggling))/C_KEV,d/C_NM);
            */
            eloss += eloss * straggling / fabs(straggling);
    } /* maximum straggling applied */


#ifdef NO_ELOSS
    eloss = 0.0;
#endif

    if (ion->scale)    /* No energy loss for scaling ions */
        eloss = 0.0;

    ion->E -= eloss;
    ion->E = max(ion->E, 0.001 * C_EV);

    ion->p.x += d * ion->opt.sin_theta * cos(ion->fii);
    ion->p.y += d * ion->opt.sin_theta * sin(ion->fii);
    ion->p.z += d * ion->opt.cos_theta;

    ion->time += d / vel;

    if (cross_layer && ion->tlayer > target->ntarget)
        ion->status = FIN_OUT_DET;

    if (cross_layer && ion->type == SECONDARY && ion->tlayer == target->ntarget)
        ion->status = FIN_RECOIL;

    ion_finished(global, ion, target);


#ifdef DEBUG
    print_ion_position(global,ion,"M",ANYSIMULATION);
    if(cross_layer)
       print_ion_position(global,ion,"L",ANYSIMULATION);
    else if(sto_dec)
       print_ion_position(global,ion,"S",ANYSIMULATION);
    else if(cross_recdist)
       print_ion_position(global,ion,"D",ANYSIMULATION);
    else if(sc == ERD_SCATTERING)
       print_ion_position(global,ion,"R",ANYSIMULATION);
    else if(sc == MC_SCATTERING)
       print_ion_position(global,ion,"M",ANYSIMULATION);
#endif

    return (sc);
}

double inter_sto(Target_sto *stop, double vel, int mode) {
/*
 *    Interpolate the electronic stopping power or straggling value.
 *    The array is assumed to have constant spacing, so we can calculate
 *    the index.
 */

    int i = 0;
    double sto;
    double d, vlow, slow, shigh;

    d = stop->stodiv;
    assert(stop->n_sto > 0);

    if (vel >= stop->vel[stop->n_sto - 1]) {
        fprintf(stderr, "Warning: ion velocity (%e) exceeds the maximum velocity of stopping power table\n", vel);
        if (mode == STOPPING)
            return (stop->sto[stop->n_sto - 1]);
        else /* straggling */
            return (stop->stragg[stop->n_sto - 1]);
    }
    if (vel <= stop->vel[0]) {
        if (mode == STOPPING)
            return (stop->sto[0]);
        else /* straggling */
            return (stop->stragg[0]);
    }

    i = (int) (vel * d);
    assert(i >= 0 && i < stop->n_sto - 1);
    vlow = stop->vel[i];
    if (mode == STOPPING) {
        slow = stop->sto[i];
        shigh = stop->sto[i + 1];
        assert(shigh > 0); /* slow can be zero (at zero velocity), but shigh should never be. */
    } else {
        slow = stop->stragg[i];
        shigh = stop->stragg[i + 1];
    }

    sto = slow + (shigh - slow) * (vel - vlow) * d;

    return (sto);

}

int ion_finished(Global *global, Ion *ion, Target *target) {
/*
 *     Check all the conditions and decide whether to finish the ion.
 */

    if (ion->E < global->emin) {
        /*fprintf(stderr, "Stopping since %g keV < %g keV\n", ion->E/C_KEV, global->emin/C_KEV);*/
        ion->status = FIN_STOP;
    }

    if (ion->type == PRIMARY && ion->tlayer >= target->ntarget)
        ion->status = FIN_TRANS;

    if (ion->tlayer < 0) {
        if (ion->type == PRIMARY)
            ion->status = FIN_BS;
        else
            ion->status = FIN_RECOIL;
    }

    if (global->rough) {
        if (ion->type == PRIMARY && ion->p.z > (target->recmaxd + target->surface.depth))
            ion->status = FIN_MAXDEPTH;
    } else {
        if (ion->type == PRIMARY && ion->p.z > target->recmaxd)
            ion->status = FIN_MAXDEPTH;
        if (ion->type == PRIMARY && ion->scale && ion->p.z > SCALE_DEPTH)
            ion->status = FIN_MAXDEPTH;
    }

    if (ion->type == SECONDARY && ion->tlayer >= target->nlayers)
        ion->status = FIN_DET;

/*
 *  Here we should also check the case that the recoil comes out of the
 *  detector layer from sides.
 */

    return (ion->status);
}

int mc_scattering(Global *global, Ion *ion, Ion *recoil, Target *target, Detector *detector,
                  Scattering **scat, SNext *snext) {
/*
 *   Normal elastic scattering from a potential during the slowing down process 
 */

    Scattering *s;
    double theta, fii, cos_theta_cm, sin_theta_cm, sin_theta, cos_theta, n, targetA, targetZ, Ef, E_recoil, theta_recoil, sin_theta_recoil, cos_theta_recoil;
    double angle;
    int recoils = FALSE;

    global->nmc++;
    if (ion->scale) {   /* This is an ion for absolute yield scaling */
        theta = 0.0;
        cos_theta = 1.0;
        sin_theta = 0.0;
        return FALSE;
    }

    targetA = target->ele[snext->natom].A; /* TODO: pick an isotope. Problem: difficult.*/
    targetZ = target->ele[snext->natom].Z;

    s = &(scat[ion->scatindex][snext->natom]); /* !!! */

    ion->opt.e = ion->E * s->E2eps;

    if (ion->E < global->emin)
        return FALSE;

    n = ion->A / targetA;

    angle = get_angle(s, ion);

    cos_theta_cm = cos(angle);
    sin_theta_cm = sin(angle);

    theta = atan2(sin_theta_cm, cos_theta_cm + n);
    cos_theta = cos(theta);
    sin_theta = sin(theta);

    if (ion->A > targetA && theta > asin(1.0 / n)) {
        fprintf(stderr, "Scattering with M1>M2 and scattering angle exceeds asin(M2/M1). Physics fails.\n");
    }

    theta_recoil = 0.5 * (PI - angle); /* This might actually be -1.0*0.5*(PI-angle) */
    if (theta_recoil > 0.5 * PI) {
        fprintf(stderr,
                "Warning: unphysically scattered backwards. angle=%g, cos(angle)=%g, sin(angle)=%g, theta=%g, cos(theta)=%g, sin(theta)=%g, theta_recoil=%g\n",
                angle, cos_theta_cm, sin_theta_cm, theta, cos_theta, sin_theta, theta_recoil);
    }
    cos_theta_recoil = cos(theta_recoil); /* But it wouldn't affect this, since cos(-x)=cos(x) */
    sin_theta_recoil = sin(theta_recoil); /* And this is not actually used! */



    Ef = ipow2((sqrt(ipow2(targetA) - ipow2(ion->A * sin_theta)) +
                ion->A * cos_theta) / (ion->A + targetA));
    E_recoil = (1.0 - Ef) * ion->E;
    if (global->cascades && E_recoil > 2.0 *
                                       global->emin) {  /* Sufficient energy to be treated as a "true" recoil, the same way as the original particle. TODO: change the factor "2.0" to something else. It's here to fix an issue when a newly created recoil is just above minimum energy and is immediately considered to be finished. */
        if (is_in_energy_detector(global, ion, target, detector,
                                  TRUE)) { /* But we treat only things happening in our detector. */
            if (!recoil) {
                fprintf(stderr,
                        "Warning: too many recoils in a cascade (recoil of %g keV only contributes as nuclear stopping. Ion trackid %"PRIu64".\n",
                        E_recoil / C_KEV, ion->trackid);
            } else {
                if (global->output_trackpoints) {
                    output_trackpoint(global, ion, target, detector,
                                      'S'); /* This is not the best place to output a track point, but since the collision is "hard" the stopping can change a lot */
                }
                recoil->Z = targetZ;
                recoil->A = targetA; /* TODO: Pick an isotope (affects targetA too, see beginning of this function) */
                recoil->E = E_recoil;
                recoil->opt.valid = FALSE; /* TODO: check that optimization values are recomputed if used */
                recoil->w = 0.0; /* We don't want these recoils to count towards actual ERD scattering events, these are merely in the detector */
                recoil->wtmp = 0.0;
                recoil->time = ion->time;
                recoil->p = ion->p;
                recoil->lab = ion->lab;
                recoil->type = SECONDARY;
                recoil->status = NOT_FINISHED;
                recoil->nsct = 0;
                recoil->virtual = ion->virtual; /* This and the few other values following are probably not that critical? Let's just copy the information. */
                recoil->scale = ion->scale;
                recoil->tlayer = ion->tlayer;
                recoil->hist.Z = targetZ;
                recoil->hist.A = targetA;
                recoil->theta = ion->theta; /*0.0;*/ /* These will be rotated to the proper angles later in this procedure */
                recoil->fii = ion->fii; /*0.0;*/
                recoil->opt.cos_theta = ion->opt.cos_theta;
                recoil->opt.sin_theta = ion->opt.sin_theta;
                recoil->E_nucl_loss_det = 0.0;
                recoil->Ed[recoil->tlayer - target->ntarget] = recoil->E;
                recoil->dt[recoil->tlayer - target->ntarget] = recoil->time;

                /*recoil->Ed = ion->Ed;
                 recoil->dt = ion->dt;
                 recoil->hit = ion->hit;*/
                /*  save_ion_history(global, ion, recoil, detector, sc_lab, sc_ion, ion->Z*/ /* NOTE: not the recoil Z, that has no information *//*, ion->A);*/
                /* TODO: not finished.  */
                recoils = TRUE;
            }
        }
    }

#ifdef SCAT_ANAL
    if(ion->nsct < NMAXSCT)
       ion->sct_angle[ion->nsct] = theta;
    else
       fprintf(stderr,"Number of scatterings exceeds the reserved space (%i)\n",NMAXSCT);
#endif

#ifdef NO_ETRANS
    Ef = 1.0;
    recoils=FALSE;
#endif

#ifdef NO_ANGLES
    cos_theta = 1.0;
    sin_theta = 0.0;
    cos_theta_recoil = 1.0; /* Momentum is not conserved in this collision... */
    sin_theta_recoil = 0.0;
#endif

#ifdef NO_PRI_ANGLES
    if(ion->type == PRIMARY){
       cos_theta = 1.0;
       sin_theta = 0.0;
    }
#endif

#ifdef NO_SEC_ANGLES
    if(ion->type == SECONDARY){
       cos_theta = 1.0;
       sin_theta = 0.0;
       cos_theta_recoil = 1.0; /* Momentum is not conserved */
       sin_theta_recoil = 0.0;
    }
#endif


#ifndef NO_RBS_SCATLIMIT
    if (global->simtype == SIM_RBS && ion->type == SECONDARY &&
        ion->tlayer < target->ntarget) { /* Limited to sample! */
        if (ion->hist.ion_recoil.theta < theta) {
            cos_theta = 1.0;
            sin_theta = 0.0;
            Ef = 1.0;
            global->nmclarge++;
/*
         fprintf(stderr,"MC scattering larger (%.1f > %.1f) than main scattering: rejected\n",theta/C_DEG,ion->hist.ion_recoil.theta/C_DEG);
*/
            theta = 0.0;
        }
    }
#endif
    ion->E_nucl_loss_det += (1.0 - Ef) *
                            ion->E; /* All nuclear energy losses are tabulated, including the prime recoil from a cascade (above) */
    ion->E *= Ef;
    ion->opt.e = ion->E * s->E2eps;

    fii = rnd(0.0, 2.0 * PI, RND_RIGHT, RND_CONT); /* random "fii", some sort of isotropic thing going on */
    ion_rotate(ion, cos_theta, fii);
    if (recoils) {
        /* fprintf(stderr, "Scattering in detector. Ef=%g, m1=%g, m2=%g, theta_cm=%g, theta=%g, cos_theta=%g, fii=%g, theta_recoil=%g, cos_theta_recoil=%g, fii_recoil=%g\n", Ef, ion->A/C_U, targetA/C_U, angle, theta, cos_theta, fii, theta_recoil, cos_theta_recoil, fmod(fii+PI, 2.0*PI)); */
        ion_rotate(recoil, cos_theta_recoil, fmod(fii + PI, 2.0 * PI));
    }
    ion->nsct++;
    return recoils;
}

double get_angle(Scattering *scat, Ion *ion) {
/*
 *   Interpolate the cosinus of the scattering angle according to the
 *   ion energy and impact parameter.
 */

    double y, e, angle, tmp0, tmp1, tmp2;
    double elow, ehigh, ylow, yhigh, angle11, angle12, angle21, angle22;
    double angle2;

    int i = 1, j = 1;

    y = ion->opt.y;
    e = ion->opt.e;

    i = (int) ((log(e) - scat->logemin) * scat->logediv) + 1;
    j = (int) ((log(y) - scat->logymin) * scat->logydiv) + 1;

    if (i > (EPSNUM - 2)) {
        fprintf(stderr, "Warning: energy exceeds the maximum of the scattering table energy\n");
        return (0.0);
    }
    if (i < 1) {
        fprintf(stderr, "Warning: energy is below the minimum of the scattering table energy\n");
        return (0.0);
    }
    if (j < 1) {
        fprintf(stderr, "Warning: impact parameter is below the minimum of the scattering table value\n");
        printf("y: %10.3f\n",y*scat->a/C_ANGSTROM);
        return (PI); /* PI or zero? */
    }
    if (j > (YNUM - 2)) {
        fprintf(stderr,"Warning: impact parameter exceeds the maximum of the scattering table value: %i\n",j);
        printf("y: %10.3f\n",y*scat->a/C_ANGSTROM);
        return (0.0);
    }

    ylow = scat->angle[0][j];
    yhigh = scat->angle[0][j + 1];
    elow = scat->angle[i][0];
    ehigh = scat->angle[i + 1][0];

    angle11 = scat->angle[i][j];
    angle12 = scat->angle[i][j + 1];
    angle21 = scat->angle[i + 1][j];
    angle22 = scat->angle[i + 1][j + 1];

    tmp0 = 1.0 / (yhigh - ylow);

    tmp1 = angle11 + (y - ylow) * (angle12 - angle11) * tmp0;
    tmp2 = angle21 + (y - ylow) * (angle22 - angle21) * tmp0;

    angle = tmp1 + (e - elow) * (tmp2 - tmp1) / (ehigh - elow);
/*
   printf("%10.6f %10.6f\n%10.6f %10.6f\n",angle21,angle22,angle11,angle12);
   printf("%10.6f %10.6f %10.6f %10.6f\n",ehigh,ylow,ehigh,yhigh);
   printf("%10.6f %10.6f %10.6f %10.6f\n",elow,ylow,elow,yhigh);
   printf("%10.6f %10.6f %10.6f\n\n",e,y,angle);
   
   angle2 = scattering_angle(scat->pot,ion);
   printf("ANG: %10.6f %10.6f\n",angle/C_DEG,angle2/C_DEG);
*/
    return (angle);

}

void ion_rotate(Ion *p, double cos_theta, double fii) {
    double x, y, z, rx, ry, rz;
    double cosa1, cosa2, cosa3, sina1, sina2, sina3;

    double sin_theta, sin_fii, cos_fii;

    /* Change this to use the general rotate subroutine!! */

    sin_theta = sqrt(fabs(1.0 - cos_theta * cos_theta));
    cos_fii = cos(fii);
    sin_fii = sin(fii);

    x = sin_theta * cos_fii;
    y = sin_theta * sin_fii;
    z = cos_theta;

    cosa1 = p->opt.cos_theta;
    sina1 = sqrt(fabs(1.0 - cosa1 * cosa1));

    sina2 = sin(p->fii + PI / 2.0);
    cosa2 = cos(p->fii + PI / 2.0);

    cosa3 = cosa2;
    sina3 = -sina2;

    rx = x * (cosa3 * cosa2 - cosa1 * sina2 * sina3) +
         y * (-sina3 * cosa2 - cosa1 * sina2 * cosa3) +
         z * sina1 * sina2;

    ry = x * (cosa3 * sina2 + cosa1 * cosa2 * sina3) +
         y * (-sina3 * sina2 + cosa1 * cosa2 * cosa3) -
         z * sina1 * cosa2;

    rz = x * sina1 * sina3 +
         y * sina1 * cosa3 +
         z * cosa1;

    p->opt.cos_theta = rz;
    p->opt.sin_theta = sqrt(fabs(1.0 - ipow2(rz)));

    p->theta = acos(rz);
    p->fii = atan2(ry, rx);
    if (p->fii < 0.0)
        p->fii += 2.0 * PI;

}

int recdist_crossing(Global *global, Ion *ion, Target *target,
                     double dist, double *dreclayer) {
/*
 *  This subroutine checks whether we are crossing a layer in the
 *  recoil material distribution. We have separate handling for 
 *  rough and non-rough surfaces.
 */

    Point p1cross, p2cross, nextp;
    double d, d1, d2;
    int scross, s1cross, s2cross, n, i = 0, i1, i2;

    n = target->nrecdist;
    d = -1.0;  /* Default is not to have layer crossing */
    i = get_reclayer(global, target, &(ion->p));

    scross = s1cross = s2cross = FALSE;
    if (global->rough) {
        nextp.x = dist * ion->opt.sin_theta * cos(ion->fii) + ion->p.x;
        nextp.y = dist * ion->opt.sin_theta * sin(ion->fii) + ion->p.y;
        nextp.z = dist * ion->opt.cos_theta + ion->p.z;
        if (i > 0) {
            target->surface.move = target->recdist[i - 1].x;
            s1cross = surface_crossing(&(target->surface), &(ion->p), &nextp, &p1cross);
        }
        if (i < n) {
            target->surface.move = target->recdist[i].x;
            s2cross = surface_crossing(&(target->surface), &(ion->p), &nextp, &p2cross);
        }
        if (s1cross) {
            scross = s1cross;
            d1 = get_distance(&(ion->p), &p1cross);
            d = d1;
        }
        if (s2cross) {
            scross = s2cross;
            d2 = get_distance(&(ion->p), &p2cross);
            d = d2;
        }
        if (s1cross && s2cross) {
            if (d1 < d2) {
                d = d1;
            } else {
                d = d2;
            }
        }
    } else {
        nextp.z = dist * ion->opt.cos_theta + ion->p.z; /* We need only z here */
        i1 = i;
        i2 = get_reclayer(global, target, &nextp);
        if (i2 > i1)
            d = (target->recdist[i1].x - ion->p.z) / ion->opt.cos_theta;
        else if (i1 > i2)
            d = (target->recdist[i1 - 1].x - ion->p.z) / ion->opt.cos_theta;
        if (i1 != i2)
            scross = TRUE;

        if (d < 0 && scross) {
            printf("d: %10.5f\n", d / C_NM);
            printf("cos_theta: %10.5e\n", ion->opt.cos_theta);
            printf("i1,i2: %i,%i\n", i1, i2);
            printf("i1.x, i1-1.x: %10.3f %10.3f\n", target->recdist[i1].x / C_NM,
                   target->recdist[i1 - 1].x / C_NM);
            printf("ion->p.z: %10.3f\n", ion->p.z / C_NM);
        }
    }

    *dreclayer = d;

    if (scross)
        return (TRUE);
    else
        return (FALSE);

}

int recdist_nonzero(Global *global, Ion *ion, Target *target, double drec) {
/*
 *  This subroutine checks whether the recoil material distribution
 *  is nonzero at point of the next recoil event. We have separate
 *  handling for rough and non-rough surfaces.
 */

    Point nextp;
    int n, i = 0;

    n = target->nrecdist;

    if (global->rough) {
        nextp.x = drec * ion->opt.sin_theta * cos(ion->fii) + ion->p.x;
        nextp.y = drec * ion->opt.sin_theta * sin(ion->fii) + ion->p.y;
        nextp.z = drec * ion->opt.cos_theta + ion->p.z;
        i = get_reclayer(global, target, &nextp);
    } else {
        nextp.z = drec * ion->opt.cos_theta + ion->p.z; /* We need only z here */
        i = get_reclayer(global, target, &nextp);
    }

    if (i == 0 || i >= n) {
        return (FALSE);
    } else {
        if (target->recdist[i - 1].y > 0 || target->recdist[i].y > 0)
            return (TRUE);
        else
            return (FALSE);
    }
}

int get_reclayer(Global *global, Target *target, Point *p) {
    double x, y, z, xtmp, ytmp;
    int i = 0, n;

    z = p->z;
    n = target->nrecdist;

    if (global->rough) {
        xtmp = p->x + target->surface.origin.x;
        ytmp = p->y + target->surface.origin.y;
        x = xtmp * target->surface.cos_r + ytmp * target->surface.sin_r;
        y = ytmp * target->surface.cos_r - xtmp * target->surface.sin_r;
        target->surface.move = 0.0;
        z -= Z(&(target->surface), x, y);
    }

    while (i < n && z > target->recdist[i].x)
        i++;

    return (i);
}
