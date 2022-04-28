#include "general.h"
#include "prototypes.h"
#include "rough_surface.h"

void output_erd(Global *global, Ion *cur_ion, Ion *ion, Target *target,
                Detector *detector) {
    double x, y, z, xtmp, ytmp;
    int i, j;
    int output = FALSE;
    /* Actually output if: a secondary particle (not the primary beam) is either:
     *  1. In the last layer
     *  2. In the energy detector or beyond (if we are outputting trackpoints and energy detector layer is given)
     *  3. The particle missed (exited detector, hit an aperture), but output misses is on.
     */
    if(cur_ion->type == SECONDARY && (cur_ion->tlayer == target->nlayers || is_in_energy_detector(global, cur_ion, target, detector, TRUE))) {
        output = TRUE;
    }
    if (cur_ion->status == FIN_OUT_DET && global->output_misses) {
        output = TRUE;
    }

    if(!output)
        return;
    int n = cur_ion->tlayer - target->ntarget; /* Foil number */
    int n_E = detector->edet[0] - target->ntarget; /* Energy detector */
    fprintf(global->master.fperd, "%c %c %c ",
            cur_ion->scale ? 'S' : 'R',
            cur_ion->virtual ? 'V' : 'R',
            global->simtype == SIM_ERD ? 'R' : 'S');
    if (global->output_trackpoints) {
        fprintf(global->master.fperd, "%12"PRIi64" ", cur_ion->trackid);
    }
    if (global->advanced_output) {
        fprintf(global->master.fperd, "%i ", cur_ion->status);
        fprintf(global->master.fperd, "%i ", n);
        fprintf(global->master.fperd, "%8.4f ", cur_ion->E / C_MEV);
        fprintf(global->master.fperd, "%8.4f ", cur_ion->hist.ion_E / C_MEV);
        fprintf(global->master.fperd, "%7.3f ", cur_ion->hist.ion_recoil.theta / C_DEG);
    }

    if (global->rough) {
        x = cur_ion->hist.tar_recoil.p.x;
        y = cur_ion->hist.tar_recoil.p.y;
        z = cur_ion->hist.tar_recoil.p.z;
        xtmp = x + target->surface.origin.x;
        ytmp = y + target->surface.origin.y;
        x = xtmp * target->surface.cos_r + ytmp * target->surface.sin_r;
        y = ytmp * target->surface.cos_r - xtmp * target->surface.sin_r;
        target->surface.move = 0.0;
        z -= Z(&(target->surface), x, y);
        /* TODO: output these values */
    }
#ifdef SCAT_ANAL
    for(i=0;i<ion[PRIMARY].nsct;i++){
        fprintf(global->master.fperd,"P %.4f\n",ion[PRIMARY].sct_angle[i]/C_DEG);
    }
    for(i=0;i<ion[SECONDARY].nsct;i++){
        fprintf(global->master.fperd,"S %.4f\n",ion[SECONDARY].sct_angle[i]/C_DEG);
    }
#endif
#ifdef REC_ANG
    fprintf(global->master.fperd," %10.3f",cur_ion->rec_ion_theta);
    fprintf(global->master.fperd," %10.3f",cur_ion->rec_ion_fii);
    fprintf(global->master.fperd," %10.3f",cur_ion->rec_lab_theta);
    fprintf(global->master.fperd," %10.3f",cur_ion->rec_lab_fii);
#endif

    fprintf(global->master.fperd, "%8.4f ", cur_ion->E / C_MEV);
    fprintf(global->master.fperd, "%3.0f ", cur_ion->hist.Z);
    fprintf(global->master.fperd, "%6.2f ", cur_ion->hist.A / C_U);
    fprintf(global->master.fperd, "%10.4f ", cur_ion->hist.tar_recoil.p.z / C_NM);
    fprintf(global->master.fperd, "%14.7e ", cur_ion->w);

    if (global->advanced_output) {
        fprintf(global->master.fperd, "%7.3f ", cur_ion->dt[detector->tdet[0] - target->ntarget] / C_NS);
        fprintf(global->master.fperd, "%7.3f ", cur_ion->dt[detector->tdet[1] - target->ntarget] / C_NS);
    } else {
        fprintf(global->master.fperd, "%10.3f ", (cur_ion->dt[detector->tdet[1] - target->ntarget] - cur_ion->dt[detector->tdet[0] - target->ntarget]) / C_NS); /* ToF */
    }
    if (!global->advanced_output) {
        fprintf(global->master.fperd, "%6.2f ", cur_ion->hit[detector->edet[0] - target->ntarget].x / C_MM);
        fprintf(global->master.fperd, "%6.2f ", cur_ion->hit[detector->edet[0] - target->ntarget].y / C_MM);
    } else {
        for (i = 0; i <= 4; i++) { /* This loop will print x, y coordinates of hits. */
            switch (i) {
                case 0: /* First aperture */
                    j = 0;
                    break;
                case 1:
                    j = detector->tdet[0] - target->ntarget; /* T1 */
                    break;
                case 2:
                    j = detector->tdet[1] - target->ntarget; /* T2 */
                    break;
                case 3:
                    j = detector->edet[0] - target->ntarget; /* Energy detector */
                    break;
                case 4:
                    j = n;
                    break;
                default: /* Not actually used, but jic. */
                    j = i;
                    fprintf(stderr, "Program execution in a strange place.\n");
                    break;
            }
            fprintf(global->master.fperd, "%6.2f ", cur_ion->hit[j].x / C_MM);
            fprintf(global->master.fperd, "%6.2f ", cur_ion->hit[j].y / C_MM);
        }
    }
    if (global->advanced_output) { /* Energy losses in layers */
        for (i = 0; i < (target->nlayers - target->ntarget - 1); i++) {
            fprintf(global->master.fperd, "%7.4f ", (cur_ion->Ed[i] - cur_ion->Ed[i + 1]) / C_MEV);
            if (((cur_ion->Ed[i] - cur_ion->Ed[i + 1]) / C_MEV < 0.0) &&
            i) { /* Energy gains in i=0 are possible due to virtual detector TODO: i AND virtual. */
                fprintf(stderr,
                        "Warning: ion (trackid %"PRIi64", ion_i=%i) GAINS energy between layers i=%i and i+1=%i\n",
                        cur_ion->trackid, cur_ion->ion_i, i, i + 1);
            }
        }
        if (cur_ion->tlayer == target->nlayers) {
            fprintf(global->master.fperd, "%7.4f",
                    (cur_ion->Ed[target->nlayers - target->ntarget - 1] - cur_ion->E) / C_MEV);
        } else {
            fprintf(global->master.fperd, "%7.4f", 0.0);
        }
    }
    fprintf(global->master.fperd, "\n");
    fflush(global->master.fperd);
}

void output_data(Global *global) {
    int nion, nmaxion;

    if (global->simstage == PRESIMULATION) {
        nion = global->cion;
        nmaxion = global->npresimu;
    } else {
        nion = global->cion - global->npresimu;
        nmaxion = global->nsimu - global->npresimu;
    }

    if (nmaxion > 100) {
        if (nion % (nmaxion / 100) == 0)
            fprintf(stderr, "Calculated %i of %i ions (%.0f%%)\n", nion, nmaxion, 100.0 * nion / nmaxion);
        fflush(stderr);
    }

}


void output_trackpoint(Global *global, Ion *cur_ion, Target *target, Detector *detector, char label) {
    Point pout, E_detector;
    pout = coord_transform(cur_ion->lab.p, cur_ion->lab.theta, cur_ion->lab.fii, cur_ion->p, FORW);
    E_detector = coord_transform(detector->foil[detector->edet[0] - target->ntarget].center, detector->angle, PI, pout,
                                 BACK);
    fprintf(global->master.fptrack, "%"PRIi64" %i %i %i %5.2lf %8.3lf %e %e %e %e %e %e %c\n",
            cur_ion->trackid,
            cur_ion->ion_i,
            cur_ion->tlayer - detector->edet[0],
            (int) (cur_ion->Z + 0.1),
            cur_ion->A / C_U, /* NOTE: Added on 2.4.2020 */
            cur_ion->E / C_KEV,
            E_detector.z / C_ANGSTROM,
            E_detector.x / C_ANGSTROM,
            E_detector.y / C_ANGSTROM,
            cur_ion->tlayer < target->nlayers ?
            inter_sto(&(target->layer[cur_ion->tlayer].sto[cur_ion->scatindex]),
                      sqrt(2.0 * cur_ion->E / cur_ion->A), STOPPING) * (C_ANGSTROM / C_EV) : 0.0,
            cur_ion->time / C_NS,
            cur_ion->E_nucl_loss_det / C_KEV,
            label);
}


