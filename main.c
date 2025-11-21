/*
 *  MCERD - a program for simulating ERD spectra using Monte Carlo methods
 *
 *   Copyright (C) 1996-2020  Kai Arstila, 2015-2025 Jaakko Julin
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include <jibal_defaults.h>
#include "general.h"
#include "prototypes.h"
#include "read_input.h"
#include "ion_stack.h"

int main(int argc, char *argv[]) {
    Global global;
    Ion *ion, *cur_ion, *previous_trackpoint_ion = NULL; /* cur_ion is a pointer to the current ion */ /* ion array contains the primary, secondary and possible unique recoils in recoil cascades */
    Ion *ions_moving; /* contains several different ions, typically the primary and secondary particles and any number of recoils  */
    Target *target = calloc(1, sizeof(Target));
    Scattering **scat; /* scat[ion][element] */
    SNext snext;
    Detector detector;
    Potential pot;
    int i, j, nscat, primary_finished;
    int64_t trackid;
    int prev_layer = 0, prev_layer_debug = 0;
    int ion_i = 0;
    int new_track = 1;
    double E_previous = 0.0, E_difference;
    fprintf(stderr, "MCERD %s compiled using JIBAL %s, current library version %s.\n", mcerd_VERSION, JIBAL_VERSION, jibal_version());
    fprintf(stderr, "MAXELEMENTS = %i, MAXISOTOPES = %i, MAXLAYERS = %i\n", MAXELEMENTS, MAXISOTOPES, MAXLAYERS);
#ifdef DEBUG
    for(i=0; i < argc; i++) {
        fprintf(stderr, "argv[%i]=%s\n", i, argv[i]);
    }
#endif
    global.jibal = jibal_init(NULL);
    if (global.jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code: %i (%s)\n", global.jibal->error,
                jibal_error_string(global.jibal->error));
        return EXIT_FAILURE;
    }
    jibal_status_print(stderr, global.jibal);
    fprintf(stderr, "Initializing parameters.\n");
    init_params(&global, target, argc, argv);

    fprintf(stderr, "Reading input files.\n");
    read_input(&global, &ion, target, &detector);

    fprintf(stderr, "Initializing output files.\n");
    init_io(&global, ion, target);

    make_screening_table(&pot);

    cascades_create_additional_ions(&global, &detector, target, &ion);
    scat = malloc(sizeof(Scattering *) * global.nions);

    fprintf(stderr, "%i ions, %i target atoms\n", global.nions, target->natoms);

#ifndef NO_GSTO_EXTRAPOLATE_ALWAYS
    global.jibal->gsto->extrapolate = TRUE; /* This needs to be before calc_stopping_and_straggling to make any difference */
#endif
    fprintf(stderr, "GSTO extrapolation turned %s.\n", global.jibal->gsto->extrapolate ? "ON" : "OFF");

    for (i = 0; i < global.nions; i++) {
        if (global.simtype == SIM_RBS && i == TARGET_ATOM)
            continue;
        ion[i].scatindex = i; /* This replaces "ion->type" in many places where ion->type was used as an index in a table. ion->type is still used to distinguish between primary and secondary in some scattering-related things */
        scat[i] = calloc(MAXELEMENTS, sizeof(Scattering));
        for (j = 0; j < target->natoms; j++) {
/*       scat[i][j].pot = &pot;  */
            fprintf(stderr, "Calculating scattering b/w ions %i (Z=%g A=%g), and target target atom %i (Z=%g A=%g)\n",
                    i, ion[i].Z, ion[i].A / C_U, j, target->ele[j].Z, target->ele[j].A / C_U);

            scattering_table(&global, &(ion[i]), target, &(scat[i][j]), &pot, j);
            calc_cross_sections(&global, &(scat[i][j]), &pot);
        }
    }


    for (j = 0; j < target->nlayers; j++) {
        target->layer[j].sto = calloc((size_t) global.nions, sizeof(Target_sto));
        for (i = 0; i < global.nions; i++) {
            if (global.simtype == SIM_RBS && i == TARGET_ATOM)
                continue;
            calc_stopping_and_straggling(&global, ion + i, target, j);
            fprintf(stderr, "Stopping calculated, scatindex=%i, layer=%i, stodiv=%g\n", ion[i].scatindex, j,
                    target->layer[j].sto[ion[i].scatindex].stodiv);
        }
    }

    jibal_gsto_print_files(global.jibal->gsto, TRUE);
    jibal_gsto_print_assignments(global.jibal->gsto);

    if (global.simtype == SIM_RBS) {
        copy_ions(ion, target, TARGET_ATOM, SECONDARY, FALSE);
        copy_ions(ion, target, SECONDARY, PRIMARY, TRUE);
    }


    init_detector(&global, &detector);

    print_data(&global, ion, target, scat, &detector);

    if (global.predata)
        init_recoiling_angle(&global, ion, target, &detector);


    trackid = ((int) ion[SECONDARY].Z) * 1000 + (global.seed % 1000);
    /*trackid *= global.nsimu;*/
    trackid *= 1000000;
    ions_moving = calloc((size_t) global.ncascades, sizeof(Ion));
    ions_moving[PRIMARY] = ion[PRIMARY]; /* These have a semi fixed position in the code */
    ions_moving[SECONDARY] = ion[SECONDARY]; /* Copy stopping too */
    if (global.simtype == SIM_RBS) {
        ions_moving[TARGET_ATOM] = ion[TARGET_ATOM];
    } /* Needed only for RBS. ERD-mode uses this slot as it pleases, RBS skips this one. */

    fprintf(stderr, "Starting simulation.\n");
    fflush(stdout);
    fflush(stderr);
    for (global.cion = 0; global.cion < global.nsimu; global.cion++) {

        output_data(&global);

        cur_ion = ions_moving + PRIMARY;
#ifdef DEBUG
        fprintf(global.master.fpdebug, "\nT Moving to primary %p, i=%i\n", cur_ion, (int)(cur_ion-ions_moving));
#endif

        create_ion(&global, cur_ion, target);
        if (global.rough) {
            move_target(target);
        }

        primary_finished = FALSE;
        while (!primary_finished) {

            next_scattering(&global, cur_ion, target, scat, &snext);
            nscat = move_ion(&global, cur_ion, target, &snext);

            if (nscat == ERD_SCATTERING) {
                if (erd_scattering(&global, ions_moving + PRIMARY, ions_moving + SECONDARY, target, &detector)) {
                    cur_ion = ions_moving + SECONDARY;
#ifdef DEBUG
                    fprintf(global.master.fpdebug, "\nT Moving (ERD) to %p, i=%i, Z=%g\n", cur_ion, (int)(cur_ion-ions_moving), cur_ion->Z);
#endif
                }
            }

            if (cur_ion->status == FIN_RECOIL || cur_ion->status == FIN_OUT_DET) {
                if (global.simstage == PRESIMULATION) {
                    finish_presimulation(&global, &detector, cur_ion);
                    cur_ion = ions_moving + PRIMARY;
#ifdef DEBUG
                    fprintf(global.master.fpdebug, "\nT Moving (OUT) to %p, i=%i, Z=%g\n", cur_ion, (int)(cur_ion-ions_moving), cur_ion->Z);
#endif
                } else {
                    move_to_erd_detector(&global, cur_ion, target, &detector);
                }
            }

            if (global.simstage == PRESIMULATION && global.cion == (global.npresimu - 1)) {
                analyze_presimulation(&global, target, &detector);
                init_recoiling_angle(&global, ion, target, &detector);
            }

            if (nscat == MC_SCATTERING && cur_ion->status == NOT_FINISHED && !global.nomc) {
                if (mc_scattering(&global, cur_ion, next_ion(&global, cur_ion, ions_moving), target, &detector, scat,
                                  &snext)) {
                    cur_ion = next_ion(&global, cur_ion, ions_moving);
                    int found = FALSE;
                    for (i = 0; i < global.nions; i++) {
                        if (i == TARGET_ATOM &&
                            global.simtype == SIM_RBS) /* Target atom doesn't move and isn't valid here. */
                            continue;
                        if (round(ion[i].Z) == round(cur_ion->Z) && round(ion[i].A / C_U) == round(cur_ion->A /
                                                                                                   C_U)) {  /* Note that target (recoiling gas atom) is an average mass. Rounding gives probably 12C and 1H. FIXME this is a horrible, horrible bit of code right here. Hurts even to think about it. */
                            /*fprintf(stderr, "Recoil cascade from Z=%g, identified scattering index %i\n", cur_ion->Z, i);*/
                            found = TRUE;
                            cur_ion->scatindex = i;
                        }
                    }
                    if (!found) {
                        fprintf(stderr,
                                "Warning: recoil cascade not possible, since recoiling ion Z=%g and A=%g u are not in ion table (and therefore not in scattering table or stopping/straggling tables). I am ignoring this recoil.\n",
                                cur_ion->Z, cur_ion->A / C_U);
                        cur_ion = prev_ion(&global, cur_ion, ions_moving);
                    } else {
                        ion_i++;
                        cur_ion->ion_i = ion_i;
                        cur_ion->trackid = trackid;
                    }
#ifdef DEBUG
                    fprintf(global.master.fpdebug, "\nT Moving (SCT) to %p, i=%i, ion_i=%i\n", cur_ion, (int)(cur_ion-ions_moving), ion_i);
#endif
                }
            }
#ifdef DEBUG
            if(cur_ion->tlayer != prev_layer_debug) {
              fprintf(global.master.fpdebug, "T Entering layer %i\n", cur_ion->tlayer);
              prev_layer_debug=cur_ion->tlayer;
            }
#endif

            if (global.output_trackpoints) {
                if (is_in_energy_detector(&global, cur_ion, target, &detector, TRUE) &&
                    !cur_ion->scale) { /* Non-scale ion is in energy detector */
                    if (new_track) {
                        fprintf(global.master.fptrack, "\n\n");
                        E_previous = 0.0;
                        trackid++;
                        ion_i = 0;
                        new_track = FALSE;
                        previous_trackpoint_ion = NULL;
                        E_previous = 0.0;
                        cur_ion->trackid = trackid;
                        cur_ion->ion_i = ion_i;
                        cur_ion->E_nucl_loss_det = 0.0;
#ifdef DEBUG
                        fprintf(global.master.fpdebug, "T Ion entered energy detector. Trackid %"PRIi64".\n", trackid);
                        fflush(global.master.fpdebug);
#endif
                    }
                    E_difference = fabs(cur_ion->E -
                                        E_previous); /* Absolute value of the energy difference since a track point was last written */
                    if (cur_ion != previous_trackpoint_ion &&
                        previous_trackpoint_ion) { /* If the "cur_ion" has changed, output a point for the previous one. Recoil cascades stay neater this way. The position and time of these should be the same and energy of the "previous_trackpoint_ion" is post-collision. */
                        output_trackpoint(&global, previous_trackpoint_ion, target, &detector, 'P');
                    }
                    if (cur_ion != previous_trackpoint_ion || cur_ion->tlayer != prev_layer ||
                        E_difference > TRACKPOINT_INTERVAL || (cur_ion->E<TRACKPOINT_LIMIT_SMALL &&
                                                                          E_difference>TRACKPOINT_INTERVAL_SMALL)) { /* At low energies, make points with smaller interval. */

                        output_trackpoint(&global, cur_ion, target, &detector,
                                          (cur_ion != previous_trackpoint_ion) ? 'N' : (cur_ion->tlayer != prev_layer
                                                                                        ? 'L' : 'E'));
                        E_previous = cur_ion->E;
                        previous_trackpoint_ion = cur_ion;
                        prev_layer = cur_ion->tlayer;
                    }
                }
            }
            if (cur_ion->type == SECONDARY && cur_ion->status != NOT_FINISHED) {
                global.finstat[SECONDARY][cur_ion->status]++;
            }

            while (ion_finished(&global, cur_ion, target)) { /* TODO: recoil cascades and finishing? */
#ifdef DEBUG
                fprintf(global.master.fpdebug, "T We are finished with this ion! Ion status is %i.\n", cur_ion->status);
                fflush(global.master.fpdebug);
#endif
                if (global.output_trackpoints) {
                    if (is_in_energy_detector(&global, cur_ion, target, &detector, TRUE) && !cur_ion->scale) {
                        output_trackpoint(&global, cur_ion, target, &detector, 'F');
                        fflush(global.master.fptrack);
                        previous_trackpoint_ion = cur_ion; /*  Upon finishing output a trackpoint. It should be a valid point since the ion hasn't (?) moved. Label 'C' for continue. */
                    }
                }
                if (new_track == FALSE) {
                    cur_ion->trackid = trackid;
                } else {
                    cur_ion->trackid = 0; /* No track was made for this one (either it didn't make it to the energy detector or was a scaling ion) */
                }

                if (cur_ion - ions_moving <= SECONDARY) { /* Either primary or the first secondary (not a recoil cascade) */
                    output_erd(&global, cur_ion, ion, target, &detector);
                    new_track = TRUE;
                }
                if (cur_ion->type == PRIMARY) {
                    primary_finished = TRUE;
                    break;
                }
                cur_ion = prev_ion(&global, cur_ion, ions_moving); /* Moving to previous ion */
                if (cur_ion->type != PRIMARY && global.output_trackpoints) {
                    output_trackpoint(&global, cur_ion, target, &detector,
                                      'C'); /* Outputting a trackpoint if this (previous) is not primary. In other words this happens only when we go from a recoil cascade back to previous ion in the cascade, all the way to the "prime" secondary. */
                    previous_trackpoint_ion = cur_ion;
                }
#ifdef DEBUG
                fprintf(global.master.fpdebug, "\nT Moving to (FIN) %p, i=%i\n", cur_ion, (int)(cur_ion-ions_moving));
#endif
            }
        }
#ifdef DEBUG
        fprintf(global.master.fpdebug, "T We are finished with this primary ion! Ion status is %i.\n", cur_ion->status);
#endif
        global.finstat[PRIMARY][cur_ion->status]++;
        finish_ion(&global, cur_ion, target, &detector);
    }

    finalize(&global);
    free(target);
    exit(0);

}
