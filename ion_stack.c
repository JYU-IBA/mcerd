#include "ion_stack.h"

Ion *next_ion(Global *global, Ion *cur_ion, Ion *ions_moving) {
    Ion *next;
    next = cur_ion + 1;
    if (next - ions_moving == TARGET_ATOM &&
        global->simtype == SIM_RBS) { /* In RBS ions_moving[TARGET_ATOM] has a special role. Skip it. */
        next++;
    }
    if (next - ions_moving >= global->ncascades)
        return NULL; /* If there is no space left. Limits the recoil cascades. */
    return next;
}

Ion *prev_ion(Global *global, Ion *cur_ion, Ion *ions_moving) {
    Ion *prev;
    prev = cur_ion - 1;
    if (prev - ions_moving == TARGET_ATOM && global->simtype == SIM_RBS) {
        prev--;
    }
    if (cur_ion - ions_moving == 0)
        return NULL; /* We tried to go to the previous ion from the first one. Not a good idea. */
    return prev;
}


void copy_ions(Ion *ion, Target *target, int dest, int src, int copy_stopping) {
    int i, j;
    ion[dest] = ion[src];
    if (copy_stopping) {
        for (i = 0; i < target->nlayers; i++) {
            target->layer[i].sto[dest].n_sto = target->layer[i].sto[src].n_sto;
            target->layer[i].sto[dest].stodiv = target->layer[i].sto[src].stodiv;
            for (j = 0; j < target->layer[i].sto[dest].n_sto; j++) {
                target->layer[i].sto[dest].vel[j] = target->layer[i].sto[src].vel[j];
                target->layer[i].sto[dest].sto[j] = target->layer[i].sto[src].sto[j];
                target->layer[i].sto[dest].stragg[j] = target->layer[i].sto[src].stragg[j];
            }
        }
    }
}

void cascades_create_additional_ions(Global *global, Detector *detector, Target *target, Ion **ion) {
    if (!global->cascades)
        return;
    int n_primary = global->nions;
    int i, j, k;
    for (i = detector->edet[0];
         i < target->nlayers; i++) { /* For all layers, in which recoil cascades can be produced */
        for (j = 0; j < target->layer[i].natoms; j++) { /* Check each possible atom */
            int i_atom = target->layer[i].atom[j];
            double Z = target->ele[i_atom].Z;
            double A = target->ele[i_atom].A;
            fprintf(stderr, "Layer %i, atom %i is %i\n", i, j, i_atom);
            fprintf(stderr, "Atom %i is Z=%g, A=%g u\n", i_atom, Z, A / C_U);
            int found = 0;
            for (k = n_primary; k < global->nions; k++) {
                if ((*ion)[k].Z == Z && (*ion)[k].A == A) {
                    fprintf(stderr, "This ion already exists.\n");
                    found = 1;
                    break;
                }
            }
            if (!found) {
                fprintf(stderr, "Adding an ion\n");
                (*ion) = realloc(*ion, sizeof(Ion) * (global->nions + 1));
                if (*ion == NULL) {
                    fprintf(stderr, "Reallocating ions failed.\n");
                    exit(EXIT_FAILURE);
                }
                memset(&((*ion)[global->nions]), 0, sizeof(Ion));
                (*ion)[global->nions].Z = Z;
                (*ion)[global->nions].A = A;
                (*ion)[global->nions].I.n = 0; /* TODO: isotopes for recoil cascades? */
                (*ion)[global->nions].type = 0;
                global->nions++;
            }
        }
    }
}
