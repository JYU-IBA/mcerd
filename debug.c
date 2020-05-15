#include "general.h"
#include "prototypes.h"

void print_ion_position(Global *global, Ion *ion, const char label[], int flag) {
    /*
     *  Print ion position in the laboratory coordinates
     */

    Point plab;

    if (global->simstage == REALSIMULATION) {
        plab = coord_transform(ion->lab.p, ion->lab.theta, ion->lab.fii, ion->p,
                               FORW);
        fprintf(global->master.fpdebug, "%s %2i %2i %9.2f %11.7f %11.7f %11.7f\n",
                label, ion->type, ion->tlayer, ion->E / C_KEV, plab.x / C_MM, plab.y / C_MM,
                plab.z / C_MM);
/*      fprintf(global->master.fpdebug," %7.2f %7.2f %7.2f\n",ion->p.x/C_NM,
              ion->p.y/C_NM,ion->p.z/C_NM); */
    }
}
