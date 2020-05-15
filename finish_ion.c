#include "general.h"
#include "prototypes.h"

void finish_ion(Global *global, Ion *ion, Target *target, Detector *detector) {
    if (ion->status == FIN_STOP)
        fprintf(global->master.fprange, "R %10.3f\n", ion->p.z / C_NM);

    if (ion->status == FIN_TRANS)
        fprintf(global->master.fprange, "T %12.6f\n", ion->E / C_MEV);
}


