#include "general.h"
#include "prototypes.h"

void finalize(Global *global) {
    int i, j;

    fprintf(global->master.fpdat, "Statistics of ion finishes: \n");

    for (i = 0; i <= SECONDARY; i++)
        for (j = 0; j < NIONSTATUS; j++)
            fprintf(global->master.fpdat, "%2i %2i: %5i\n", i, j, global->finstat[i][j]);

    if (global->simtype == SIM_RBS)
        fprintf(global->master.fpdat, "%.2f %% of MC scatterings rejected\n",
                100.0 * global->nmclarge / ((double) global->nmc));
    fclose(global->master.fperd);
    fclose(global->master.fpout);
    fclose(global->master.fpdat);
    if (global->output_trackpoints) {
        fclose(global->master.fptrack);
    }
#ifdef DEBUG
    fclose(global->master.fpdebug);
#endif
}
   
