#include "general.h"
#include "prototypes.h"

void check_ioerror(int c, int n, char *);

void read_detector_file(char *fname, Global *global, Detector *detector, Target *target) {
    FILE *fp;
    char dtype[10], ftype[20], foil_file[LINE], buf[LINE];
    int c, n = 0;

    global->virtualdet = FALSE;

    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "Could not open the ERD-detector file %s\n", fname);
        exit(10);
    }

    c = fscanf(fp, "Detector type: %s\n", dtype);
    check_ioerror(c, 1, fname);

    if (strcmp(dtype, "TOF") == 0) {
        detector->type = DET_TOF;
    } else if (strcmp(dtype, "GAS") == 0) {
        detector->type = DET_GAS;
    } else {
        fprintf(stderr, "Only the TOF and gas detector supported\n");
        exit(12);
    }

    c = fscanf(fp, "Detector angle: %lf\n", &(detector->angle));
    check_ioerror(c, 1, fname);
    detector->angle *= C_DEG;

    c = fscanf(fp, "Virtual detector size: %lf %lf\n",
               &(detector->vsize[0]), &(detector->vsize[1]));
    check_ioerror(c, 2, fname);

    c = fscanf(fp, "Timing detector numbers: %i %i\n",
               &(detector->tdet[0]), &(detector->tdet[1]));
    check_ioerror(c, 2, fname);

    if (global->output_trackpoints) {
        c = fscanf(fp, "Energy detector layer: %i\n",
                   &(detector->edet[0]));
        if (c != 1) {
            fprintf(stderr, "Give energy detector layer if you want to output trackpoints.\n");
        }
        check_ioerror(c, 1, fname);
    }
    c = fscanf(fp, "Description file for the detector foils: %s\n", foil_file);
    check_ioerror(c, 1, fname);

    while (fgets(buf, LINE, fp) != NULL) {
        c = fscanf(fp, "Foil type: %s\n", ftype);
        check_ioerror(c, 1, fname);
        if (strcmp(ftype, "circular") == 0)
            detector->foil[n].type = FOIL_CIRC;
        else if (strcmp(ftype, "rectangular") == 0)
            detector->foil[n].type = FOIL_RECT;
        else {
            fprintf(stderr, "Detector foil type %s not supported\n", ftype);
            exit(13);
        }
        if (detector->foil[n].type == FOIL_CIRC) {
            c = fscanf(fp, "Foil diameter: %lf\n", &(detector->foil[n].size[0]));
            check_ioerror(c, 1, fname);
            detector->foil[n].size[0] *= 0.5 * C_MM; /* We use radius! */
        } else if (detector->foil[n].type == FOIL_RECT) {
            c = fscanf(fp, "Foil size: %lf %lf\n",
                       &(detector->foil[n].size[0]), &(detector->foil[n].size[1]));
            check_ioerror(c, 2, fname);
            detector->foil[n].size[0] *= 0.5 * C_MM;
            detector->foil[n].size[1] *= 0.5 * C_MM;
        }
        c = fscanf(fp, "Foil distance: %lf\n", &(detector->foil[n].dist));
        check_ioerror(c, 1, fname);
        detector->foil[n].dist *= C_MM;
        n++;
    }
    detector->nfoils = n;
    if (detector->vsize[0] > 1.0 && detector->vsize[1] > 1.0) {
        global->virtualdet = TRUE;
    }

    read_target_file(foil_file, global, target);

    fclose(fp);

}

void check_ioerror(int c, int n, char *s) {
    if (c != n) {
        fprintf(stderr, "Problem with the input data in file %s\n", s);
        exit(2);
    }
}
