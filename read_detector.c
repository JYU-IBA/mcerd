#include "general.h"
#include "prototypes.h"

void check_ioerror(int c, int n, char *);

void read_detector_file(char *fname, Global *global, Detector *detector, Target *target) {
    FILE *fp;
    char dtype[10], ftype[20], foil_file[LINE], buf[LINE];
    int c, n = 0;
    foil_file[0] = '\0';
    global->virtualdet = FALSE;
    memset(detector, 0, sizeof(Detector));
    fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open the ERD-detector file %s\n", fname);
        exit(10);
    }

    c = fscanf(fp, "Detector type: %s\n", dtype);
    check_ioerror(c, 1, fname);

    if (strcmp(dtype, "TOF") == 0) {
        detector->type = DET_TOF;
    } else if (strcmp(dtype, "ENERGY") == 0) {
        detector->type = DET_ENERGY;
    } else {
        fprintf(stderr, "Only TOF and ENERGY detector supported\n");
        exit(12);
    }

    c = fscanf(fp, "Detector angle: %lf\n", &(detector->angle));
    check_ioerror(c, 1, fname);
    detector->angle *= C_DEG;

    c = fscanf(fp, "Virtual detector size: %lf %lf\n",
               &(detector->vsize[0]), &(detector->vsize[1]));
    check_ioerror(c, 2, fname);


    c = fscanf(fp, "Timing detector numbers: %i %i\n", &(detector->tdet[0]), &(detector->tdet[1]));
    if(detector->type == DET_TOF) {
        check_ioerror(c, 2, fname);
    }

    c = fscanf(fp, "Energy detector number: %i\n", &(detector->edet[0]));
    if (detector->type == DET_ENERGY || global->output_trackpoints) {
        check_ioerror(c, 1, fname);
    }

    if(fgets(buf, LINE, fp)) { /* Shitty workaround to support spaces in filenames, because this code used to rely on fscanf */
        static const char *df_line = "Description file for the detector foils: ";
        size_t df_line_len = strlen(df_line);
        if(strncmp(buf, df_line, df_line_len) == 0) {
            strncpy(foil_file, buf + df_line_len, LINE);
            foil_file[strcspn(foil_file, "\r\n")] = 0; /* Strips newlines */
        } else {
            fprintf(stderr, "Detector foil description file not given (got \"%s\", expected something starting with \"%s\")\n", buf, df_line);
            exit(2);
        }
    } else {
        exit(2);
    }

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
