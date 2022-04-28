#ifdef Linux
#include <fpu_control.h>
#endif

#include "general.h"
#include "rough_surface.h"

#ifndef C_NM
#define C_NM 1.0e-9
#endif

int surface_crossing(Surface *surface, Point *P1, Point *P2, Point *pc) {
    SPoint *p1, *p2, *p3;
    double ky, by, kz, bz, Dx, Dy, Dz, Dc, nx, ny, nc, sign, dend, s, cos_rot, sin_rot;
    int end, type, cross = FALSE, qtest, off_surf;

    p1 = (SPoint *) malloc(sizeof(SPoint));
    p2 = (SPoint *) malloc(sizeof(SPoint));
    p3 = (SPoint *) malloc(sizeof(SPoint));

    /*
     *  Make the coordinate transform from the target coordinates to the
     *  rough surface coordinates (to randomize the surface handling)
     */

    p3->x = P1->x + surface->origin.x;
    p3->y = P1->y + surface->origin.y;

    cos_rot = surface->cos_r;
    sin_rot = surface->sin_r;

    p1->x = p3->x * cos_rot + p3->y * sin_rot;
    p1->y = p3->y * cos_rot - p3->x * sin_rot;
    p1->z = P1->z - surface->move;

    p3->x = P2->x + surface->origin.x;
    p3->y = P2->y + surface->origin.y;

    p2->x = p3->x * cos_rot + p3->y * sin_rot;
    p2->y = p3->y * cos_rot - p3->x * sin_rot;
    p2->z = P2->z - surface->move;

    /*
     * We project the point to the depth field of the surface in order
     * to speed up the cases where points are very far away
     */

    off_surf = project_to_surface(p1, p2, 0.0, surface->depth);

    if (off_surf) {
        free(p1);
        free(p2);
        free(p3);

        return (FALSE);
    }

    s = surface->step;

    ky = by = kz = bz = 0.0;

    p1->x0 = floor(p1->x / s) * s;
    p1->y0 = floor(p1->y / s) * s;
    p2->x0 = floor(p2->x / s) * s;
    p2->y0 = floor(p2->y / s) * s;

    p1->dx = p1->x - p1->x0;
    p1->dy = p1->y - p1->y0;
    p2->dx = p2->x - p2->x0;
    p2->dy = p2->y - p2->y0;

    if (fabs(p1->x0 - p2->x0) < 0.1 * s && fabs(p1->y0 - p2->y0) < 0.1 * s) { /* same square */
        qtest = quick_crossing(surface, p1, p2); /* Quick test for trivial cases */
        if (qtest) {
            free(p1);
            free(p2);
            free(p3);
            return (FALSE);
        }
        if ((p1->dy - p1->dx) * (p2->dy - p2->dx) > 0.0) { /* same triangle */
            cross = test_crossing(surface, p1, p2, pc);
        } else { /* same square, different triangle */
            if (p2->x == p1->x) {
                p3->x = p1->x;
                p3->y = p1->x0 + p1->y0;
                p3->z = p1->z + (p3->x - p1->x) * (p2->z - p1->z) / (p2->x - p1->x);
            } else {
                ky = (p2->y - p1->y) / (p2->x - p1->x);
                by = p1->dy - ky * p1->dx;
                p3->x = p1->x0 + by / (1.0 - ky);
                p3->y = p1->y0 + p3->x - p1->x0;
                p3->z = p1->z + (p3->x - p1->x) * (p2->z - p1->z) / (p2->x - p1->x);
            }
            cross = test_crossing(surface, p1, p3, pc);
            if (!cross)
                cross = test_crossing(surface, p3, p2, pc);
        }
    } else {                /* Different square */
        if (p2->x == p1->x) {  /* line between the points is vertical */
            if (p2->y > p1->y)
                sign = 1.0;
            else
                sign = -1.0;
            Dy = sign * s;
            Dc = sign * s;
            Dx = BIG;
            nx = BIG;
            first_vertical_cross(p1, p2, s, &ny, &nc);
            dend = p2->y - p1->y;
            Dz = Dy * (p2->z - p1->z) / (p2->y - p1->y);
            type = VERTICAL;
        } else {
            if (p2->x > p1->x)
                sign = 1.0;
            else
                sign = -1.0;
            ky = (p2->y - p1->y) / (p2->x - p1->x);
            by = p1->y - ky * p1->x;

            kz = (p2->z - p1->z) / (p2->x - p1->x);
            bz = p1->z - kz * p1->x;

            Dx = sign * s;
            if (ky == 0)
                Dy = BIG;
            else
                Dy = sign * s / fabs(ky);
            if (ky == 1.0)
                Dc = BIG;
            else
                Dc = sign * s / (1.0 - ky);

            nx = first_cross(p1, p2, s, ky, by, R_X);
            ny = first_cross(p1, p2, s, ky, by, R_Y);
            nc = first_cross(p1, p2, s, ky, by, R_C);

            Dz = Dx * (p2->z - p1->z) / (p2->x - p1->x);

            dend = p2->x - p1->x;

            type = NON_VERTICAL;
        }

        *p3 = *p1;

        end = FALSE;

        while (!end && !cross) {
            get_next_cross(p3, &nx, &ny, &nc, Dx, Dy, Dc, &dend, &end, type);
            if (type == NON_VERTICAL) {
                p3->y = ky * p3->x + by;
            }
            p3->z = kz * p3->x + bz;
            cross = test_crossing(surface, p1, p3, pc);
            *p1 = *p3;
        }
    }

    /*
     * We have to translate back the crossing point to the original
     *  coordinates
     */

    if (cross) {
        cos_rot = surface->cos_r;
        sin_rot = -surface->sin_r;

        p3->x = pc->x * cos_rot + pc->y * sin_rot;
        p3->y = pc->y * cos_rot - pc->x * sin_rot;
        pc->x = p3->x - surface->origin.x;
        pc->y = p3->y - surface->origin.y;
        pc->z += surface->move;
    }

    free(p1);
    free(p2);
    free(p3);

    return (cross);
}

int project_to_surface(SPoint *p1, SPoint *p2, double low, double high) {
    SPoint *pp1, *pp2;
    double dx, dy, dz, ddz;

    if (p1->z < p2->z) {
        pp1 = p1;
        pp2 = p2;
    } else {
        pp2 = p1;
        pp1 = p2;
    }

    if (pp1->z >= low && pp2->z <= high)
        return (FALSE); /* Both points in the field of the surface heights */

    if (pp2->z < low || pp1->z > high)
        return (TRUE); /* Both points on the same side of the surface */

    dx = pp2->x - pp1->x;
    dy = pp2->y - pp1->y;
    dz = pp2->z - pp1->z;

    if (pp2->z > high) {
        ddz = (high - pp1->z) / dz;
        pp2->x = pp1->x + ddz * dx;
        pp2->y = pp1->y + ddz * dy;
        pp2->z = pp1->z + ddz * dz;
    }
    if (pp1->z < low) {
        ddz = (low - pp1->z) / dz;
        pp1->x = pp1->x + ddz * dx;
        pp1->y = pp1->y + ddz * dy;
        pp1->z = pp1->z + ddz * dz;
    }

    return (FALSE);
}

int quick_crossing(Surface *surface, SPoint *p1, SPoint *p2) {
    /*
     *   We check here if the points are off from the surface square
     */

    return (FALSE);
}

int test_crossing(Surface *surface, SPoint *p1, SPoint *p2, Point *pc) {
    double p, z1, z2, kp, ks, dz, r, d;
    int c = FALSE;

    z1 = Z(surface, p1->x, p1->y);
    z2 = Z(surface, p2->x, p2->y);
    dz = p2->z - p1->z;

    p = (p1->z - z1) * (p2->z - z2);

    if ((p1->x == p2->x) && (p1->y == p2->y)) {
        if (fabs(z1 - p1->z) <= fabs(dz) && fabs(z2 - p1->z) <= fabs(dz)) {
            /* surface is between the points */
            pc->x = p1->x;
            pc->y = p1->y;
            pc->z = z1;
            return (TRUE);
        } else
            return (FALSE);
    }

    if (p < 0.0) {
        c = TRUE;
        r = 1.0 / sqrt(ipow2(p2->x - p1->x) + ipow2(p2->y - p1->y));
        kp = r * dz;
        ks = r * (z2 - z1);
        if (kp == ks)         /* points and surface are parallel */
            return (FALSE);
        d = (z1 - p1->z) / (kp - ks);
        pc->z = d * kp + p1->z;       /* coordinates of the crossing point */
        pc->x = d * r * (p2->x - p1->x) + p1->x;
        pc->y = d * r * (p2->y - p1->y) + p1->y;
    } else {
        return (FALSE);
    }

    return (c);
}

void get_next_cross(SPoint *pc, double *nx, double *ny, double *nc,
                    double Dx, double Dy, double Dc,
                    double *dend, int *end, int type) {
    double d;

    if (fabs(*nx) < fabs(*ny)) {
        if (fabs(*nc) < fabs(*nx)) { /* nc smallest */
            d = *nc;
            *nc += Dc;
        } else {                   /* nx smallest */
            d = *nx;
            *nx += Dx;
        }
    } else {
        if (fabs(*nc) < fabs(*ny)) { /* nc smallest */
            d = *nc;
            *nc += Dc;
        } else {                   /* ny smallest */
            d = *ny;
            *ny += Dy;
        }
    }

    *nx -= d;
    *ny -= d;
    *nc -= d;

    if (fabs(*dend) < fabs(d)) {   /* end point is the closest */
        *end = TRUE;
        if (type == VERTICAL)
            pc->y += *dend;
        else
            pc->x += *dend;
    } else {
        *end = FALSE;
        if (type == VERTICAL)
            pc->y += d;
        else
            pc->x += d;
    }

    *dend -= d;

}

double first_cross(SPoint *p1, SPoint *p2, double s, double k, double b, int type) {
    double value = 0, y, c;

    switch (type) {
        case R_X:
            if (p1->x == p2->x)
                value = BIG;
            else if (p1->x < p2->x) /* going to the right */
                value = p1->x0 + s;
            else
                value = p1->x0;
            break;
        case R_Y:
            if (p1->y == p2->y)
                value = BIG;
            else {
                if (p1->y < p2->y) /* going up */
                    y = p1->y0 + s;
                else              /* going down */
                    y = p1->y0;
                value = (y - b) / k;
            }
            break;
        case R_C:
            if (k == 1)  /* parallel with the diagonals */
                value = BIG;
            else {
                if (p1->dx > p1->dy) {  /* lower triangle */
                    if ((p1->x - p2->x + p2->y) < p1->y) {  /* going to the right */
                        c = p1->y0 - (p1->x0 + s);
                    } else { /* going to the left */
                        c = p1->y0 - p1->x0;
                    }
                } else { /* upper triangle */
                    if ((p1->x - p2->x + p2->y) < p1->y) {  /* going to the right */
                        c = p1->y0 - p1->x0;
                    } else { /* going to the left */
                        c = p1->y0 + s - p1->x0;
                    }
                }
                value = (b - c) / (1 - k);
            }
            break;
    }

    value -= p1->x;

    return (value);
}

void first_vertical_cross(SPoint *p1, SPoint *p2, double s,
                          double *ny, double *nc) {
    if (p2->y > p1->y) /* upwards */
        *ny = p1->y0 + s;
    else
        *ny = p1->y0;

    if (p1->dx > p1->dy) {  /* lower triangle */
        if (p2->y > p1->y) /* upwards */
            *nc = p1->x + p1->y0 - p1->x0;
        else
            *nc = p1->x + p1->y0 - s - p1->x0;
    } else {
        if (p2->y > p1->y) /* upwards */
            *nc = p1->x + p1->y0 + s - p1->x0;
        else
            *nc = p1->x + p1->y0 - p1->x0;
    }

    *ny -= p1->y;
    *nc -= p1->y;
}

double

Z(Surface *surface, double x, double y) {
    Point p1, p2;
    double z1, z2, z, x0, y0;
    int i, j;

    i = ifloor(x / surface->step);
    j = ifloor(y / surface->step);

    x0 = i * surface->step;
    y0 = j * surface->step;

    /* p1 is a point on the diagonal corresponding x,y */

    p1.x = x0 + (y - y0);
    p1.y = y;
    z1 = dp(surface, i, j);
    z2 = dp(surface, i + 1, j + 1);
    p1.z = z1 + (p1.x - x0) * (z2 - z1) / surface->step;

    if (x > p1.x) {     /* x in the lower right triangle */
        p2.x = x0 + surface->step;
        p2.y = y;
        z1 = dp(surface, i + 1, j);
        z2 = dp(surface, i + 1, j + 1);
        p2.z = z1 + (y - y0) * (z2 - z1) / surface->step;
    } else {          /* x in the upper left triangle */
        p2 = p1;
        p1.x = x0;
        p1.y = y;
        z1 = dp(surface, i, j);
        z2 = dp(surface, i, j + 1);
        p1.z = z1 + (y - y0) * (z2 - z1) / surface->step;
    }

    if (p2.x == p1.x)  /* x,y in the corner */
        z = p1.z;
    else
        z = p1.z + (x - p1.x) * (p2.z - p1.z) / (p2.x - p1.x);

    return (z);

}

double dp(Surface *surface, int i, int j) {
    int x, y;

    x = mirror(i, surface->nsize);
    y = mirror(j, surface->nsize);

    return (surface->z[x][y]);
}

int mirror(int i, int nsize) {
    int x, f, p = 0;
    nsize--;
    if (i > 0 && i < nsize) {
        x = i;
    } else {
        if (i < 0)
            p = 1;
        f = ifloor(((double) i + p - 0.5) / nsize);
        if (f % 2 == 0) {
            x = i - f * nsize;
        } else {
            x = (f + 1) * nsize - i;
        }
    }
    return (x);
}

int read_afm(char *fname, Surface *surface) {
    FILE *fp;
    double *data, min = 1e20, max = -1e20;
    int i, j, n;

    data = (double *) malloc(NDATA * NDATA * sizeof(double));


    fprintf(stderr, "Opening AFM-file %s\n", fname);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", fname);
        exit(12);
    }

    i = 0;
    while (i < NDATA * NDATA && fscanf(fp, "%lf", data + i) == 1)
        i++;

    if (i == NDATA * NDATA) {
        fprintf(stderr, "Not all data points could be read, exiting..\n");
        exit(11);
    }

    n = ((int) sqrt((double) i) + 0.5);

    if (n * n != i) {
        fprintf(stderr, "Number of data points is not square of an integer\n");
        exit(13);
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            surface->z[i][j] = -data[i + j * n] * C_NM; /* we have to invert the data */
            if (surface->z[i][j] < min)
                min = surface->z[i][j];
            if (surface->z[i][j] > max)
                max = surface->z[i][j];
        }
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            surface->z[i][j] += (0.001 * C_NM - min);
        }
    }

    surface->nsize = n;

    surface->depth = (max - min + 0.01 * C_NM);

    free(data);

    return (1);

}

#if 0
double ipow2(double x)
{
   return(x*x);
}
#endif

int ifloor(double x) {
    double f;
    int i;
    f = floor(x);
    if (f < -0.5) {
        i = (int) (f - 0.5);
    } else {
        i = (int) (f + 0.5);
    }
    return (i);
}
