/*  
 *  Subroutine for calculating random numbers with a given distribution.
 *  
 *  make_rand_dist() prepares a data structure according to given
 *  distribution in x[] and y[] with n points.
 *
 *  rand_dist() calculates the random numbers themselves.
 */

#include "general.h"
#include "prototypes.h"


void make_rand_dist(RDist *D, double *x, double *y, int n) {
    int i, j;
    double c = 0.0, A, B, C, k;

    D->x = (double *) malloc(sizeof(double) * n);
    D->y = (double *) malloc(sizeof(double) * n);
    D->a = (double *) malloc(sizeof(double) * n);
    D->b = (double *) malloc(sizeof(double) * n);
    D->c = (double *) malloc(sizeof(double) * n);
    D->d = (int *) malloc(sizeof(int) * n);
    D->n = n;

    D->x[0] = 0.0;
    D->y[0] = x[0];

/* First we calculate integral function */

    for (i = 1; i < n; i++) {
        c += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
        D->x[i] = c;
        D->y[i] = x[i];    /* to get an inverse function we swap x and y */
    }

    if (c == 0.0) {
        fprintf(stderr, "Could not calculate cumulative function,\n");
        fprintf(stderr, "integral of the density function is zero\n");
        exit(1);
    }

/* Integral function is scaled to 1.0 */

    for (i = 0; i < n; i++) {
        D->x[i] /= c;
        y[i] /= c;
    }

/* 
 *  Here we calculate constants for interpolating the integral function.
 *  The function is of the second order between the points.
 *
 *  We also must be careful with the regions where the density function
 *  is zero, since the function must not give random number in these
 *  areas. They may have severe consequences elsewhere (eg. data without
 *  physical meaning.
 */

    for (i = 0; i < (n - 1); i++) {
        k = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        A = 0.5 * k;
        B = y[i] - k * x[i];
        C = 0.5 * k * x[i] * x[i] - y[i] * x[i];
        if (k == 0.0) {        /* Density function is constant */
            if (y[i] == 0.0) {  /* Density function is zero */
                D->d[i] = 2;
            } else {
                D->a[i] = 1.0 / B;
                D->b[i] = -(C + D->x[i]) / B;
                D->d[i] = 0;    /* constant slope */
            }
        } else {
            D->a[i] = 1.0 / A;
            D->b[i] = ipow2(B / (2.0 * A)) - (C + D->x[i]) / A;
            D->c[i] = -B / (2.0 * A);
            if (y[i + 1] >= y[i])
                D->d[i] = 1;    /* ascending slope */
            else
                D->d[i] = -1;      /* descending slope */
        }
    }

    D->x[n - 1] = 1.0;
    j = 0;
    for (i = 0; i < NINDEX; i++) {
        C = 1.0 * i / NINDEX;
        while (C >= D->x[j] && j < n)
            j++;
        j--;
        D->t[i] = j;
    }

    free(D->y);
}

double rand_dist(RDist *D) {

    int i, j;
    double value, r;

    r = (rand() + 1.0) / (RAND_MAX + 2.0);  /* random number ]0,1[ */

    i = (int) (r * NINDEX);
    j = D->t[i];

    while (r >= D->x[j])
        j++;
    j--;

    if (D->d[j] == 0)
        value = D->a[j] * r + D->b[j];
    else
        value = D->d[j] * sqrt(D->a[j] * r + D->b[j]) + D->c[j];

    return (value);

}
