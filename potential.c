#include <math.h>
#include <stdlib.h>

#define NPOINTS 50
#define POTMIN 1.0e-10
#define POTERR 1.0e-3
#define MAXPOINTS 30000
#define XMAX 1000.0

typedef struct {
    double x, y;
} Point2;

typedef struct {
    int n;
    int d;
    Point2 *u;
} Potential;

double U(double);
void make_screening_table(Potential *);
double get_max_x(void);
int get_npoints(double);

void make_screening_table(Potential *pot) {
    double x, xstep, xmax;
    int i;

    xmax = get_max_x();
    pot->n = get_npoints(xmax);

    xstep = xmax / (pot->n - 1);
    pot->d = 1.0 / xstep;

    pot->u = (Point2 *) malloc(pot->n * sizeof(Point2));

    x = 0;
    for (i = 0; i < pot->n; i++) {
        pot->u[i].x = x;
        pot->u[i].y = U(x);
        x += xstep;
    }

}

double get_max_x(void) {
    double x0, x = 1.0;

    while (U(x) > POTMIN && x < XMAX)
        x *= 2.0;

    x0 = x / 2.0;
    while (U(x) < POTMIN)
        x = 0.5 * (x0 + x);

    return (x);
}

int get_npoints(double xmax) {
    double x, maxdiff, xlow, xhigh, xstep;
    double pred1, pred2, accur1, accur2, err1, err2;
    int i, npoints;

    maxdiff = 0.0;
    x = 0.0;
    xlow = 0.0;
    xhigh = xmax;
    xstep = xmax / NPOINTS;

    for (i = 0; i < NPOINTS; i++) {
        pred1 = 0.5 * (U(x) + U(x + xstep));
        accur1 = U(x + xstep / 2.0);
        err1 = fabs(1.0 - pred1 / accur1);
        if (err1 > maxdiff) {
            xlow = x;
            xhigh = x + xstep;
            maxdiff = err1;

        }
    }

    do {
        xstep = (xhigh - xlow) / 2.0;
        pred1 = 0.5 * (U(xlow) + U(xlow + xstep));
        pred2 = 0.5 * (U(xhigh) + U(xhigh - xstep));
        accur1 = U(xlow + xstep / 2.0);
        accur2 = U(xhigh - xstep / 2.0);
        err1 = fabs(1.0 - pred1 / accur1);
        err2 = fabs(1.0 - pred2 / accur2);
        if (err1 > err2)
            xhigh = xlow + xstep;
        else
            xlow = xhigh - xstep;
        xstep /= 2.0;
    } while (!(err1 < POTERR && err2 < POTERR));

    npoints = xmax / (xhigh - xlow) + 1;

    return (npoints);
}

double U(double x) {
    double value;

    value = 0.18175 * exp(-3.1998 * x) + 0.50986 * exp(-0.94229 * x) +
            0.28022 * exp(-0.4029 * x) + 0.028171 * exp(-0.20162 * x);

    return (value);
}
