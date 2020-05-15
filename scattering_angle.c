#include "general.h"
#include "prototypes.h"

#define EPS 1.0e-5
#define MAXSTEP 100
#define DEPS 1.0e-7
#define DISTEPS 1.0e-6

typedef struct {
    double x0;
    double i_y2;
    double i_ey2;
    double tmp2;
    double e;
    double y;
} Opt;

double Angint(double, Potential *, Opt *);
double U(double);
double Ut(Potential *, double);
double trapezoid(double, double, Potential *, int, Opt *);
double simpson(double, double, Potential *, Opt *);
double mindist(Potential *, Opt *);
double Psi(double, Potential *, Opt *);
double ipow2(double);

double scattering_angle(Potential *p, Ion *ion) {
    double theta;
    Opt s;

    s.e = ion->opt.e;
    s.y = ion->opt.y;
    s.i_y2 = 1.0 / ipow2(s.y);
    s.i_ey2 = 1.0 / (s.e * ipow2(s.y));

    s.x0 = mindist(p, &s);

    s.tmp2 = ipow2(s.x0) / (s.y * s.y * s.e);

    theta = C_PI - 4.0 * simpson(0.0 + DEPS, 1.0 - DEPS, p, &s);

    if (theta < 0.0)
        theta = fabs(theta);
    if (theta > C_PI)
        theta -= C_PI;

    return (theta);
}

double Angint(double u, Potential *p, Opt *s) {
    double value, u2, tmp0, tmp1, tmp2, tmp3;

    u2 = u * u;
    tmp0 = s->x0 / (1.0 - u2);
    tmp1 = Ut(p, s->x0) / s->x0 - Ut(p, tmp0) / tmp0;
/*   tmp2 = (s->x0*s->x0)/(s->y*s->y*u*u*s->e); */
    tmp2 = s->tmp2 / (u2);
    tmp3 = 2.0 - u2 + tmp2 * tmp1;

    value = 1.0 / sqrt(tmp3);

    return (value);
}

double Ut(Potential *p, double x) {
    double value, xstep;
    double xlow, ylow, yhigh;
    int i;

    if (x < 0)
        return (p->u[0].y);

    if (x > p->u[p->n - 2].x)
        return (p->u[p->n - 2].y);

    xstep = p->u[1].x - p->u[0].x;


    i = (int) (x * p->d);

    xlow = p->u[i].x;
    ylow = p->u[i].y;
    yhigh = p->u[i + 1].y;

    value = ylow + (yhigh - ylow) * (x - xlow) * p->d;

    return (value);
}

double trapezoid(double a, double b, Potential *p, int n, Opt *s) {
    static double value;
    static int nextn = 1;
    double h, x, sum = 0.0;
    int i;

    if (n == 1) {
        value = 0.5 * (b - a) * (Angint(a, p, s) + Angint(b, p, s));
        nextn = 1;
    } else {
        h = (b - a) / nextn;
        x = a + 0.5 * h;
        for (i = 1; i <= nextn; i++) {
            sum += Angint(x, p, s);
            x += h;
        }
        value = 0.5 * (value + (b - a) * sum / nextn);
        nextn *= 2;
    }

    return (value);
}

double simpson(double a, double b, Potential *p, Opt *stmp) {
    double s, st, os, ost, tmp;
    int i;

    os = ost = -1.0e30;

    for (i = 1; i <= MAXSTEP; i++) {
        st = trapezoid(a, b, p, i, stmp);
        s = (4.0 * st - ost) / 3.0;
/*
      tmp = fabs((s - os)/(PI/4.0 - s));
*/
        tmp = fabs((s - os) / os);
        if ((i > 2) && (tmp < EPS))
            break;
        os = s;
        ost = st;
    }
    if (i > MAXSTEP)
        fprintf(stderr, "Warning: Maximum number of integration points reached\n");

    return (s);
}

double mindist(Potential *p, Opt *s) {
    double x1, x2, diff = 1e6, diffold;

    x1 = (1.0 + sqrt(1.0 + 4.0 * s->y * s->y * s->e * s->e)) / (2.0 * s->e);

    while (Psi(x1, p, s) < 0.0)
        x1 *= 2;

    x2 = x1;

    do {
        x1 = x2;
        x2 = x1 - DEPS / (Psi(x1 + DEPS, p, s) / Psi(x1, p, s) - 1.0);
        diffold = diff;
        diff = fabs(x2 - x1);
    } while (diff > DISTEPS && diff < diffold);
/*   
   if(diff > DISTEPS)
   	printf("convergence stopped\n");
*/
    return (x2);
}

double Psi(double x, Potential *p, Opt *s) {
    double value;
/*
   value = (x*x)/(s->y*s->y) - (x*U(x))/(s->e*s->y*s->y) - 1.0;
*/
    value = x * x * s->i_y2 - x * Ut(p, x) * s->i_ey2 - 1.0;

    return (value);
}
