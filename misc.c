/*
 *      Legendre polynomials for mps_erd program
 */

#include "general.h"
#include "prototypes.h"

double P0(double x) {
    return (1);
}

double P2(double x) {
    double p2;

    p2 = 0.5 * (3 * x * x - 1);

    return (p2);
}

double P4(double x) {
    double p4;

    p4 = 0.125 * (35 * x * x * x * x - 30 * x * x + 3);

    return (p4);
}

double ipow2(double x) {
    return (x * x);
}

double ipow(double x, int a) {
    int i;
    double value = 1.0;

    for (i = 0; i < a; i++)
        value *= x;

    return (value);
}

void fatal_error(const char *s) {
    fprintf(stderr, "Fatal error: %s\n", s);
    exit(10);
}

Point coord_transform(Point porig, double theta, double fii, Point pin, int flag) {
/*
 *   This routine will calculate the cartesian coordinates of point pin
 *   in another coordinate system. The origin of pin system in this
 *   other system is given by point porig. The rotation angles of pin
 *   system is given by theta and fii.
 *   Flag says which way the conversion is done.
 */


    Point pout;
    double r, in_theta, in_fii, out_theta, out_fii;

    if (flag == BACK) {
        pin.x -= porig.x;
        pin.y -= porig.y;
        pin.z -= porig.z;
    }

    r = sqrt(ipow2(pin.x) + ipow2(pin.y) + ipow2(pin.z));

    if (r == 0.0) {
        pout = porig;
        return (pout);
    }

    in_theta = acos(pin.z / r);
    in_fii = atan2(pin.y, pin.x);

    rotate(theta, fii, in_theta, in_fii, &out_theta, &out_fii);

    pout.x = r * sin(out_theta) * cos(out_fii);
    pout.y = r * sin(out_theta) * sin(out_fii);
    pout.z = r * cos(out_theta);


    if (flag == FORW) {
        pout.x += porig.x;
        pout.y += porig.y;
        pout.z += porig.z;
    }

    return (pout);
}
