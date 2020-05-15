#include <math.h>

#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif

#define max(A, B)  ((A) > (B)) ? (A) : (B)
#define min(A, B)  ((A) < (B)) ? (A) : (B)

void rotate(double, double, double, double, double *, double *);

void rotate(double theta2, double fii2, double theta1, double fii1,
            double *theta, double *fii) {
/*
 *   Definition of the angles: theta2,fii2 is the angle of the
 *   second coordinate system in the first coordinate system.
 *   Theta1,fii1 is the specific direction in the second coordinate 
 *   system, which direction in the first system is calculated
 *   in this routine (theta,fii).
 *
 *   This routine cannot be explained easily. Read eg. Goldstein
 *   about the Euler angles and coordinate transforms. Typical
 *   time for understanding this is three days ;-)
 */

    double cos_theta;
    double x, y, z, rx, ry, rz;
    double cosa1, cosa2, cosa3, sina1, sina2, sina3;

    double sin_theta, sin_fii, cos_fii;

    cos_theta = cos(theta1);
    sin_theta = sin(theta1);
    cos_fii = cos(fii1);
    sin_fii = sin(fii1);

    x = sin_theta * cos_fii;
    y = sin_theta * sin_fii;
    z = cos_theta;

    cosa1 = cos(theta2);
    sina1 = sin(theta2);

    sina2 = sin(fii2 + PI / 2.0);
    cosa2 = cos(fii2 + PI / 2.0);

    cosa3 = cosa2;
    sina3 = -sina2;

    rx = x * (cosa3 * cosa2 - cosa1 * sina2 * sina3) +
         y * (-sina3 * cosa2 - cosa1 * sina2 * cosa3) +
         z * sina1 * sina2;

    ry = x * (cosa3 * sina2 + cosa1 * cosa2 * sina3) +
         y * (-sina3 * sina2 + cosa1 * cosa2 * cosa3) -
         z * sina1 * cosa2;

    rz = x * sina1 * sina3 +
         y * sina1 * cosa3 +
         z * cosa1;

    rz = max(min(rz, 1.0), -1.0);

    *theta = acos(rz);
    if (rx != 0.0) {
        *fii = atan2(ry, rx);
    } else {
        *fii = 0.0;
    }
    if (fabs(*fii) > 2.0 * PI) {
        *fii = fmod(*fii, 2.0 * PI);
    }
    if (*fii < 0.0) {
        *fii += 2.0 * PI;
    }
}
