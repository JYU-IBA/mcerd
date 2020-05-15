#include "general.h"
#include "prototypes.h"

double rnd(double low, double high, int period, long seed) {
    double value = 0.0, length;

    length = high - low;

    if (length < 0.0)
        fatal_error("Random length negative or zero\n");

    switch (period) {
        case RND_CLOSED:
            do {
                value = length * RANDOM_GENERATOR(&seed) + low;
            } while (value < low || value > high);
            break;
        case RND_OPEN:
            do {
                value = length * RANDOM_GENERATOR(&seed) + low;
            } while (value <= low || value >= high);
            break;
        case RND_LEFT: /* left closed, right open */
            do {
                value = length * RANDOM_GENERATOR(&seed) + low;
            } while (value < low || value >= high);
            break;
        case RND_RIGHT: /* right closed, left open */
            do {
                value = length * RANDOM_GENERATOR(&seed) + low;
            } while (value <= low || value > high);
            break;
        case RND_SEED:
            RANDOM_GENERATOR(&seed);
            break;
        default:
            value = low - high;
            break;
    }


    if (value >= low && value <= high)
        return (value);
    else {
        fprintf(stderr, "%f %f %f %i %li\n", low, high, value, period, seed);
        fatal_error("Illegal value from random number generator\n");
        return (-10.0);      /* Actually we should never get here */
    }

}

double gaussian(long seed) {
    double value;

    value = gasdev(&seed);

    return (value);
}

double gasdev(long *idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (iset == 0) {
        do {
            v1 = 2.0 * RANDOM_GENERATOR(idum) - 1.0;
            v2 = 2.0 * RANDOM_GENERATOR(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}

