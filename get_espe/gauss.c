#include <math.h>
#include "gauss.h"

double gaussian() {
    double value;
    static pcg32_random_t rng = {0, 0};
    if (rng.state == 0) {
        pcg32_srandom_r(&rng, 12345, 1); /* 12345 guaranteed to be a random number */
    }
    value = gasdev_pcg(&rng);
    return (value);
}

double gasdev_pcg(pcg32_random_t *rng) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (iset == 0) {
        do {
            v1 = 2.0 * ldexp(pcg32_random_r(rng), -32) - 1.0;
            v2 = 2.0 * ldexp(pcg32_random_r(rng), -32) - 1.0;
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
