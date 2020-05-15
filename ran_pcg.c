#include "ran_pcg.h"
#include <math.h>

double ran_pcg(long *idum) {
    static pcg32_random_t rng; /* Assuming we don't do MT */
    if (*idum < 0) { /* Old convention - seed when called with a negative number */
        uint64_t seed = (*idum < 0 ? -*idum : *idum);
        pcg32_srandom_r(&rng, seed, 1); /* Note: using hard coded 1 for initseq. */
    }
    double d = ldexp(pcg32_random_r(&rng), -32); /* d is [0,1) that has been rounded down to the nearest multiple of
 * 1/2^32 */
    return d;
}
