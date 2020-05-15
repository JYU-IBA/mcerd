#include <stdio.h>
#include <stdlib.h>

#define NINDEX 100

typedef struct {
    int t[NINDEX];
    double *x;
    double *y;
    double *a;
    double *b;
    double *c;
    int *d;
    int n;
} RDist;

void make_rand_dist(RDist *, double *, double *, int);
double rand_dist(RDist *);
