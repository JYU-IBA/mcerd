#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef PI
#ifdef M_PI
#define PI M_PI
#endif
#endif

#define NAMELEN 50

#define NDATA 520

#define SMALL_DIST 0.01*1e-9

#define TRUE  1
#define FALSE 0


#define R_X 0
#define R_Y 1
#define R_C 2

#define VERTICAL     1
#define NON_VERTICAL 2

#define BIG 1e100

#if 0
typedef struct {
   double x,y,z;
} Point;

typedef struct {
   double z[NDATA][NDATA];    /* Depth data leveled below zero */
   int nsize;        /* Number of data points per side */
   double size;      /* Physical length of the side */
   double depth;     /* Maximum depth in the data */
   double step;      /* size/(nsize - 1) */
   Point origin;     /* (Random) origin of the surface x-y -data */
   double cos_r,sin_r; /* cosinus and sinus of the rotation (for speed) */
   double move;      /* Distance to move the surface */
} Surface;
#endif

typedef struct {
    double x, y, z;
    double x0, y0;  /* limits in the step width in x- and y-direction */
    double dx, dy;  /* distances to the next lower step */
} SPoint;   /* Surface point */


int quick_crossing(Surface *, SPoint *, SPoint *);
int project_to_surface(SPoint *, SPoint *, double, double);
int read_afm(char *, Surface *);
double Z(Surface *, double, double);
double dp(Surface *, int, int);
int mirror(int, int);
int surface_crossing(Surface *, Point *, Point *, Point *);
void get_next_cross(SPoint *, double *, double *, double *, double, double, double, double *, int *, int);
double first_cross(SPoint *, SPoint *, double, double, double, int);
void first_vertical_cross(SPoint *, SPoint *, double, double *, double *);
int test_crossing(Surface *surface, SPoint *p1, SPoint *p2, Point *pc);
int test_crossing2(Surface *surface, SPoint *p1, SPoint *p2, Point *pc);
double ipow2(double x);
int ifloor(double);
