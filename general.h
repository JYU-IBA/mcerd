#ifndef _MCERD_GENERAL_H_
#define _MCERD_GENERAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <ctype.h>
#include <jibal.h>

#include "rand_dist.h"
#include "symbols.h"

/* Definitions for the use of the random number generator */
#define RND_CLOSED 0
#define RND_OPEN   1
#define RND_LEFT   2
#define RND_RIGHT  3
#define RND_SEED   4    /* Set the seed number for random number generator */

#define RND_CONT   1    /* Do not reset random number generator */

#define SMALLNUM   1e-10

#define RANDOM_GENERATOR(a) ran_pcg(a) /* type of the random number generator */

/* Different ion values for ion status */

#define NOT_FINISHED  0  /* Ion continues */
#define FIN_STOP      1  /* Ion stopped */
#define FIN_TRANS     2  /* Ion transmitted through the sample */
#define FIN_BS        3  /* Ion backscattered from the sample */
#define FIN_RECOIL    4  /* Recoil came out of the sample */
#define FIN_NO_RECOIL 5  /* Could not recoil atom (because of the kinematics) */
#define FIN_MISS_DET  6  /* Recoil didn't hit one of the detector aperture */
#define FIN_OUT_DET   7  /* Recoil came out in the middle of a detector layers */
#define FIN_DET       8  /* Recoil came out from the last detector layer */
#define FIN_MAXDEPTH  9  /* Primary ion reached maximum recoil depth */
#define FIN_RECTRANS 10  /* Recoil was transmitted through the sample */
#define NIONSTATUS   11  /* Number of different ion statuses */

/*      Constants needed in program execution, which may have effects to the */
/*      physics of simulation */

#define MAXANGLE    177      /* Minimum of maximum scattering angles */
/* This is really quite small, so don't */
/* calculate backscattering with this ;-) */

/*      Constants needed with arrays etc. Shouldn't have physical meaning */

#define EPSNUM      50    /* Number of energies in scattering angle table */
#define YNUM        100   /* Number of impact factors in scattering angle table */
#define EPSIMP      200   /* Number of points in maximum impact factor array */
#define N_RAND_DIST 500   /* Number of points in random distribution */
#define N_INPUT     20000 /* Maximum number of lines in input file */
#define LINE        1024  /* Maximum length of line in input files */
#define PATH        ""    /* Path for input files */
#define NFILE       512   /* Maximum name length of the input file */
#define MAXSTO      500   /* Maximum number of values in stopping power array */
#define MAXELOSS    0.01  /* Maximum relative change in electonic energy loss */
#define MAXENEFF    2000   /* Maximum number of point in the energy eff. curve */
#define NSURFDATA   520   /* Maximum width of the surface topography image */

#define NRECDIST 500

#define MAXLAYERS 50    /* Number of different target layers */
#define MAXATOMS 20     /* Number of different atoms in the target layer */
#define MAXELEMENTS 20  /* Maximum total number of different atoms in
                           the target */
#define NANGDIST 200    /* Number of points in reaction angular distribution */

#define MAXISOTOPES 20

#define PRIMARY     0
#define SECONDARY   1
#define TARGET_ATOM 2
#define GAS_RECOIL 3

#define FORW 1
#define BACK 2

#define STOPPING   0
#define STRAGGLING 1

#define NO_SCATTERING  0
#define MC_SCATTERING  1
#define ERD_SCATTERING 2

#define SIM_ERD        1
#define SIM_RANGE      2
#define SIM_RBS        3

#define ANYSIMULATION   0
#define PRESIMULATION   1
#define REALSIMULATION  2
#define SCALESIMULATION 3  /* This is actually not used */

#define REC_NARROW 1    /* Recoiling to a narrow solid angle */
#define REC_WIDE   0    /* Recoiling to a wide solid angle */

#define SCALE_DEPTH (10.0*C_NM) /* All scaling particles are recoiled below this */

#define PRESIMU_LEVEL 0.01 /* Ratio of ions for defining recoil solid angle */

#define BEAM_NONE  0    /* No beam divergence */
#define BEAM_GAUSS 1    /* Beam profile gaussian */
#define BEAM_FLAT  2    /* Beam profile flat */
#define BEAM_DIST  3    /* Beam profile a given distribution */

/*      Different target types */

#define TARGET_FILM       0
#define TARGET_BULK       1
#define TARGET_ABSORBER   2

/*      Logical constants */

#define TRUE  1
#define FALSE 0

/*      Definitions for line and plane calculations */

#define GENERAL_PLANE 0
#define X_PLANE       1
#define Z_PLANE       2

#define L_GENERAL 0
#define L_XYPLANE 1
#define L_YAXIS   2

#define CROSS    TRUE
#define NO_CROSS FALSE

#define NSENE    2000 /* Maximum number of energies in the cross section table */
#define NSANGLE 200  /* Maximum number of angle in the cross section table */

#define TRACKPOINT_INTERVAL (10.0*C_KEV)
#define TRACKPOINT_LIMIT_SMALL (50.0*C_KEV)
#define TRACKPOINT_INTERVAL_SMALL (1.0*C_KEV)

#define MAXSTOFILEPREFIXLEN 200

#ifdef SCAT_ANAL
#define NMAXSCT 10000
#endif

/*      Macros  */

#define max(A, B)  (((A) > (B)) ? (A) : (B))
#define min(A, B)  (((A) < (B)) ? (A) : (B))

typedef struct {
    double x, y, z;      /* Structure for the 3D point */
} Point;

typedef struct {
    double x, y;        /* Structure for the 2D point */
} Point2;

/* Structure for saving the analysis of the presimulation, */
/* especially the recoiling solid angle as function of depth */

typedef struct {
    double depth;
    double angle;
    int layer;
} Presimu;

/*
 *    Structure for the general parameters needed only by the master
 *    process in the MPI-simulation
 */
typedef struct {
    char **argv;         /* Command line arguments */
    int argc;            /* Number of command line items */
    char fdata[NFILE];   /* File name of the data file */
    FILE *fperd;         /* File pointer for the ERD output data */
    FILE *fpdat;         /* File pointer for the formatted output data */
    FILE *fpout;         /* File pointer for a general output */
    FILE *fpdebug;       /* File pointer for debugging output */
    FILE *fprange;       /* File pointer for range or transmission output */
    FILE *fptrack;       /* File pointer for ion tracks in (gaseous) energy detector */
} Master;

/*    Structure for general parameters of the simulation      */

typedef struct {
    int mpi;             /* Boolean for parallel simulation */
    double E0;           /* Energy of the primary ion beam */
    int nions;           /* Number of different ions in simulation, 1 or 2 */
    int ncascades;
    int nsimu;           /* Total number of the simulated ions */
    double emin;         /* Minimum energy of the simulation */
    double ionemax;      /* Maximum possible ion energy in simulation */
    double minangle;     /* Minimum angle of the scattering */
    unsigned int seed;   /* Seed number of the random number generator */
    int cion;            /* Number of the current ion */
    int simtype;         /* Type of the simulation */
    double beamangle;    /* Angle between target surface normal and beam */
    Point2 bspot;        /* Size of the beam spot */
    int simstage;        /* Presimulation or real simulation */
    int npresimu;        /* Number of simulated ions the  presimulation */
    int nscale;          /* Number of simulated ions per scaling ions */
    int nrecave;        /* Average number of recoils per primary ion */
    int cpresimu;        /* Counter of simulated ions the the presimulation */
    Presimu *presimu;    /* Data structure for the presimulation */
    int predata;         /* Presimulation data given in a file */
    Master master;       /* Data structure for the MPI-master */
    double frecmin;
    double costhetamax;
    double costhetamin;
    int recwidth;    /* Recoiling angle width type */
    int virtualdet;      /* Do we use the virtual detector */
    char basename[NFILE];
    int finstat[SECONDARY + 1][NIONSTATUS];
    int beamdiv;        /* Angular divergence of the beam, width or FWHM */
    int beamprof;    /* Beam profile: flat, gaussian, given distribution */
    int rough;        /* Rough or non-rough sample surface */
    int nmclarge;        /* Number of rejected (large) MC scatterings (RBS) */
    int nmc;             /* Number of all MC-scatterings */
    int output_trackpoints; /* Bool TRUE/FALSE */
    int output_misses; /* Bool TRUE/FALSE */
    int cascades; /* Bool TRUE/FALSE */
    int advanced_output; /* Bool TRUE/FALSE */
    jibal jibal;
} Global;

/*      Structures for the properties of the moving ion */

typedef struct {
    int valid;           /* Boolean for the validity of these opt-variables */
    double cos_theta;    /* Cosinus of laboratory theta-angle [-1,1] */
    double sin_theta;    /* Sinus of laboratory theta-angle [-1,1] */
    double e;            /* Dimensionless energy for the next scattering */
    double y;            /* Dimensionless impact parameter  -"- */
} Ion_opt;

typedef struct {
    Point p;          /* Cartesian coordinates of the object */
    double theta, fii; /* Direction of the movement of the object */
    double v;         /* Object velocity */
} Vector;

/* Structure for saving different variables connected to the recoil history */

typedef struct {
    Vector tar_recoil; /* Recoil vector in target coord. at recoil moment */
    Vector ion_recoil; /* Recoil vector in prim. ion coord. at recoil moment */
    Vector lab_recoil; /* Recoil vector in lab. coord. at recoil moment */
    Vector tar_primary; /* Primary ion vector in target coordinates */
    Vector lab_primary; /* Primary ion vector in lab coordinates */
    double ion_E;
    double recoil_E;
    int layer;          /* Recoiling layer */
    int nsct;           /* Number of scatterings this far */
    double time;        /* Recoiling time */
    double w;           /* Statistical weight at the moment of the recoiling */
    double Z;           /* Atomic number of the secondary atom */
    double A;           /* Mass of the secondary atom in the scattering */
} Rec_hist;

typedef struct {
    double A[MAXISOTOPES]; /* List of the isotope masses in kg */
    double c[MAXISOTOPES]; /* Isotope concentrations normalized to 1.0 */
    double c_sum;          /* Sum of concentrations, should be 1.0. */
    int n;                 /* Number of different natural isotopes */
    double Am;             /* Mass of the most abundant isotope */
} Isotopes;

typedef struct {
    double Z;      /* Atomic number of ion (non-integer for compat.) */
    double A;      /* Mass of ion in the units of kg */
    double E;      /* Energy of the ion */
    Isotopes I;    /* Data structure for natural isotopes */
    Point p;       /* Three dimensional position of ion */
    double theta;  /* Laboratory theta-angle of the ion [0,PI] */
    double fii;    /* Laboratory fii-angle of the ion [0,2*PI] */
    int nsct;      /* Number of scatterings of the ion */
    int status;    /* Status value for ion (stopped, recoiled etc.) */
    Ion_opt opt;   /* Structure for the optimization-variables */
    double w;      /* Statistical weight of the ion */
    double wtmp;   /* Temporary statistical weight of the ion */
    double time;   /* Time since the creation of the ion */
    int tlayer;    /* Number of the current target layer */
    Vector lab;    /* Translation and rotation of the current coordinate */
    /* system in the laboratory coordinate system */
    int type;      /* Primary, secondary etc. */
    Rec_hist hist; /* Variables saved at the moment of the recoiling event */
    double dist;   /* Distance to the next ERD-scattering point */
    int virtual;   /* Did we only hit the virtual detector area */
    Point hit[MAXLAYERS];    /* Hit points to the detector layers */
    double Ed[MAXLAYERS]; /* Ion energy in the detector layers */
    double dt[MAXLAYERS];  /* Passing times in the detector layers */
    int scale;      /* TRUE if we have a scaling ion */
    double effrecd; /* Parameter for scaling the effective thickness of recoil material */
    int64_t trackid;
    int scatindex; /* Index of the scattering table for this ion */
    int ion_i; /* "i"th recoil in the track */
#ifdef SCAT_ANAL
    double sct_angle[NMAXSCT];
#endif
#ifdef REC_ANG
    double rec_ion_theta;
    double rec_ion_fii;
    double rec_lab_theta;
    double rec_lab_fii;
#endif
    double E_nucl_loss_det;
} Ion;

/*      Structure for impact parameter variables      */

typedef struct {
    double emin;
    double emax;
    double estep;
    double *b;
} Cross_section;

/*      Structure for potential variables      */

typedef struct {
    int n;
    int d;
    Point2 *u;
} Potential;

/*      Structures for the scattering parameters     */

typedef struct {
    double angle[EPSNUM][YNUM]; /* 2D array for scattering angles */
    Cross_section cross;        /* Data for maximum impact parameters */
    double logemin;             /* Logarithm for minimum reduced energy */
    double logymin;             /* Logarithm for minimum reduced impact parameter  */
    double logediv;             /* Logarithm for difference of reduced energies */
    double logydiv;             /* Logarithm for difference reduced impact parameters  */
    double a;                   /* Reduced unit for screening length */
    double E2eps;               /* Constant for changing energy unit to reduced energy */
/*   Potential *pot; */
} Scattering;

typedef struct {
    double d;                   /* Distance to the next scattering point */
    int natom;                  /* Number of the scattering atom */
    int r;                      /* Boolean for true mc-scattering */
    /*   (instead of eg. interface crossing ) */
} SNext;

/*      Structures for target material parameters      */

typedef struct {
    double z[NSURFDATA][NSURFDATA];    /* Depth data leveled below zero */
    int nsize;        /* Number of data points per side */
    double size;      /* Physical length of the side */
    double depth;     /* Maximum depth in the data */
    double step;      /* size/(nsize - 1) */
    Point origin;     /* (Random) origin of the surface x-y -data */
    double cos_r, sin_r; /* cosinus and sinus of the rotation (for speed) */
    double move;
} Surface;

/*      Structures for target material parameters      */

typedef struct {
    double Z;
    double A;
} Target_ele;

typedef struct {
    double vel[MAXSTO];
    double sto[MAXSTO];
    double stragg[MAXSTO];
    double stodiv;
    int n_sto;
} Target_sto;

typedef struct {
    int natoms;                /* Number of different elements in this layer */
    double dlow;               /* lower depth of the target layer */
    double dhigh;              /* higher depth of the target layer */
    int atom[MAXATOMS];        /* Array of the indices of the target elements */
    double N[MAXATOMS];        /* Array of the atomic densities */
    double Ntot;               /* Total atomic density in the layer */
    Target_sto *sto;           /* Electronic stopping for different ions */
    int type;                  /* Type of the target layer */
    int gas;              /* Whether the target layer is gas or not */
    char stofile_prefix[MAXSTOFILEPREFIXLEN];
} Target_layer;

typedef struct {
    double a, b, c;
    int type;
} Plane;

typedef struct {
    double minN;                 /* Minimum atomic density in the target */
    Target_ele ele[MAXELEMENTS]; /* Array of different elements in the target */
    Target_layer layer[MAXLAYERS]; /* Array of the target layers */
    int nlayers;                   /* Total number of target layers */
    int ntarget;         /* Number of the layers in the target itself */
    int natoms;          /* Total number of atoms in different target layers */
    Point2 recdist[NRECDIST]; /* Recoil material distribution in the target */
    int nrecdist;        /* Number of points in recoil material distribution */
    double effrecd;      /* Effective thickness of recoil material (nonzero values) */
    double recmaxd;      /* Maximum depth of the recoiling material */
    Plane plane;         /* Target surface plane */
    double efin[NRECDIST];  /* Energies below which the recoil can't reach
                              the surface */
    Point2 recpar[MAXLAYERS];  /* Recoiling solid angle fitted parameters */
    double angave;       /* Average recoiling half-angle */
    Surface surface;     /* Data structure for the rough surface */
    double cross[NSENE][NSANGLE]; /* Cross section table for scattering relative
                                   to the Rutherford cross sections as a
                        function of ion lab. energy and lab scattering angle */
    int table;  /* True if cross section table is available */
} Target;

/*      Structures for detectors variables       */

typedef struct {
    double a, b, c, d;
    int type;
} Line;

typedef struct {
    Point p[4];
} Rect;

typedef struct {
    Point p;
    double r;
} Circ;

#define DET_TOF  1
#define DET_GAS  2
#define DET_FOIL 3

#define MAXFOILS 10

#define FOIL_CIRC 0
#define FOIL_RECT 1

typedef struct {
    int type;        /* Rectangular or circular */
    int virtual;     /* Is this foil virtual */
    double dist;     /* Distance from the target center */
    double angle;    /* Angle of the foil */
    double size[2];  /* Diameter for circular, width and height for rect. */
    Plane plane;     /* Plane of the detector foil */
    Point center;    /* Point of the center of the foil */
} Det_foil;

typedef struct {
    int type;         /* TOF, GAS or FOIL */
    double angle;     /* Detector angle relative to the beam direction */
    int nfoils;       /* Number of foils in the detector */
    int virtual;      /* Do we have the virtual detector */
    double vsize[2];  /* Size of the virtual detector relative to the real */
    int tdet[2]; /* Layer numbers for the timing detectors */
    int edet[MAXLAYERS]; /* Layer number(s) for energy detector (layers) */
    double thetamax;  /* maximum half-angle of the detector (spot size incl.) */
    double vthetamax; /* maximum half-angle of the virtual detector */
    Det_foil foil[MAXFOILS];   /* Detector foil data structure */
    Det_foil vfoil;   /* Virtual detector data structure */
} Detector;

/*
 Detector angle:
TOF:
 Front foil distance:
 Front foil diameter:
 Front foil material:
 Front foil thickness:
 Back foil distance:
 Back foil diameter:
 Back foil material:
 Back foil thickness:
 Energy detector distance:
 Energy detector area:
GAS:
 Gas type:                      isobutane
 Gas pressure:                  50 mbar
 Entrance foil distance:        20.0 cm
 Entrance foil shape:           circular
 Entrance foil size:            22 mm
 Entrance foil material:        mylar
 Entrance foil thickness:       100 ug/cm2
 Ionization detector distance:  22.0 cm
 Number of electrode fields:    3
 Lengths of electrode fields    4 cm 8 cm 18 cm
 Backgammon position:           anode
 Backgammon distance:           26 cm
 Backgammon length:             8 cm
 Backgammon shape length:       1 cm
 Energy detector distance:      55 cm
 Energy detector area:          300 mm2
*/

#endif /* _MCERD_GENERAL_H_ */
