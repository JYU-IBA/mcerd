#include "symbols.h"
#include "mcerd_config.h"

#define MAXUNITSTRING 20
#define MAXLEN 200
#define ERR_INPUT 100

#define I_TYPE       0
#define I_ION        1
#define I_ENERGY     2
#define I_TARGET     3
#define I_DETECTOR   4
#define I_RECOIL     5
#define I_RECDIST    6
#define I_TANGLE     7
#define I_SPOTSIZE   8
#define I_MINANGLE   9
#define I_MINENERGY 10
#define I_NIONS     11
#define I_NPREIONS  12
#define I_RECAVE    13
#define I_SEED      14
#define I_RECWIDTH  15
#define I_PREDATA   16
#define I_MINSCAT   17
#define I_BDIV      18
#define I_BPROF     19
#define I_SURFACE   20
#define I_SURFSIZE  21
#define I_NSCALE    22
#define I_TRACKP    23
#define I_MISSES    24
#define I_CASCADES  25
#define I_NCASCADES 26
#define I_ADVANCED  27

#define NINPUT      28

typedef struct {
    char unit[MAXUNITSTRING];
    double value;
    int flag;
} Units;

typedef struct textline {
    char *text;
    struct textline *next;
    int read;
} TextLine;

typedef struct {
    double a[N_INPUT];
    double b[N_INPUT];
    double c[N_INPUT];
} Fvalue;

static const char *inlines[] = {
        "Type of simulation:",
        "Beam ion:",
        "Beam energy:",
        "Target description file:",
        "Detector description file:",
        "Recoiling atom:",
        "Recoiling material distribution:",
        "Target angle:",
        "Beam spot size:",
        "Minimum angle of scattering:",
        "Minimum energy of ions:",
        "Number of ions:",
        "Number of ions in the presimulation:",
        "Average number of recoils per primary ion:",
        "Seed number of the random number generator:",
        "Recoil angle width (wide or narrow):",
        "Presimulation result file:",
        "Minimum main scattering angle:",
        "Beam divergence:",
        "Beam profile:",
        "Surface topography file:",
        "Side length of the surface topography image:",
        "Number of real ions per each scaling ion:",
        "Trackpoint output:",
        "Output misses:",
        "Recoil cascades:",
        "Number of recoils in a cascade:",
        "Advanced output:"
};

#define C_DEFAULT 1.0
#define C_V0     2187691.42    /* Bohr velocity to m/s */

static Units units[] = {
        {"rad",     C_DEFAULT,                   0},
        {"radians", C_DEFAULT,                   0},
        {SYM_DEG,   C_DEG,                       0},
        {"degree",  C_DEG,                       0},
        {"deg",     C_DEG,                       0},

        {"J",       C_DEFAULT,                   0},
        {"eV",      C_EV,                        0},
        {"keV",     C_KEV,                       0},
        {"MeV",     C_MEV,                       0},

        {"kg",      C_DEFAULT,                   0},
        {"u",       C_U,                         0},

        {"s",                 1.0,               0},
        {"ms",                1.0e-3,            0},
        {SYM_US,              1.0e-6,            0},
        {"ns",                1.0e-9,            0},
        {"ps",                1.0e-12,           0},
        {"fs",                1.0e-15,           0},
        {"as",                1.0e-18,           0},

        {"m",                 1.0,               0},
        {"km",                1.0e3,             0},
        {"cm",                1.0e-2,            0},
        {"mm",                1.0e-3,            0},
        {SYM_UM,              1.0e-6,            0},
        {"nm",                1.0e-9,            0},
        {"pm",                1.0e-12,           0},
        {"fm",                1.0e-15,           0},
        {SYM_A,               1.0e-10,           0},

        {"m/s",               1.0,               0},
        {"km/s",              1000.0,            0},
        {"cm/s",              1.0e-2,            0},
        {"v0",      C_V0,                        0},
        {"c",       C_C,                         0},
        {"%c",                0.01 * C_C,        0},

        {"N",                 1.0,               0},
        {"J/m",               1.0,               0},
        {"keV/nm",  (C_KEV/C_NM) ,                    0},
        {SYM_KEV_UM,          1000.0 * (C_KEV/C_NM) , 0},
        {"keV/um",            1000.0 * (C_KEV/C_NM) , 0},
        {SYM_EV_A,            100.0 * (C_KEV/C_NM) ,  0},
        {"eV/A",              100.0 * (C_KEV/C_NM) ,  0},
        {"MeV/mm",            1000.0 * (C_KEV/C_NM) , 0},
        {SYM_KEV_UG_CM2,      1.0,               1},
        {"keV/(ug/cm2)",      1.0,               1},
        {SYM_EVCM2_1E15ATOMS, 1.0,               1},
        {SYM_MEV_MG_CM2,      1.0,               1},
        {"MeV/(mg/cm2)",      1.0,               1},
        {SYM_EV_UG_CM2,       1.0,               1},
        {"eV/(ug/cm2)",       1.0,               1},

        {SYM_G_CM3,           1000.0,            0},
        {"g/cm3",             1000.0,            0},
        {SYM_KG_M3, C_DEFAULT,                   0},
        {"kg/m3",   C_DEFAULT,                   0},

        {"END",               -1.0,              -1}
};

char *trim_space(char *);
char *get_string(char *);
int mc_isnumber(char);
char *get_unit_value(char *, double *, double);
char *get_number(char *, double *);
char *get_word(char *, int *);


void get_ion(const jibal_element *elements, char *c, double *z, double *m, Isotopes *I);
void get_atom(const jibal_element *elements, char *c, double *z);
Fvalue *read_file(char *, int, int *);


