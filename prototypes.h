#include <jibal.h>
#include "ran_pcg.h"

double P0(double);
double P2(double);
double P4(double);
double calc_cross(double, double, Scattering *, Potential *);
void calc_cross_sections(Global *, Scattering *, Potential *);
Point coord_transform(Point, double, double, Point, int);
void create_ion(Global *, Ion *, Target *);
void move_to_erd_detector(Global *, Ion *, Target *, Detector *);
void erd_detection(Global *, Ion *, Detector *);
void erd_init(Global *, Detector *);
int erd_scattering(Global *, Ion *, Ion *, Target *, Detector *);
void fatal_error(const char *);
void finalize(Global *);
void finish_ion(Global *, Ion *, Target *, Detector *);
void fpu(void);
void fpu_mask(void);
double gasdev(long *);
double gaussian(long);
double get_angle(Scattering *, Ion *);
double get_angle_between_points(Point *, Point *, Point *);
Rect get_beam_rect(double, Point2, Point2);
double get_cross(Ion *, Scattering *);
double get_distance(Point *, Point *);
Vector get_lab_angle(Point *, Point *);
Line get_line_params(Vector *, Point *);
int get_line_plane_cross(Line *, Plane *, Point *);
Rect get_mask_rect(double, double, Point2, double, double, double);
Plane get_plane_params(Point *, Point *, Point *);
Plane get_circ_plane(Circ *, double);
Vector get_recoiling_dir(Global *, Ion *, Ion *, Target *, Detector *);
void init_detector(Global *, Detector *);
void init_io(Global *, Ion *, Target *);
void init_params(Global *, Target *, int, char *[]);
void init_recoiling_angle(Global *, Ion *, Target *, Detector *);
double inter_sto(Target_sto *, double, int);
void ion_rotate(Ion *, double, double);
double ipow(double, int);
double ipow2(double);
int is_in_rect(Point *, Rect *);
int is_in_foil(Point *, Det_foil *);
int is_on_right_side(Point *, Vector *, Point *);
void make_screening_table(Potential *);
int mc_scattering(Global *, Ion *, Ion *, Target *, Detector *, Scattering **, SNext *);
int move_ion(Global *, Ion *, Target *, SNext *);
void next_scattering(Global *, Ion *, Target *, Scattering **, SNext *);
void output_erd(Global *, Ion *, Ion *, Target *, Detector *);
void output_data(Global *);
void print_data(Global *, Ion *, Target *, Scattering **, Detector *);
void print_ion_position(Global *, Ion *, const char *, int);
void read_detector_file(char *, Global *g, Detector *, Target *);
void read_input(Global *, Ion **, Target *, Detector *);
void read_target_file(char *, Global *, Target *);
double rnd(double, double, int, long);
void rotate(double, double, double, double, double *, double *);
double scattering_angle(Potential *, Ion *);
void scattering_table(Global *, Ion *, Target *, Scattering *, Potential *, int);
int ion_finished(Global *, Ion *, Target *);

void finish_presimulation(Global *, Detector *, Ion *);
void analyze_presimulation(Global *, Target *, Detector *);
void save_ion_history(Global *, Ion *, Ion *, Detector *, Vector, Vector, double, double);

void calc_stopping_and_straggling(Global *, Ion *, Target *, int);

void hit_virtual_detector(Global *, Ion *, Target *, Detector *, Point *, Point *);

void move_target(Target *);

int get_reclayer(Global *, Target *, Point *);

void output_trackpoint(Global *, Ion *, Target *, Detector *, char);
int is_in_energy_detector(Global *, Ion *, Target *, Detector *, int);


Ion *next_ion(Global *, Ion *, Ion *);
Ion *prev_ion(Global *, Ion *, Ion *);
void copy_ions(Ion *, Target *, int, int, int);
