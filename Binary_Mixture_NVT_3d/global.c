#include "global.h"
#include "subroutine.h"

/* ------------ For LJ Pot ----------- */
double sigma[2][2];
double epsilon[2][2];


/* ------------ For Cell list ----------- */
int     list_counter;
int     n_cellx, n_celly, n_cellz, N_cell;               
int     *cell_map[13];                               
double  max_disp2;
double  cell_cut, cell_size;


/* ------- particle variables ------- */
double eta[M_chain], p_eta[M_chain];
double Q[M_chain];
double H_eff;

/* ------- particle variables ------- */
double *x, *y, *z;
double *xu, *yu, *zu;
double *vx, *vy, *vz;
double *fx, *fy, *fz;
NEIB   *neighbour;
int    *pt_kind;



/* ------------ Imp. variable ------------ */
double sigma2, sigma_by_r2, sigma_by_r6;
double C0, C2, L, half_L;
double PE, KE, T_inst;
double pressure, virial;
int *save_time;

