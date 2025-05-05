#include "global.h"
// #include "subroutine.h"


int list_counter;

/* ------------ For Cell list ----------- */
int n_cellx, n_celly, n_cellz, N_cell;               
double cell_cut, cell_size;                          
int *cell_map[13];                               



/* ------- particle variables ------- */
double *x, *y, *z;
double *vx, *vy, *vz;
double *fx, *fy, *fz;
NEIB *neighbour;


/* ------------ Imp. variable ------------ */
double sigma2, sigma_by_r2, sigma_by_r6;
double C0, C2;
double L, PE, KE, T_inst;