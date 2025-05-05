#ifndef GLOBAL_H
    #define GLOBAL_H
    

    /* -------- Libaries -------- */
    #include<stdio.h>
    #include<math.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    /* -------------------------- */


    /* -------- System Variable -------- */
    #define     N              128                       /* Total Number of Particles */
    #define     max_step       (int) 1e6                 /* Step upto Simulation runs */
    #define     rho            1.20                      /* Density of the system --- */
    #define     vol            (double) (N / rho)        /* volume of the system ---- */
    #define     dt             0.001                     /* Time step for integration */
    #define     half_dt        (0.50) * dt               /* half time step for int -- */
    #define     init_temp      0.50                      /* Initial Temp. of system - */
    /* --------------------------------- */


    /* -------- For LJ Potential -------- */
    #define     epsilon        1.0                       /* Interation strength of LJ Pot */                 
    #define     sigma          1.0                       /* effective radius of particles */
    #define     r_c            2.5                       /* Cut-Off for LJ Potential ---- */
    #define     sigma_by_rc    (double) (sigma/r_c)      /* sigma to r_cutoff ratio ----- */
    /* ---------------------------------- */


    /* -------- For Verlet list -------- */
    #define     max_neigh      50                        /* strictly more than int(PI * rl^2 * rho) */
    #define     skin_depth     0.3 * sigma               /* Skin Width for Verlet List ------------ */
    #define     rl             (r_c + skin_depth)        /* Radius for NBR list (r_c + skin) -----  */

    typedef struct {                                     /* NEIB data type has two branch */
        int list[max_neigh];                             /* one array branch  */
        int max_counter;                                 /* one scalar branch */
    }NEIB;

    extern int list_counter;                             /* counter to control after how many steps list will update */ 
    /* --------------------------------- */



    /* ------------ For Cell list ----------- */
    extern  int         n_cellx;                         /* Numbers of cell in x direction */
    extern  int         n_celly;                         /* Numbers of cell in y direction */               
    extern  int         N_cell;                          /* Total numbers of cell of sys.  */
    extern  double      cell_cut;                        /* cell size we want ------------ */
    extern  double      cell_size;                       /* true size of cell ------------ */
    extern  int         *cell_map[4];                    /* 4 neigh cells for 2d --------- */
    /* -------------------------------------- */


    /* ------- particle variables ------- */
    extern  double  *x, *y;                              /* Position arrays ----- */
    extern  double  *vx, *vy;                            /* Velocity arrays ----- */
    extern  double  *fx, *fy;                            /* Force arrays -------- */
    extern  NEIB    *neighbour;                          /* Paricle Neighbour --- */
    /* ---------------------------------- */


    /* ------------ Imp. variable ------------ */
    extern  double   sigma2;                             /* sigma square -----*/
    extern  double   sigma_by_r2;                        /* (sigma/dr)^2 -----*/
    extern  double   sigma_by_r6;                        /* (sigma/dr)^6 -----*/
    extern  double   C0, C2;                             /* LJ modifier ------*/
    extern  double   L;                                  /* Length of box ----*/
    extern  double   T_inst;                             /* Instaneous temp.  */         
    extern  double   PE;                                 /* Potential Energy -*/
    extern  double   KE;                                 /* Kinetic Energy ---*/
    /* --------------------------------------- */


#endif
