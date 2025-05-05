#ifndef GLOBAL_H
    #define GLOBAL_H
    

    /* -------- Libaries -------- */
    #include<stdio.h>
    #include<math.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include <unistd.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    #include <errno.h>
    /* -------------------------- */


    /* -------- System Variable -------- */
    #define     N              1000                      /* Total Number of Particles */
    #define     NA             (int) (0.80 * N)          /* Number of A type Particle */
    #define     NB             (int) (0.20 * N)          /* Number of B type Particle */
    #define     max_step       (int) 1e6               /* Step upto Simulation runs */
    #define     rho            1.20                      /* Density of the system --- */
    #define     vol            (double) (N / rho)        /* volume of the system ---- */
    #define     dt             0.005                     /* Time step for integration */
    #define     half_dt        (0.50) * dt               /* half time step for int -- */
    #define     init_temp      0.50                      /* Initial Temp. of system - */
    /* --------------------------------- */


    /* -------- For LJ Potential -------- */
    #define     EPSILON        1.0                       /* Interation strength of LJ Pot */                 
    #define     SIGMA          1.0                       /* effective radius of particles */
    #define     r_c            2.5                       /* Cut-Off for LJ Potential ---- */
    #define     r_c2           (r_c * r_c)               /* Cut-Off square -------------- */
    #define     sigma_by_rc    (SIGMA/r_c)               /* Sigma to cutoff ratio ------- */

    extern  double      sigma[2][2];                     /* For Binary we have different  */
    extern  double      epsilon[2][2];                   /* types of sigma and epsilon -- */
    /* ---------------------------------- */


    /* -------- For Verlet list -------- */
    #define     max_neigh      200                       /* strictly more than int(PI * rl^3 * rho) */
    #define     skin_depth     0.30 * SIGMA              /* Skin Width for Verlet List ------------ */
    #define     skin_depth2    (skin_depth * skin_depth) /* Skin Width square --------------------- */
    #define     rl             (r_c + skin_depth)        /* Radius for NBR list (r_c + skin) -----  */

    typedef struct {                                     /* NEIB data type has two branch */
        int list[max_neigh];                             /* one array branch  */
        int max_counter;                                 /* one scalar branch */
    }NEIB;

    extern int     list_counter;                         /* Counter to control after how many steps list will update */ 
    extern double  max_disp2;                            /* Max. Displacement of a particle in a step ---------------*/
    /* --------------------------------- */



    /* ------------ For Cell list ------------ */
    extern  int         n_cellx;                         /* Numbers of cell in x direction */
    extern  int         n_celly;                         /* Numbers of cell in y direction */
    extern  int         n_cellz;                         /* Numbers of cell in z direction */               
    extern  int         N_cell;                          /* Total numbers of cell of sys.  */
    extern  double      cell_cut;                        /* cell size we want ------------ */
    extern  double      cell_size;                       /* true size of cell ------------ */
    extern  int         *cell_map[13];                   /* 13 neigh cells for 3d -------- */
    /* --------------------------------------- */


    /* ---------- Nose-Hoover Stuff ---------- */
    #define     M_chain     3                            /* M chain Nose-Hoover -------*/
    #define     DOF         (int) (3 * N)                /* DOF of the system ---------*/
    #define     T_res       1.10                         /* Reservoid Temperature -----*/ 
    #define     tau         0.01                         /* Time constant of NH -------*/ 

    extern  double      eta[M_chain], p_eta[M_chain];    /* Extended dimensions -------*/
    extern  double      Q[M_chain];                      /* Effective mass ------------*/
    extern  double      H_eff;                           /* Effective Hamiltonian -----*/
    /* --------------------------------------- */


    /* ------- particle variables ------- */
    extern  double  *x, *y, *z;                          /* Position arrays ----- */
    extern  double  *xu, *yu, *zu;                       /* Unwrapped Position --- */
    extern  double  *vx, *vy, *vz;                       /* Velocity arrays ----- */
    extern  double  *fx, *fy, *fz;                       /* Force arrays -------- */
    extern  NEIB    *neighbour;                          /* Paricle Neighbour --- */
    extern  int     *pt_kind;                            /* Type of the particle  */
    /* ---------------------------------- */


    /* ------------ Imp. variable ------------ */
    extern  double   sigma2;                             /* SIGMA square -----*/
    extern  double   sigma_by_r2;                        /* (sig_ij/dr)^2 ----*/
    extern  double   sigma_by_r6;                        /* (sig_ij/dr)^6 ----*/
    extern  double   C0, C2;                             /* LJ Modifier ------*/
    extern  double   L, half_L;                          /* Length of box ----*/
    extern  double   T_inst;                             /* Instaneous temp.  */         
    extern  double   PE;                                 /* Potential Energy -*/
    extern  double   KE;                                 /* Kinetic Energy ---*/
    extern  double   pressure, virial;                   /* Pressure terms ---*/
    extern  int     *save_time;                      /* For saving data - */
    /* --------------------------------------- */


#endif
