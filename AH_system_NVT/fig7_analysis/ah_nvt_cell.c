/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* NVT Simulation of AH system with PBC */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ----- Libaries ----- */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
/* -------------------- */


/* ---- System specification ----- */
#define     N           2000                           /* Total Number of Paricles ----- */
#define     rho         0.20                           /* density of system ------------ */
#define     L           sqrt((double) N / rho)         /* Lenght of box ---------------- */
#define     L_by2       0.50 * L                       /* Half Length of box ----------- */

const int    max_step = (int) 1e4;                     /* No. of steps simulation run -- */
const double T_res = 1.0;                              /* Reserver Temperature --------- */
const double dt = 0.001;                               /* time step for integration ---- */
const double half_dt = 0.5 * dt;                       /* half time step for ----------- */
const double K = 1.0;                                  /* Velocity spin coupling ------- */
const double J = 1.0;                                  /* NN coupling ------------------ */
const double h = 0.00001;                              /* h of Step function ----------- */
/* ------------------------------- */


/* -------- Potential -------- */
const double epsilon = 1.0;
const double sigma = 1.0;
const double sigma2 = 1.0;
const double r_wca = 1.122462048 * sigma;     /* Cut-off for WCA potential */
const double r_c = 1.5 * sigma;               /* Cut-off for Net potential */
/* -------------------------- */


/* ---- Particles variables ---- */
double *x, *y, *theta;
double *px, *py, *omega;
double *fx, *fy, *omega_dot;
/* ---------------------------- */


/* -------- For Verlet list -------- */
#define max_neigh  30                      
const double skin_depth = 0.3 * sigma;
const double rl = (r_c + skin_depth);
int list_counter;                           /* counter to control after how many steps list will update */
typedef struct
{
    int list[max_neigh];
    int max_counter;
}NEIB;

NEIB *neighbour;
/* --------------------------------- */


/* ------------ For Cell list ----------- */
#define max_part_in_cell 30
const double cell_cut = (r_c + skin_depth);     /* size of a cell we want -------------------------- */
int n_cellx, n_celly, N_cell;                   /* numbers of cell --------------------------------- */
double cell_size;                               /* true size of cell ------------------------------- */
int *cell_map[4];                               /* 4 cells I choose for 2d ------------------------- */
typedef struct                                  /* list of particles in a cell with max. numbers --- */
{
    int list[max_part_in_cell];
    int max_counter;
}CLIST;    
/* -------------------------------------- */



/* ----- Useful variables ----- */
double sigma_by_r2, sigma_by_r6;
double p_s = 0.0, s = 1.0;                        /* --- Extended dimensions --- */
double T1, T2, T3, T;                             /* --- Temperature ----------- */
double virial, pressure;                          /* --- Pressure terms -------- */
double KE, PE, RE, E_s, E;                        /* --- Energies -------------- */
double H_nose0;                                   /* --- H_Nose at t=0 --------- */
const double Q = 1.0;                             /* --- Extended mass --------- */
const double Intertia = 1.0;                      /* --- Moment of Intertia ---- */
const double pi = 3.14159265359;
int run = 0;                                      /* --- for ensemble ---------- */
/* ---------------------------- */



// Allocating Arrays for Use
void allocate_arrays()
{
    x     = malloc(sizeof(double) * N);
    y     = malloc(sizeof(double) * N);
    theta = malloc(sizeof(double) * N);

    if(x == NULL || y == NULL || theta == NULL)
    {
        printf("Allocation request denied for positions. \n");
    }


    px    = malloc(sizeof(double) * N);
    py    = malloc(sizeof(double) * N);
    omega = malloc(sizeof(double) * N);
    if(px == NULL || py == NULL || omega == NULL)
    {
        printf("Allocation request denied for velocity. \n");
    }


    fx        = malloc(sizeof(double) * N);
    fy        = malloc(sizeof(double) * N);
    omega_dot = malloc(sizeof(double) * N);
    if(fx == NULL || fy == NULL || omega_dot == NULL)
    {
        printf("Allocation request denied for Forces. \n");
    }

    neighbour = (NEIB *)malloc(sizeof(NEIB) * N);
    if(neighbour == NULL) printf("Allocation request denied for neighbour. \n");

    
    for (int i = 0; i < 4; i++)
    {
        cell_map[i] = malloc(sizeof(int) * N_cell );
        if (cell_map[i] == NULL) printf("Allocation request denied for cell_map \n");
    }
}


/* --- Load initial position --- */
void load_initial_position()
{ 
    // No of particle in x, y direction
    int n = (int) ceil( sqrt((double) N) );
        
    // Distances b/w particles
    // double a = L / (double) n;
    double a = 1.20;

    int pt_num = 0;

    int exit_flag = 0;

    for(int i=0; i<n; i++)             // x-loop
     {
        for(int j=0; j<n; j++)         // y-loop
         {
            if(pt_num > N-1)
            {
                exit_flag = 1;
                break;
            }
            x[pt_num] = i * a + 3.5 * a;
            y[pt_num] = j * a + 3.5 * a;

            pt_num++ ;
        }
        if(exit_flag) break;
    }

    if(pt_num < N)
    {
        printf("Warning ! All particles are NOT placed. \n");
    } 

}


/* ---- Make Centered hexagonal Lattice ---- */
void load_triangular_lattice()
{
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
    * H_n = nth centered hexagonal number            *
    * n = number requried for getting H_n            *
    * source: https://stackoverflow.com/a/14283364   *
    *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    int n, H_n = N;
    n = ceil(0.5 + sqrt(12*H_n-3)/6.0); 

    int pt_num = 0, exit_flag = 0;
    double x_pos, y_pos;
    double sqrt3_by2 = sqrt(3)/2.0;
    double a = 1.123;

    for(int i=0; i < n; i++)
    {
        y_pos = sqrt3_by2 * a * i;
        for (int j = 0; j < (2*n-1-i); j++)
        {
            x_pos = (-(2*n-i-2)*a)/2.0 + j*a;

            if (pt_num > N-1){
                exit_flag = 1;
                break;
            }

            x[pt_num] = x_pos;
            y[pt_num++] = y_pos;
            
            if (y_pos != 0)
            {
                if (pt_num > N-1){
                    exit_flag = 1;
                    break;
                }
                x[pt_num] = x_pos;
                y[pt_num++] = -y_pos;
            }
        }
        if(exit_flag) break;
    }

    // Bring to center of box
    for (int i = 0; i < N; i++)
    {
        x[i] = x[i] + L_by2;
        y[i] = y[i] + L_by2;
    }
    
    if(pt_num < N)
    printf("Warning ! All particles are NOT placed. \n");
}


/* --- Load initial velocity --- */
void load_initial_velocity()
{
    double sumpx = 0.0, sumpy = 0.0;

    // Initize the seed 
    srand48(time(NULL));

    for(int i=0; i<N; i++)
    {
        px[i] = 2.0 * drand48() - 1.0;
        py[i] = 2.0 * drand48() - 1.0;

        sumpx = sumpx + px[i];
        sumpy = sumpy + py[i];
    }

    sumpx = sumpx / (double) N;
    sumpy = sumpy / (double) N;

    for(int i=0; i<N; i++)
    {
        px[i] = px[i] - sumpx;
        py[i] = py[i] - sumpy;
    }

}


/* --- Load initial spins --- */
void load_initial_spin()
{
    double sum_omega = 0.0;
    for (int i = 0; i < N; i++)
    {
        theta[i] = 0.0;
        omega[i] = 2.0 * drand48() - 1.0;
        sum_omega = sum_omega + omega[i];
    }
    sum_omega = sum_omega / (double) N;

    for (int i = 0; i < N; i++) omega[i] -= sum_omega;
}


/* --- Load final config of a previous run --- */
void load_final_config()
{
    FILE *fid_load;

    fid_load = fopen("final_config", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config.txt\n");
        exit(1);
    }

    fscanf(fid_load, "%lf %lf", &s, &p_s);
    for ( int i = 0; i < N; i++ )
    {
        fscanf(fid_load, "%lf %lf %lf %lf %lf %lf",
                &x[i], &y[i], &theta[i], &px[i], &py[i], &omega[i]);
    }

    fclose(fid_load);
}



/* mode:1 leads to make a new config, * 
 * mode:0 leads to read old cofigure  */
void ProgramInitializer(int mode)
{
    n_cellx = n_celly = (int) floor(L / cell_cut);
    N_cell = n_cellx * n_celly;
    cell_size = L / (double)n_cellx;   /* (cell_size > cell_cut) */

    allocate_arrays();
    printf("Arrays are allocated. \n");
    printf(" -------------------------- \n");

    if (mode)
    {
        // load_initial_position();
        load_triangular_lattice();
        printf("Initial position is loaded. \n");
        printf(" ---------------------------- \n");
    
        load_initial_velocity();
        printf("Initial Momentum is given. \n");
        printf(" ---------------------------- \n");
    
        load_initial_spin();
        printf("Initial Spin and Ang. Momentum is given. \n");
        printf(" --------------------------------------- \n");
    }
    else
    {
        load_final_config();
        printf("Final configaration of previous run reads successfully. \n");
        printf(" ------------------------------------------------------ \n");
    }
}


/* ---- Calculate KE & Rotational Energy of the system ---- */
/* --- Calculate temperature of system at current state --- */
void Temperature()
{
    double sump2 = 0.0, sum_omega2 = 0.0;
    double sum_pA = 0.0;
    double temp, s_inv, s_inv2;

    s_inv = 1.0 / s;
    s_inv2 = s_inv * s_inv;

    for (int i = 0; i < N; i++)
    {
        sump2      += px[i]*px[i] + py[i]*py[i];
        sum_omega2 += omega[i]*omega[i];
        sum_pA     += px[i]*cos(theta[i]) + py[i]*sin(theta[i]);
    }
    
    // T1 = \frac{1}{3NK_B} [ p_i^2/s^2]
    T1 = s_inv2 * sump2 / (3.0 * N);

    // T2 = - \frac{1}{3NK_B} p.A/m 
    T2 = - s_inv * K * sum_pA / (3.0 * N);

    // T3 = \frac{1}{3NK_B} [ w^2/I]
    T3 = s_inv2 * sum_omega2 / ( 3.0 * Intertia * N); 
    
    T = T1 + T2 + T3;

    // KE = sum [ p_i^2 / 2ms^2 - p.A/m ] + NK^2/2m
    KE = 0.50 * s_inv2 * sump2 - K * s_inv * sum_pA + 0.5 * (N * K * K);

    // Rotational Energy = w^2/2I
    RE = 0.50 * s_inv2 * sum_omega2 / Intertia;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * calculate distance between i and j particles   *
 * Input: Particle index i, j                     *
 * Output: Distance b/w the particles.            *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double particle_distance(int i, int j)
{
    double dx, dy, dr2;
    
    dx = x[i] - x[j];
    dy = y[i] - y[j];

    // minimum image condition
    if (dx > L_by2) dx = dx - L;
    else if (dx < -L_by2) dx = dx + L;

    if (dy > L_by2) dy = dy - L;
    else if (dy < -L_by2) dy = dy + L;

    dr2 = dx*dx + dy*dy;

    return dr2;
}


/* ------ give cell index for (i,j) pair box ------- */
int getcell_index(int i, int j)
{
    // apply pbc box
    if (i < 0) i = i + n_cellx;
    else if (i >= n_cellx) i = i - n_cellx;

    if (j < 0) j = j + n_celly;
    else if (j >= n_celly) j = j - n_celly;

    int cellno = i + j * n_cellx;
    return cellno;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Give appropiate neighbours(L shape) cell num for cell list *
 * For a sample output for n_cellx = 3, see comment at end of *
 * the program.                                               * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void maps()
{
    int imap;

    for (int i = 0; i < n_cellx; i++)
    {
        for (int j = 0; j < n_celly; j++)
        {
            imap = getcell_index(i, j);
            
            cell_map[0][imap] = getcell_index(i+1, j);      // east(1,0)
            cell_map[1][imap] = getcell_index(i+1, j+1);    // north-east(1,1)
            cell_map[2][imap] = getcell_index(i, j+1);      // north(0,1)
            cell_map[3][imap] = getcell_index(i-1, j+1);    // north-west(-1,1)
        }
        
    }
    
}


/* ------- construct neighbour list using cell list ------- */
void construct_cell_list()
{
    int cell_index, L_index, k;         // neigh cell index L type 
    int nbc;                            // neighbour cell (1 to 4)
    int par_index1, par_index2;         // particle index

    double dr2;

    CLIST cell_list[N_cell];


    // restart counter from 0
    list_counter = 0;

    // --> Initialize Neighbour list <--
    for (int i = 0; i < N; i++)
    {
        neighbour[i].max_counter = 0;
        for (int j = 0; j < max_neigh; j++)
        {
            neighbour[i].list[j] = 0;
        }
    }

    // --> Initialize cell list <--
    for (int i = 0; i < N_cell; i++)
    {
        cell_list[i].max_counter = 0;
        for (int j = 0; j < max_part_in_cell; j++)
        {
            cell_list[i].list[j] = 0;
        }
    }


    /* ~~~~~~~~~ distribute particles to cell ~~~~~~~~~ */
    for (int i = 0; i < N; i++)
    {
        cell_index = floor( x[i]/cell_size ) + floor ( y[i]/cell_size ) *  n_cellx;
        /* cell_index(k) = i + j * n_cellx */

        k = cell_index;
        
        // Debug
        if ( k < 0 || k >= N_cell)
        {
            printf("Something went wrong. \n");
            exit(0);
        }
        

        cell_list[k].list[cell_list[k].max_counter] = i;
        cell_list[k].max_counter ++;
    }


    /* ~~~~~ generate neighbour list using cell structure ~~~~~ */
    for (int i = 0; i < N_cell; i++)                                
    {
        k = cell_list[i].max_counter;       // max particle in a particular cell

        /****** Intra cell neighbour *******/
        for (int j = 0; j < k-1; j++)
        {
            par_index1 = cell_list[i].list[j];
        
            for (int m = j+1; m < k; m++)       
            {
                par_index2 = cell_list[i].list[m];
                dr2 = particle_distance(par_index1, par_index2);
                if (dr2 < rl * rl)
                {
                    neighbour[par_index1].list[neighbour[par_index1].max_counter] = par_index2;
                    neighbour[par_index1].max_counter++ ;
                }
                
            }
        }


        /****** Inter cell neighbour *******/
        for (int j = 0; j < k; j++)
        {
            par_index1 = cell_list[i].list[j];

            for (nbc = 0; nbc < 4; nbc++)       // 4 neighbour box for 2d
            {
                L_index = cell_map[nbc][i];     // select index of a neighbour box

                for (int m = 0; m < cell_list[L_index].max_counter; m++)
                { 
                    // select another particles of that neighbour box
                    par_index2 = cell_list[L_index].list[m];        

                    dr2 = particle_distance(par_index1, par_index2);
                    if (dr2 < rl * rl)
                    {
                        neighbour[par_index1].list[neighbour[par_index1].max_counter] = par_index2;
                        neighbour[par_index1].max_counter ++;
                    } 
                }
            }
            
        }  
        
    }
    
}


void Cell_Initializer()
{

    maps();
    printf("Cell divided successfully. \n");
    printf(" --------------------------------------- \n");

    construct_cell_list();
    printf("Cell list construct sucessfully. \n");
    printf(" --------------------------------------- \n");
}


/* ---- Function for select nearest neighboours ---- */
double step_function(double dr)
{
    double t = dr - r_c;
    double t2 = t * t;
    double t4 = t2 * t2;
    return (t4/(h + t4));
}


void force_calculation()
{
    double dx, dy, theta_ij;
    double dr, dr2;
    double t0, f0, grij;
    double r_dist;
    int j, tot_neigh;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    for (int i = 0; i < N; omega_dot[i] = 0.0, i++);
    
    PE = 0.0;        /* PE initialize */
    virial = 0.0;    /* Pressure initilize */

    for (int i = 0; i < N; i++)
    {
        tot_neigh = neighbour[i].max_counter;
        for (int k = 0; k < tot_neigh; k++)
        {
            j = neighbour[i].list[k];
            dx = x[i] - x[j];
            dy = y[i] - y[j];

            // Minimum Image condition
            if (dx > L_by2) dx = dx - L;
            else if (dx < - L_by2) dx = dx + L;

            if (dy > L_by2) dy = dy - L;
            else if (dy < - L_by2) dy = dy + L;

            dr2 = (dx*dx + dy*dy);

            if (dr2 < sigma2/100.0)
            {
                printf("Is this physical? %.14lf \n", dr2);
                exit(1);
            }
            
            if (dr2 < r_c*r_c)
            {
                dr = sqrt(dr2);

                /* ------------- For Omega_dot ------------- */
                theta_ij = theta[i] - theta[j];
                grij = step_function(dr);
                t0 = - J * s * grij * sin(theta_ij);
                omega_dot[i] = omega_dot[i] + t0;
                omega_dot[j] = omega_dot[j] - t0;
                /* ----------------------------------------- */

                
                /* -------------- For fx and fy -------------- */
                r_dist = pow((dr - r_c), 5);
                t0 = (grij * grij) / (dr * r_dist);
                t0 = t0 * (4 * J * s * h) * cos(theta_ij);

                fx[i] = fx[i] + t0 * dx;
                fy[i] = fy[i] + t0 * dy;

                fx[j] = fx[j] - t0 * dx;
                fy[j] = fy[j] - t0 * dy;
                /* ------------------------------------------- */


                PE     -= J * grij * cos(theta_ij);
                virial += t0 * dr2 / s;

                if (dr < r_wca )
                {
                    sigma_by_r2 = sigma2/dr2;
                    sigma_by_r6 = sigma_by_r2 * sigma_by_r2 * sigma_by_r2;
                    f0 = 48 * epsilon * sigma_by_r6 * (sigma_by_r6 - 0.5 ) / dr2;
                    f0 = f0 * s;
                    
                    fx[i] = fx[i] + f0 * dx;
                    fy[i] = fy[i] + f0 * dy;

                    fx[j] = fx[j] - f0 * dx;
                    fy[j] = fy[j] - f0 * dy;

                    PE     += 4.0 * epsilon * sigma_by_r6 * (sigma_by_r6 - 1.0) + epsilon;
                    virial += f0 * dr2 / s; 
                }

            }

        }

    }

    // debug
    // for (int i = 0; i < N; i++)
    // {
    //     printf("fx[i] = %.12lf\n", fx[i]);
    // }
    

}


/* ---- Update Torque ----*/
void update_torque()
{
    for (int i = 0; i < N; i++)
    {
        omega_dot[i] += K * (py[i]*cos(theta[i]) - px[i]*sin(theta[i]));
    }
    
}


/*--- Canonical Momentum Updatation (omega) ---*/
void update_cano_ang_momenta()
{
    for (int i = 0; i < N; i++)
    {
        omega[i] += omega_dot[i] * half_dt;
    }
    
}


/*--- Canonical Momentum Updatation (px, py) ---*/
void update_cano_lin_momenta()
{
    for (int i = 0; i < N; i++)
    {
        px[i] = px[i] + fx[i] * half_dt;
        py[i] = py[i] + fy[i] * half_dt;
    }

}


/*----- Update p_s step: 1 -----*/
void update_ps_step1()
{
    double M1, M2, M3, M;
    M1 = M2 = M3 = M = 0.0;

    double dx, dy, theta_ij, dr2;

    for (int i = 0; i < N; i++)
    {
        M1 += px[i]*px[i] + py[i]*py[i] + (omega[i] * omega[i]) / Intertia ; 
    }

    M1 = M1 / (2 * s * s) - 0.50 * (N * K * K);

    M2 = - PE;

    M3 = - 3.0 * N * T_res * (1.0 + log(s));
    
    M = M1 + M2 + M3 + H_nose0;

    p_s = p_s + half_dt * M;
    
}


/* ----- Update p_s step: 2 ----- */
void update_ps_step2()
{
    double tt = 1.0 + (0.5 * half_dt / Q) * p_s;    /* 0.5*(half_dt) = (0.25)*dt */
    p_s = p_s / tt;
}


/* ----- Update s ----- */
void update_s()
{
    s = s * exp( half_dt * p_s / Q);
}


/* ----- Update theta ----- */
void update_theta()
{
    double s_inv = 1.0 / s;

    for (int i = 0; i < N; i++)
    {
        theta[i] += omega[i] * half_dt * s_inv / Intertia;

        // wrap angle b/w [0, 2pi)
        if (theta[i] >= 2 * pi) theta[i] -= 2.0 * pi;
        else if (theta[i] < 0.0) theta[i] += 2.0 * pi;
    }
    
}



void integrate()
{
    double dx, dy, dr2;
    double s_inv, max_disp2 = 0.0;

    /* -------- half-kick --------- */
    update_cano_ang_momenta();
    update_cano_lin_momenta();
    update_ps_step1();
    update_ps_step2();
    update_s();
    update_theta();
    /* ---------------------------- */

    s_inv = 1.0 / s;

    for (int i = 0; i < N; i++)
    {
        dx = ( px[i] * s_inv - K * cos(theta[i]) ) * dt;
        x[i] = x[i] + dx;
        if (x[i] >= L) x[i] = x[i] - L;
        else if (x[i] < 0) x[i] = x[i] + L;

        dy = ( py[i] * s_inv - K * sin(theta[i]) ) * dt;
        y[i] = y[i] + dy;
        if (y[i] >= L) y[i] = y[i] - L;
        else if (y[i] < 0) y[i] = y[i] + L;

        dr2 = dx*dx + dy*dy;
        if (dr2 > max_disp2) max_disp2 = dr2;
    }

    list_counter ++;
    if ((2.0 * sqrt(max_disp2) * list_counter) > skin_depth)
    {
        // printf("Clist constructed after step: %d \n", list_counter);
        construct_cell_list();
    }
    

    /* -------- end-kick --------- */
    update_theta();
    update_s();
    force_calculation();        // Recalculate forces 
    update_ps_step2();
    update_ps_step1();
    force_calculation();
    update_cano_lin_momenta();
    update_torque();            // Recalculate omega_dot
    update_cano_ang_momenta();
    /* ---------------------------- */

}



void save_data_infile(int step)
{
    char filename[100];
    sprintf(filename, "movie_data/step_%06d.txt", step); // "step_0001.txt", "step_0002.txt", ...

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < N; i++)
    fprintf(file, "%.8lf %.8lf %.8lf\n", x[i], y[i], theta[i]);


    fclose(file);
}


void show_progress (int current)
{
    int percentage = (current * 100) / max_step;
    
    // \r returns cursor to start of line
    printf("\rProgram completed: %d%%", percentage);
    
    // Ensure immediate output
    fflush(stdout);
}


/* Rem. to free memory at end */
void deallocate_arrays()
{
    free(x); free(y); free(theta);

    free(px); free(py); free(omega);

    free(fx); free(fy); free(omega_dot);

    free(neighbour);

    for (int i = 0; i < 4; i++) free(cell_map[i]);

}


int main()
{
    ProgramInitializer(1); // switch to control , 1 (new), 0 (old)

    Cell_Initializer();

    Temperature();
    printf("Temperature of system : %f \n", T);
    printf(" --------------------------------------- \n");
    
    force_calculation();
    printf("Initial Force calculation is done. \n");
    printf(" --------------------------------------- \n");
    
    update_torque();

    E_s = 0.5 * (p_s*p_s) / Q + 3.0 * N * T_res * log(s);
    H_nose0 = KE + RE + PE + E_s;

    // Equilibriate the system
    for (int i = 0; i < (int)1e5; i++)
    {
        // debug
        double sumpx = 0.0, sumpy = 0.0, sump;
        for (int i = 0; i < N; i++) sumpx += px[i];
        for (int i = 0; i < N; i++) sumpy += py[i];
        sump = sumpx*sumpx + sumpy*sumpy;
        if ( sump > 1e-12){
            printf("\nmistake detected = %.12lf\n", sqrt(sump));
            exit(1);
        }

        integrate();
    }
    
    // save data from here
    for (int step = 0; step < max_step; step++)
    {
        // debug
        double sumpx = 0.0, sumpy = 0.0, sump;
        for (int i = 0; i < N; i++) sumpx += px[i];
        for (int i = 0; i < N; i++) sumpy += py[i];
        sump = sumpx*sumpx + sumpy*sumpy;
        if ( sump > 1e-12){
            printf("\nmistake detected = %.12lf\n", sqrt(sump));
            exit(1);
        }

        Temperature();
        pressure = rho * T + 0.5 * virial / (L*L);

        E_s = 0.5 * (p_s*p_s) / Q + 3.0 * N * T_res * log(s);   
        E = KE + RE + PE;
    
        integrate();
    }

    /* ------------ Saving avg -------------- */
    char make_dir[128], filename4[256];
    sprintf(make_dir, "mkdir -p velocity_data/run_%02d", run);
    system(make_dir);
    sprintf(filename4, "velocity_data/run_%d/N_%d_T_%.1f", run, N, T_res);
    FILE *fid4;
    fid4 = fopen(filename4, "w");
    if (!fid4) printf("Error! while opening fid4");

    double vx, vy;
    double s_inv = 1.0 / s;

    for (int i = 0; i < N; i++)
    {
        vx = px[i] * s_inv - K * cos(theta[i]);
        vy = py[i] * s_inv - K * sin(theta[i]);

        fprintf(fid4, "%.12lf\t %.12lf\n", vx, vy);
    }
    fclose(fid4);
    /* -------------------------------------- */

    deallocate_arrays();
    
    printf("\n ===========================================\n");
    printf("    ✅ SUCCESS! Everything ran smoothly.     \n");
    printf(" ===========================================\n\n");
    

    return 0;
}
