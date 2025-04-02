/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* NVE Simulation of AH system with PBC */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ----- Libaries ----- */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
/* -------------------- */

 
/* ---- System specification ----- */
#define     N           128                            /* Number of particles ---------- */
#define     rho         0.20                           /* density of system ------------ */
#define     L           sqrt((double) N / rho)         /* Lenght of box ---------------- */

const double dt = 0.001;                               /* time step for integration ---- */
const int    max_step = (int) 1e6;                     /* No. of steps simulation run -- */
const double K = 1.0;                                  /* spin velocity coupling ------- */
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
int list_counter;                           // counter to control after how many steps list will update
typedef struct
{
    int list[max_neigh];
    int max_counter;
}NEIB;

NEIB *neighbour;
/* --------------------------------- */


/* ------------ For Cell list ----------- */
#define max_part_in_cell 30
const double cell_cut = (r_c + skin_depth);     // size of a cell we want
int n_cellx, n_celly, N_cell;                   // numbers of cell
double cell_size;                               // true size of cell
int *cell_map[4];                               // 4 cells I choose for 2d
typedef struct                                  // list of particles in a cell with max. numbers
{
    int list[max_part_in_cell];
    int max_counter;
}CLIST;    
/* -------------------------------------- */


/* ----- Useful variables ----- */
double sigma_by_r2, sigma_by_r6;
double T1, T2, T3, T;                              /* --- Temperature ------ */
double KE, PE, RE, E;                              /* Energies ------------- */
const double Intertia = 1.0;                       /* Moment of Intertia --- */
const double pi = 3.14159265359;
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


/* --- Load final config of a previous run --- */
/* --- Also tranform from prime to normal ---- */
void load_final_config()
{
    double s, p_s, s_inv;

    FILE *fid_load;

    fid_load = fopen("final_config_NVT", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config.txt");
        exit(1);
    }

    fscanf(fid_load, "%lf %lf", &s, &p_s);
    s_inv = 1 / s;

    for ( int i = 0; i < N; i++ )
    {
        fscanf(fid_load, "%lf %lf %lf %lf %lf %lf",
        &x[i], &y[i], &theta[i], &px[i], &py[i], &omega[i]);

        px[i] = px[i] * s_inv;
        py[i] = py[i] * s_inv;
        omega[i] = omega[i] * s_inv;
    }
    

    fclose(fid_load);
}



/* ---- Calculate KE & Rotational Energy of the system ---- */
/* --- Calculate temperature of system at current state --- */
void Temperature()
{
    double sump2 = 0.0, sum_omega2 = 0.0;
    double sum_pA = 0.0;

    for (int i = 0; i < N; i++)
    {
        sump2      += px[i]*px[i] + py[i]*py[i];
        sum_omega2 += omega[i]*omega[i];
        sum_pA     += px[i]*cos(theta[i]) + py[i]*sin(theta[i]);
    }
    
    // T1 = \frac{1}{3NK_B} [ p_i^2/s^2]
    T1 = sump2 / (3.0 * N);

    // T2 = - \frac{1}{3NK_B} p.A/m 
    T2 = - K * sum_pA / (3.0 * N);

    // T3 = \frac{1}{3NK_B} [ w^2/I]
    T3 = sum_omega2 / ( 3.0 * Intertia * N); 
    
    T = T1 + T2 + T3;

    // KE = sum [ p_i^2 / 2m - p.A/m + w^2/2I] + NK^2/2m
    KE = 0.50 * sump2 - K * sum_pA + 0.5 * (N * K * K);
    KE = KE / (double) N;

    RE = 0.50 * sum_omega2 / Intertia;
    RE = RE / (double) N;

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
    if (dx > 0.50 * L) dx = dx - L;
    else if (dx < -0.50 * L) dx = dx + L;

    if (dy > 0.50 * L) dy = dy - L;
    else if (dy < -0.50 * L) dy = dy + L;

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


/* --------- construct neighbour list using cell list --------- */
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

        if (abs(k) >= N_cell)
        {
            printf("Break at debug point 1 \n");
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


/* ---- Function for select nearest neighboours ---- */
double step_function(double dr)
{
    double t = dr - r_c;
    double t2 = t * t;
    double t4 = t2 * t2;
    return t4/(h + t4);
}


/* ---- Calculate force using cell lists ---- */
void force_calculation()
{
    double dx, dy, theta_ij;
    double dr, dr2;
    double t0, f0, int_strength;
    double r_dist;
    int j, tot_neigh;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    for (int i = 0; i < N; omega_dot[i] = 0.0, i++);
    
    PE = 0.0;                         /* PE initialize */

    for (int i = 0; i < N; i++)
    {
        tot_neigh = neighbour[i].max_counter;
        for (int k = 0; k < tot_neigh; k++)
        {
            j = neighbour[i].list[k];
            dx = x[i] - x[j];
            dy = y[i] - y[j];

            // Minimum Image condition
            if (dx > 0.5 * L) dx = dx - L;
            else if (dx < - 0.5 * L) dx = dx + L;

            if (dy > 0.5 * L) dy = dy - L;
            else if (dy < - 0.5 * L) dy = dy + L;

            dr2 = (dx*dx + dy*dy);
            
            if (dr2 < sigma2/100.0)
            {
                printf("Is this physical? %f \n", dr2);
                exit(1);
            }

            
            if (dr2 < (r_c * r_c) )
            {
                dr = sqrt(dr2);
                theta_ij = theta[i] - theta[j];
                int_strength = step_function(dr);
                t0 = - J * int_strength * sin(theta_ij);
                omega_dot[i] = omega_dot[i] + t0;
                omega_dot[j] = omega_dot[j] - t0;

                r_dist = pow((dr - r_c), 5);
                t0 = (int_strength * int_strength * cos(theta_ij))/(dr * r_dist);
                t0 = t0 * (4 * J * h);
                fx[i] = fx[i] + t0 * dx;
                fy[i] = fy[i] + t0 * dy;

                fx[j] = fx[j] - t0 * dx;
                fy[j] = fy[j] - t0 * dy;

                PE += - J * int_strength * cos(theta_ij);

                if ( dr < r_wca )
                {
                    sigma_by_r2 = sigma2/dr2;
                    sigma_by_r6 = sigma_by_r2 * sigma_by_r2 * sigma_by_r2;
                    f0 = 48 * epsilon * sigma_by_r6 * (sigma_by_r6 - 0.5 ) / dr2;
                    
                    fx[i] = fx[i] + f0 * dx;
                    fy[i] = fy[i] + f0 * dy;

                    fx[j] = fx[j] - f0 * dx;
                    fy[j] = fy[j] - f0 * dy;

                    PE += 4.0 * epsilon * sigma_by_r6 * (sigma_by_r6 - 1.0) + epsilon;
                }

            }

        }
        
    }
    PE = PE / (double) N;

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
        omega[i] += omega_dot[i] * 0.50 * dt;
    }
    
}


/*--- Canonical Momentum Updatation (px, py) ---*/
void update_cano_lin_momenta()
{
    for (int i = 0; i < N; i++)
    {
        px[i] = px[i] + fx[i] * 0.50 * dt;
        py[i] = py[i] + fy[i] * 0.50 * dt;
    }

}



/*----- Update theta -----*/
void update_theta()
{
    for (int i = 0; i < N; i++)
    {
        theta[i] = theta[i] + omega[i] * 0.5 * dt / Intertia;

        // wrap angle b/w [0, 2pi)
        if (theta[i] >= 2 * pi) theta[i] -= 2.0 * pi;
        else if (theta[i] < 0.0) theta[i] += 2.0 * pi;
    }
    
}


void integrate()
{
    double dx, dy;
    double dr2, max_disp2 = 0.0;

    /* --------- half-kick ---------- */
    update_cano_ang_momenta();
    update_cano_lin_momenta();
    update_theta();
    /* ------------------------------ */

    for (int i = 0; i < N; i++)
    {
        dx = (px[i] - K * cos(theta[i])) * dt;
        x[i] = x[i] + dx;
        if (x[i] >= L) x[i] = x[i] - L;
        else if (x[i] < 0) x[i] = x[i] + L;

        dy = (py[i] - K * sin(theta[i])) * dt;
        y[i] = y[i] + dy;
        if (y[i] >= L) y[i] = y[i] - L;
        else if (y[i] < 0) y[i] = y[i] + L;

        dr2 = dx*dx + dy*dy;
        if (dr2 > max_disp2) max_disp2 = dr2;
    }

    
    list_counter ++;
    if ( (2.0 * sqrt(max_disp2) * list_counter) > skin_depth )
    {        
        // printf("clist constructed after step: %d \n", list_counter);
        construct_cell_list();
    }
    

    /* --------- end-kick ----------- */
    update_theta();
    force_calculation();
    update_cano_lin_momenta();
    update_torque();
    update_cano_ang_momenta();
    /* ------------------------------ */
    
}


void save_data_infile(int step, double x[], double y[], double theta[])
{
    char filename[100];
    sprintf(filename, "movie_data/step_%06d.txt", step); // "step_0001.txt", "step_0002.txt", ...

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        fprintf(file, "%.12lf %.12lf %.12lf\n", x[i], y[i], theta[i]);
    }

    fclose(file);
}


void show_progress (int current, int total)
{
    int percentage = (current * 100) / total;
    
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
    /* For cell list */
    n_cellx = n_celly = (int) floor(L / cell_cut);
    N_cell = n_cellx * n_celly;
    cell_size = L / (double)n_cellx;   /* cell_size > cell_cut */
    

    allocate_arrays();
    printf("Arrays are allocated. \n");
    printf(" -------------------------- \n");
    

    load_final_config();
    printf("Final configaration of NVT reads successfully. \n");
    printf(" ------------------------------------------------------ \n");
    
    printf("Variable changes from prime to normal. \n");
    printf("------------------------------------------ \n");

    Temperature();
    printf("Temperature of system : %f \n", T);
    printf(" --------------------------------------- \n");

    maps();
    printf("Cell divided successfully. \n");
    printf(" --------------------------------------- \n");

    construct_cell_list();
    printf("Cell list construct sucessfully. \n");
    printf(" --------------------------------------- \n");

    force_calculation();
    printf("Initial Force calculation is done. \n");
    printf(" --------------------------------------- \n");
    
    update_torque();

    FILE *fid1, *fid2;
    fid1 = fopen("Energy_data_NVE_N128_rho1.20_K1.0_dt0.001", "w");
    fid2 = fopen("Temperature_NVE_N128_rho1.20_K1.0_dt0.001", "w");
    if ( !fid1 || !fid2 ) return 1;
 
    for (int step = 0; step < max_step; step++)
    {
        // if (step % 10 == 0) save_data_infile(step, x, y, theta);
        

        Temperature();
        // Save thermodynamics property
        if (step % 100 == 0)
        {
            fprintf(fid1, "%.12lf\t %.12lf\t %.12lf\t %.12lf\n ", 
                         KE, RE, PE, (KE+RE+PE) );
            fprintf(fid2, "%.12lf\t %.12lf\t %.12lf\t %.12lf\n", 
                    T1, T2, T3, T );
        }
        integrate();

        show_progress(step, max_step);

    }
    
    fclose(fid1);
    fclose(fid2);



    // Change to true to run
    if (0)
    {
        FILE *fid3;
        fid3 = fopen("final_config", "w");
        if( fid3 == NULL  ) return(1);
        for( int i = 0; i < N; i++ )
        {
            fprintf(fid3, "%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf \n",
                x[i], y[i], theta[i], px[i], py[i], omega[i]);
        }
        fclose(fid3); 
    }

    deallocate_arrays();

    printf("\n");
    return 0;
}










/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * This is NVE simulation of the Actice Hamiltonian system presented  *
 * in the paper. As it is NVE system we don't need s, p_s variables.  *
 * We only try to plot energy of real hamiltonian which should be     *
 * conserved in this case. Along with flactuation dt^2.               *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
