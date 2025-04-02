/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* NVT Simulation of AH system with PBC */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ----- Libaries ----- */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
/* -------------------- */


/* ----- System specification ----- */
#define     N           2                            /* Number of particles ---------- */
#define     rho         0.20                           /* density of system ------------ */
#define     L           sqrt((double) N / rho)         /* Lenght of box ---------------- */

const int    max_step = (int) 2;                     /* No. of steps simulation run -- */
const double T_res = 0.15;                             /* Relax Temperature for NVT ---- */
const double dt = 0.001;                               /* time step for integration ---- */
const double K = 0.0;                                  /* spin velocity coupling ------- */
const double J = 1.0;                                  /* NN coupling ------------------ */
const double h = 0.00001;                              /* h of Step function ----------- */
/* ------------------------------- */


/* -------- Potential -------- */
const double epsilon = 1.0;
const double sigma = 1.0;
const double sigma2 = 1.0;
const double r_wca = 1.122462048 * sigma;     /* Cut-off for WCA potential */
const double sigma_by_rw = sigma/r_wca;
const double r_c = 1.5 * sigma;               /* Cut-off for Net potential */
/* -------------------------- */


/* ---- Particles variables ---- */
double *x, *y, *theta;
double *px, *py, *omega;
double *fx, *fy, *omega_dot;
/* ---------------------------- */


/* ----- Useful variables ----- */
double sigma_by_r2, sigma_by_r6;
double p_s = 0.0, s = 1.0;                         /* Extended dimensions --- */
double T1, T2, T3, T;                              /* Temperature ----------- */
double KE, PE, RE, E_s, E;                         /* Energies -------------- */
double H_nose0;                                    /* H_Nose at t=0 --------- */
const double Q = 1.0;                              /* Extended mass --------- */
const double Intertia = 1.0;                       /* Moment of Intertia ---- */
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
}


/* --- Load initial position --- */
void load_initial_position()
{ 
    // No of particle in x, y direction
    int n = (int) ceil( sqrt((double) N) );
        
    // Distances b/w particles
    // double a = L / (double) n;
    double a = 1.0;

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
            x[pt_num] = i * a + 5 * a;
            y[pt_num] = j * a + 5 * a;

            pt_num++ ;
        }
        if(exit_flag) break;
    }

    if(pt_num < N)
    {
        printf("Warning ! All particles are NOT placed. \n");
    } 

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
        printf("Error: Unable to open file: final_config.txt");
        exit(1);
    }
    // fscanf(fid_load, "%lf %lf", &s, &p_s);
    for ( int i = 0; i < N; i++ )
    {
        fscanf(fid_load, "%lf %lf %lf %lf %lf %lf",
                &x[i], &y[i], &theta[i], &px[i], &py[i], &omega[i]);
    }

    fclose(fid_load);
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


/* ---- Function for select only nearest neighboours ---- */
double step_function(double dr)
{
    double t = dr - r_c;
    double t2 = t * t;
    double t4 = t2 * t2;

    return t4/(h + t4);
}


void force_calculation()
{
    double dx, dy, theta_ij;
    double dr, dr2;
    double t0, f0, NeighborSelect;
    double r_dist;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    for (int i = 0; i < N; omega_dot[i] = 0.0, i++);
    
    PE = 0.0;     /* PE initialize */
    // E_s = 0.5 * (p_s*p_s) / Q + 3.0 * N * Temperature() * log(s);

    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
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
                printf("Is this physical? %.12lf \n", dr2);
                exit(1);
            }
            
            if (dr2 < r_c*r_c)
            {
                dr = sqrt(dr2);

                /* ------------- For Omega_dot ------------- */
                theta_ij = theta[i] - theta[j];
                NeighborSelect = step_function(dr);
                t0 = - J * s * NeighborSelect * sin(theta_ij);
                omega_dot[i] = omega_dot[i] + t0;
                omega_dot[j] = omega_dot[j] - t0;
                /* ----------------------------------------- */

                /* -------------- For fx and fy -------------- */
                r_dist = pow((dr - r_c), 5);
                t0 = (NeighborSelect * NeighborSelect) / (dr * r_dist);
                t0 = t0 * (4 * J * s * h) * cos(theta_ij);

                fx[i] = fx[i] + t0 * dx;
                fy[i] = fy[i] + t0 * dy;

                fx[j] = fx[j] - t0 * dx;
                fy[j] = fy[j] - t0 * dy;
                /* ------------------------------------------- */


                PE -= J * NeighborSelect * cos(theta_ij);

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

                    PE += 4.0 * epsilon * sigma_by_r6 * (sigma_by_r6 - 1.0) + epsilon;
                }

            }

        }
        
    }

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

    p_s = p_s + (0.5 * dt) * M;
    
}


/* ----- Update p_s step: 2 ----- */
void update_ps_step2()
{
    double tt = 1.0 + (0.25 * dt / Q) * p_s;
    p_s = p_s / tt;
}


/* ----- Update s ----- */
void update_s()
{
    s = s * exp( (0.5*dt) * p_s / Q);
}


/* ----- Update theta ----- */
void update_theta()
{
    double s_inv = 1.0 / s;

    for (int i = 0; i < N; i++)
    {
        theta[i] += omega[i] * (0.5*dt) * s_inv / Intertia;

        // wrap angle b/w [0, 2pi)
        if (theta[i] >= 2 * pi) theta[i] -= 2.0 * pi;
        else if (theta[i] < 0.0) theta[i] += 2.0 * pi;
    }
    
}



void integrate()
{
    double xdot, ydot;
    double s_inv;

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
        xdot = px[i] * s_inv - K * cos(theta[i]);
        x[i] = x[i] + xdot * dt;
        if (x[i] >= L) x[i] = x[i] - L;
        else if (x[i] < 0) x[i] = x[i] + L;

        ydot = py[i] * s_inv - K * sin(theta[i]);
        y[i] = y[i] + ydot * dt;
        if (y[i] >= L) y[i] = y[i] - L;
        else if (y[i] < 0) y[i] = y[i] + L;
    }

    /* -------- end-kick --------- */
    update_theta();
    update_s();
    force_calculation();        // Recalculate forces 
    update_ps_step2();
    update_ps_step1();
    update_cano_lin_momenta();
    update_torque();            // Recalculate omega_dot
    update_cano_ang_momenta();
    /* ---------------------------- */

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



/* Rem. to free memory at end */
void deallocate_arrays()
{
    free(x); free(y); free(theta);

    free(px); free(py); free(omega);

    free(fx); free(fy); free(omega_dot);

}


int main()
{
    int swit = 1;   // switch to control whether final config reads/write
    double sumpx;

    allocate_arrays();
    printf("Arrays are allocated. \n");
    printf(" -------------------------- \n");
    
    if (swit)
    {
        load_initial_position();
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


    Temperature();
    printf("Temperature of system : %f \n", T);
    printf(" --------------------------------------- \n");

    force_calculation();
    printf("Initial Force calculation is done. \n");
    printf(" --------------------------------------- \n");
    
    update_torque();

    E_s = 0.5 * (p_s*p_s) / Q + 3.0 * N * T_res * log(s);
    H_nose0 = KE + RE + PE + E_s;

    FILE *fid1, *fid2, *fid3;
    fid1 = fopen("Energy_data_NVT_N128_rho1.20_K0.0_dt0.001", "w");
    fid2 = fopen("Temperature_data_NVT_N128_rho1.20_K0.0_dt0.001", "w");
    fid3 = fopen("Extended_var_NVT_N128_rho1.20_K0.0_dt0.001", "w");

    if (!fid1 || !fid2 || !fid3)
    {
        printf("Error! while opening files.\n");
        return 1;
    }
 
    for (int step = 0; step < max_step; step++)
    {
        if (step % 1 == 0)
        {
            save_data_infile(step, x, y, theta);
        }
        

        Temperature();
        E_s = 0.5 * (p_s*p_s) / Q + 3.0 * N * T_res * log(s);
        E = KE + RE + PE + E_s; 

        // Save thermodynamics property
        fprintf(fid1, "%d\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\n ", 
                        step, s*KE, s*RE, s*PE, s*E_s, s*(E-H_nose0) );
        fprintf(fid2, "%d\t %.12lf\t %.12lf\t %.12lf\t %.12lf\n", 
                        step, T1, T2, T3, T);
        fprintf(fid3, "%d\t %.12lf\t %.12lf\n", step, s, p_s);
        
        
        // Lets forward in time
        integrate();

        // debug
        sumpx = 0.0;
        for (int i = 0; i < N; i++)
        {
            sumpx = sumpx + py[i];
        }
        if (abs(sumpx) > 1e-8)
        {
            printf("Net motion in system: %lf \n", sumpx);
        }
        
        

    }
    
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);


    // Change to true to run
    if (swit)
    {
        FILE *fid4;
        fid4 = fopen("final_config", "w");
        if( fid4 == NULL  ) return(1);

        // fprintf(fid4, "%.12lf %.12lf \n", s, p_s);
        for( int i = 0; i < N; i++ )
        {
            fprintf(fid4, "%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf \n",
                x[i], y[i], theta[i], px[i], py[i], omega[i]);
        }
        fclose(fid4); 
    }

    deallocate_arrays();

    return 0;

}


