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
#define     N           2                            /* Number of particles ---------- */
#define     rho         0.50                           /* density of system ------------ */
#define     L           sqrt((double) N / rho)         /* Lenght of box ---------------- */

const double dt = 0.001;                               /* time step for integration ---- */
const int    max_step = 1;                          /* No. of steps simulation run -- */
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
double p_s = 0.0, s = 1.0;                         /* extended dimensions ---*/
double KE, PE, RE, E;                                  /* Energies ------------- */
double H_nose0;                                    /* H_Nose at t=0 -------- */
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
}


/* --- Load initial position --- */
void load_initial_position()
{ 
    // No of particle in x, y direction
    int n = (int) ceil( sqrt((double) N) );
        
    // Distances b/w particles
    double a = L / (double) n;

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
            x[pt_num] = i * a + a / 2.0 ;
            y[pt_num] = j * a + a / 2.0 ;

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


    for(int i=0; i<N; i++)
    {
        px[i] = px[i] - sumpx / (double) N;
        py[i] = py[i] - sumpy / (double) N;
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

    fid_load = fopen("final_config.txt", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config.txt");
        exit(1);
    }

    for ( int i = 0; i < N; i++ )
    {
        fscanf(fid_load, "%lf %lf %lf %lf %lf %lf",
                &x[i], &y[i], &theta[i], &px[i], &py[i], &omega[i]);
    }

    fclose(fid_load);
}


/* ---- Calculate KE of the system ---- */
/* --- Calculate temperature of system at current state --- */
double Temperature()
{
    double sump2 = 0.0, sum_omega2 = 0.0;
    double sum_pA = 0.0;
    double temp;

    for (int i = 0; i < N; i++)
    {
        sump2      += px[i]*px[i] + py[i]*py[i];
        sum_omega2 += omega[i]*omega[i];
        sum_pA     += px[i]*cos(theta[i]) + py[i]*sin(theta[i]);
    }
    
    // T = \frac{1}{3NK_B} [ p_i^2/m - p.A/m + w^2/I]
    temp = sump2 + sum_omega2 / Intertia - K * sum_pA; 
    temp = temp / (3.0 * N);

    // KE = sum [ p_i^2 / 2m - p.A/m + w^2/2I] + NK^2/2m
    KE = 0.50 * sump2 - K * sum_pA + 0.50 * (N * K * K);
    KE = KE / (double) N;

    RE = 0.50 * sum_omega2 / ( N * Intertia);       // Rotational Energy per particles

    return temp;
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
    double t0, f0, int_strength;      /* Interaction strength */
    double r_dist;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    for (int i = 0; i < N; omega_dot[i] = 0.0, i++);
    
    PE = 0.0;                         /* PE initialize */


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
                printf("Is this physical? %f \n", dr2);
                exit(1);
            }
            
            if (dr2 < r_c*r_c)
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

                if (dr < r_wca )
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

    // check omega_dot
    double s_omega_dot = 0.0;
    for (int i = 0; i < N; i++)
    {
        s_omega_dot += omega_dot[i];
    }
    printf("sum over omega = %lf \n", s_omega_dot);
    
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
        omega[i] = omega[i] + 0.50 * omega_dot[i] * dt;
    }
    
}


/*--- Canonical Momentum Updatation (px, py) ---*/
void update_cano_lin_momenta()
{
    for (int i = 0; i < N; i++)
    {
        px[i] = px[i] +  0.50 * fx[i] *dt;
        py[i] = py[i] +  0.50 * fy[i] * dt;
    }

}


/*----- Update theta -----*/
void update_theta()
{
    for (int i = 0; i < N; i++)
    {
        theta[i] = theta[i] + 0.50 * omega[i] * dt / Intertia;

        // wrap angle b/w [0, 2pi)
        if (theta[i] >= 2 * pi) theta[i] -= 2.0 * pi;
        else if (theta[i] < 0.0) theta[i] += 2.0 * pi;
    }
    
}

/* ------ Update positions ------- */
void update_position()
{
    double x_dot, y_dot;
    for (int i = 0; i < N; i++)
    {
        x_dot = px[i] - K * cos(theta[i]);
        x[i] = x[i] + x_dot * dt;
        if (x[i] >= L) x[i] = x[i] - L;
        else if (x[i] < 0) x[i] = x[i] + L;

        y_dot = py[i] - K * sin(theta[i]);
        y[i] = y[i] + y_dot * dt;
        if (y[i] >= L) y[i] = y[i] - L;
        else if (y[i] < 0) y[i] = y[i] + L;
    }
    
}


void integrate()
{
    /* --------- half-kick ---------- */
    update_cano_ang_momenta();
    update_cano_lin_momenta();
    update_theta();
    /* ------------------------------ */

    update_position();

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
    sprintf(filename, "movie_data/step_%04d.txt", step); // "step_0001.txt", "step_0002.txt", ...

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        fprintf(file, "%lf %lf %lf\n", x[i], y[i], theta[i]);
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
    double temp;

    allocate_arrays();
    printf("Arrays are allocated. \n");

    load_initial_position();
    printf("Initial position is loaded. \n");

    load_initial_velocity();
    printf("Initial Momentum is given. \n");

    load_initial_spin();
    printf("Initial Spin and Ang. Momentum is given. \n");

    printf("Temperature of system : %f \n", Temperature());

    force_calculation();
    printf("Initial Force calculation is done. \n");
    update_torque();

    for (int i = 0; i < N; i++)
    {
        printf("fx = %.14lf, fy = %.14lf, torque = %.14lf \n", fx[i], fy[i], omega_dot[i]);
    }
    
    // exit(0);

    FILE *fid;
    fid = fopen("T_KE_PE_E_data", "w");
    if ( fid == NULL ) return 1;
 
    for (int step = 0; step < max_step; step++)
    {
        save_data_infile(step, x, y, theta);

        temp = Temperature();
        // Save thermodynamics property
        fprintf(fid, "%d\t %.10lf\t %.10lf\t %.10lf\t %.10lf\t %.10lf\n ", 
                    step, temp, KE, RE, PE, (KE+RE+PE) );
        
        integrate();

    }
    
    fclose(fid);



    // Change to true to run
    if (0)
    {
        FILE *fid3;
        fid3 = fopen("final_config.txt", "w");
        if( fid3 == NULL  ) return(1);
        for( int i = 0; i < N; i++ )
        {
            fprintf(fid3, "%lf %lf %lf %lf %lf %lf \n",
                x[i], y[i], theta[i], px[i], py[i], omega[i]);
        }
        fclose(fid3); 
    }

    deallocate_arrays();

    return 0;
}










/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * This is NVE simulation of the Actice Hamiltonian system presented  *
 * in the paper. As it is NVE system we don't need s, p_s variables.  *
 * We only try to plot energy of real hamiltonian which should be     *
 * conserved in this case. Along with flactuation dt^2.               *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/