/*-----------------------------------------------------*
 *           NVE simulation with PBC                   * 
 *-----------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

// For LJ
const double epsilon = 1.0;
const double sigma = 1.0;
const double r_c = 2.5;
const double sigma_by_rc = sigma/r_c;


// System variables
const int N = 100;
const int max_step = 1000;
const double rho = 1.20;
const double vol = (double) N / rho;
const double dt = 0.001;
const double init_temp = 1.0;

// Particles variables
double *x, *y;
double *vx, *vy;
double *fx, *fy;

// Req. later
double sigma2, sigma6, sigma12;
double C0, C2;
double L;

double PE, KE;

const double PI = 3.141592653589;

// Global Variable initilizer
void global_variable_initializer()
{
    sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;
    sigma6 = sigma4 * sigma2;
    sigma12 = sigma6 * sigma6;

    C0 = 4.0 * pow(sigma_by_rc,6.0) - 7.0 * pow(sigma_by_rc,12.0);
    C2 = 6.0 * pow(sigma_by_rc,14.0) - 3.0 * pow(sigma_by_rc,8.0);

    L = pow(vol, 0.50); 
}


// Allocating Arrays for Use
void allocate_arrays()
{
    x = (double *)malloc(sizeof(double) * N);
    y = (double *)malloc(sizeof(double) * N);

    if(x == NULL || y == NULL)
    {
        printf("Allocation request denied for positions. \n");
    }


    vx = (double *)malloc(sizeof(double) * N);
    vy = (double *)malloc(sizeof(double) * N);
    if(vx == NULL || vy == NULL)
    {
        printf("Allocation request denied for velocity. \n");
    }


    fx = (double *)malloc(sizeof(double) * N);
    fy = (double *)malloc(sizeof(double) * N);
    if(fx == NULL || fy == NULL)
    {
        printf("Allocation request denied for Forces. \n");
    }
}


void load_initial_position()
{ 
    // No of particle in x, y, z direction
    int n = (int) round(pow( (double) N, 0.50));
        
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


void load_initial_velocity()
{
    double sumvx = 0.0, sumvy = 0.0;
    double sumv2 = 0.0, sfactor;

    // Initize the seed 
    srand48(time(NULL));

    for(int i=0; i<N; i++)
    {
        vx[i] = 2.0 * drand48() - 1.0;
        vy[i] = 2.0 * drand48() - 1.0;

        sumvx = sumvx + vx[i];
        sumvy = sumvy + vy[i];
    }


    for(int i=0; i<N; i++)
    {
        vx[i] = vx[i] - sumvx / (double) N;
        vy[i] = vy[i] - sumvy / (double) N;

        sumv2 = sumv2 + vx[i]*vx[i] + vy[i]*vy[i];
    }

    sfactor = sqrt(2.0* N * init_temp / sumv2);

    for(int i = 0; i < N; i++)
    {
        vx[i] = vx[i] * sfactor;
        vy[i] = vy[i] * sfactor;
    }

}


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
        fscanf(fid_load, "%lf %lf %lf %lf",
                &x[i], &y[i], &vx[i], &vy[i]);
    }

    fclose(fid_load);
}


double Temperature()
{
    double T, sumv2 = 0;

    for(int i=0; i<N; i++)
    {
        sumv2 = sumv2 + vx[i]*vx[i] + vy[i]*vy[i];
    }
    
    T = sumv2 / (2.0 * N);

    return T;
}


void force_calculation()
{
    double dx, dy;
    double dr2, dr4, dr8;
    double f0;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    
    PE = 0.0;

    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            dx = x[i] - x[j];
            dy = y[i] - y[j];

            // minimum image condition
            if(dx > 0.50 * L) dx = dx - L;
            if(dx < -0.50 * L) dx = dx + L;

            if(dy > 0.50 * L) dy = dy - L;
            if(dy < -0.50 * L) dy = dy + L;

            dr2 = dx*dx + dy*dy;

            if ( dr2 < sigma2 / 100.0 ) 
            {
                printf("Particles came so close.");
                exit(0);
            }

            if ( sqrt(dr2) < r_c )
            {
                dr4 = dr2 * dr2;
                dr8 = dr4 * dr4;

                f0 = 8.0 * epsilon * ( (6.0 * sigma12)/(dr8 * dr4 * dr2) - 
                        (3.0 * sigma6)/dr8 - C2/(sigma2) );
                
                fx[i] = fx[i] + f0 * dx;
                fx[j] = fx[j] - f0 * dx;

                fy[i] = fy[i] + f0 * dy;
                fy[j] = fy[j] - f0 * dy;

                PE = PE + 4.0 * epsilon * ( (sigma12/(dr8 * dr4)) -
                 (sigma6 * dr2 / dr8) + C0 + C2 * dr2 / sigma2 );

            }
         }   
     }

    PE = PE / (double) N;   //PE per particle

} 



void integrate()
{
    for (int i = 0; i < N; i++)
    {
        // velocity-verlet: step1
        vx[i] = vx[i] + fx[i] * (dt/2.0);
        vy[i] = vy[i] + fy[i] * (dt/2.0);


        x[i] = x[i] + vx[i] * dt;
        y[i] = y[i] + vy[i] * dt;
        
        // applying PBC
        if ( x[i] >= L )  x[i] = x[i] - L;
        if ( x[i] < 0.0 ) x[i] = x[i] + L;

        if ( y[i] >= L )  y[i] = y[i] - L;
        if ( y[i] < 0.0 ) y[i] = y[i] + L;
    }
    
    // Calculate new force
    force_calculation();

    for( int i = 0; i < N; i++)
    {
       vx[i] = vx[i] + fx[i] * (dt/2.0);
       vy[i] = vy[i] + fy[i] * (dt/2.0);
    }

}


void deallocate_arrays()
{
    free(x);
    free(y);

    free(vx);
    free(vy);

    free(fx);
    free(fy);
}


int main()
{
    
    double temp;
    FILE *fid1, *fid2;

    // Assign global variables
    global_variable_initializer();


    allocate_arrays();
    printf("Arrays are allocated. \n");
    printf("----------------------- \n");

    // load_initial_position();
    // printf("Initial Position is ready. \n");
    // printf("--------------------------- \n");

    // load_initial_velocity();
    // printf("Initial velocity is given. \n");
    // printf("-------------------------------- \n");

    load_final_config();

    temp = Temperature();
    printf("Initial temperature of the system: %8.4f \n", temp);

    force_calculation();
    printf("Initial force given to system. \n");
    printf("------------------------------------ \n");

    fid1 = fopen("animation1.xyz", "w");
    if( fid1 == NULL ) return(1);

    fid2 = fopen("T_PE_KE_E_data.txt", "w");
    if( fid2 == NULL ) return(1);


    printf("simulation started ... \n");
    

    for(int i = 0; i < max_step; i++)
    {
        // writing positions
        fprintf(fid1, "%d \n", N);
        fprintf(fid1, "%s \n", " ");

        for( int j = 0; j < N; j++ )
        {
            fprintf(fid1, "%s \t %f \t %f \t %f \n",
                    "Ar", x[j], y[j], 0.0);
        }

        // writing T, PE, KE, E
        temp = Temperature();
        
        KE = temp;    // KE = N*KB*T
        fprintf(fid2, "%d \t %f \t %f \t %f \t %f \n",
                i, temp, KE, PE, (PE+KE));

        integrate();
    }

    // Change to true to run
    if (0)
    {
        FILE *fid3;
        fid3 = fopen("final_config.txt", "w");
        if( fid3 == NULL  ) return(1);
        for( int i = 0; i < N; i++ )
        {
            fprintf(fid3, "%f\t%f\t%f\t%f\n", 
                    x[i], y[i], vx[i], vy[i]);
        }
        fclose(fid3); 
    }
    

    deallocate_arrays();

    fclose(fid1);
    fclose(fid2);

    return 0;

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! In this program we implement a md simulation for 2d system !
! with PBC. Three importent benchmark has done :             !
! 1. Total energy is conserved upto 3 decimal points         !
! 2. Energy error flactuates as dt^2                         !
! 3. Force function checked extensively by a benchmark set   !
!    by Antik. Force Disbalenced F=sqrt(Fx^2 + Fy^2) is      !
!    matched with his code result. (0.0000000036204378)      !
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
