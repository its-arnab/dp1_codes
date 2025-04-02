#include "global.h"
#include "subroutine.h"


/* ---------- Initialize global variables -------------- */
void global_variable_initializer()
{
    // general variables 
    L = sqrt(vol);
    sigma2 = sigma * sigma;

    C0 = 4.0 * pow(sigma_by_rc, 6.0) - 7.0 * pow(sigma_by_rc, 12.0);
    C2 = 6.0 * pow(sigma_by_rc, 14.0) - 3.0 * pow(sigma_by_rc, 8.0);

    // For cell list 
    cell_cut = rl;
    n_cellx = n_celly = (int) floor(L / cell_cut);
    N_cell = n_cellx * n_celly;
    cell_size = L / (double)n_cellx;   /* cell_size > cell_cut */

    if (N_cell < 9)
    printf("Cell List is NOT possible. As no of cell = %d\n", N_cell);

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Allocate arrays for storing Positions,           *
 * Velocites, Forces, Neighbours of particles       *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void allocate_arrays()
{
    x = malloc(sizeof(double) * N);
    y = malloc(sizeof(double) * N);
    if(x == NULL || y == NULL) printf("Allocation request denied for positions. \n");


    vx = malloc(sizeof(double) * N);
    vy = malloc(sizeof(double) * N);
    if(vx == NULL || vy == NULL) printf("Allocation request denied for velocity. \n");


    fx = malloc(sizeof(double) * N);
    fy = malloc(sizeof(double) * N);
    if(fx == NULL || fy == NULL) printf("Allocation request denied for Forces. \n");


    neighbour = (NEIB *)malloc(sizeof(NEIB) * N);
    if(neighbour == NULL) printf("Allocation request denied for neighbour. \n");

    
    for (int i = 0; i < 4; i++)
    {
        cell_map[i] = (int *)malloc(sizeof(int) * N_cell );
        if (cell_map[i] == NULL) printf("Allocation request denied for cell_map \n");
    }
    
}


/* ------- Setup SC configarion -------- */
void load_initial_position()
{ 
    // No of particle in x, y direction
    int n = ceil( sqrt((double) N) );
        
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
    printf("Warning ! All particles are NOT placed. \n");

}


/* Gives random no.s between 0 and 1 (both inclusive) */
double rand_01()                                   
{
    return (double)rand() / (double)RAND_MAX ;
}


/* ------- Uniform Random Velocity -------- */
void load_initial_velocity()
{
    double sumvx = 0.0, sumvy = 0.0;
    double sumv2 = 0.0, sfactor;

    // Initize the seed 
    srand(time(NULL));

    for(int i=0; i<N; i++)
    {
        vx[i] = 2.0 * rand_01() - 1.0;
        vy[i] = 2.0 * rand_01() - 1.0;

        sumvx = sumvx + vx[i];
        sumvy = sumvy + vy[i];
    }

    sumvx = sumvx / (double) N;
    sumvy = sumvy / (double) N;

    for(int i=0; i<N; i++)
    {
        vx[i] = vx[i] - sumvx;
        vy[i] = vy[i] - sumvy;

        sumv2 += vx[i]*vx[i] + vy[i]*vy[i];
    }

    // scale factor to rescale velocity to get desire temperature
    sfactor = sqrt(2.0* N * init_temp / sumv2);

    for(int i = 0; i < N; i++)
    {
        vx[i] = vx[i] * sfactor;
        vy[i] = vy[i] * sfactor;
    }

}


/* ------- Read a final config from file --------- */
void load_final_config()
{
    FILE *fid_load;

    fid_load = fopen("data_files/final_config", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config.txt");
        exit(1);
    }

    for ( int i = 0; i < N; i++ )
    {
        if (fscanf(fid_load, "%lf %lf %lf %lf", &x[i], &y[i], &vx[i], &vy[i]) != 4) {
            fprintf(stderr, "Error reading data from file at line %d\n", i);
            exit(EXIT_FAILURE);
        }        
    }

    fclose(fid_load);
}


/* ------ Program Initilizer(NCM) ------ */
void Initilize_Program(int mode)
{
    global_variable_initializer();
    allocate_arrays();

    printf(" Variables & Arrays are ready. \n");
    printf("-------------------------------\n");

    if (mode==1)
    {
        load_initial_position();
        load_initial_velocity();

        printf(" Initial Position & Velocity given properly. \n");
        printf("---------------------------------------------\n");
    }
    else
    {
        load_final_config();
        printf(" Previous configuration reads successfully. \n");
        printf("---------------------------------------------\n");
    }
    
}



/* ------ free heap memory(NCM) ------ */
void deallocate_arrays()
{
    free(x); free(y);

    free(vx); free(vy);

    free(fx); free(fy);

    free(neighbour);

    for (int i = 0; i < 4; i++) free(cell_map[i]);
    
}


// NCM = Need to Call in Main();
