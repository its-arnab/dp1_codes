#include "global.h"
#include "subroutine.h"


/* ---------- Initialize global variables -------------- */
void global_variable_initializer()
{
    // general variables 
    L = pow(vol, 1.0/3.0);
    sigma2 = sigma * sigma;

    C0 = 4.0 * pow(sigma_by_rc, 6.0) - 7.0 * pow(sigma_by_rc, 12.0);
    C2 = 6.0 * pow(sigma_by_rc, 14.0) - 3.0 * pow(sigma_by_rc, 8.0);

    // For cell list 
    cell_cut = rl;
    n_cellx = (int) floor(L / cell_cut);
    n_celly = n_cellz = n_cellx;
    N_cell = n_cellx * n_celly * n_cellz;
    cell_size = L / (double)n_cellx;   /* cell_size > cell_cut */

    if (N_cell < 27)
    printf("Cell List is NOT Possible as no of cell = %d \n", N_cell);

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Allocate arrays for storing Positions,           *
 * Velocites, Forces, Neighbours of particles       *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void allocate_arrays()
{
    x = malloc(sizeof(double) * N);
    y = malloc(sizeof(double) * N);
    z = malloc(sizeof(double) * N);
    if(x == NULL || y == NULL || z == NULL)
    printf("Allocation request denied for positions. \n");


    vx = malloc(sizeof(double) * N);
    vy = malloc(sizeof(double) * N);
    vz = malloc(sizeof(double) * N);
    if(vx == NULL || vy == NULL || vz == NULL)
    printf("Allocation request denied for velocity. \n");


    fx = malloc(sizeof(double) * N);
    fy = malloc(sizeof(double) * N);
    fz = malloc(sizeof(double) * N);
    if(fx == NULL || fy == NULL || fz == NULL)
    printf("Allocation request denied for Forces. \n");


    neighbour = (NEIB *)malloc(sizeof(NEIB) * N);
    if(neighbour == NULL) printf("Allocation request denied for neighbour. \n");

    
    for (int i = 0; i < 13; i++)
    {
        cell_map[i] = (int *)malloc(sizeof(int) * N_cell );
        if (cell_map[i] == NULL) printf("Allocation request denied for cell_map \n");
    }
    
}


/* ------- Setup SC configarion -------- */
void load_initial_position()
{ 
    // No of particle in x, y direction
    int n = ceil( pow((double) N, 1.0/3.0) );
        
    // Distances b/w particles
    double a = L / (double) n;

    int pt_num = 0;

    int exit_flag = 0;

    for (int i = 0; i < n; i++)             // x-loop
    {
        for (int j = 0; j < n; j++)         // y-loop
        {
            for (int k = 0; k < n; k++)
            {
                if(pt_num > N-1)
                {
                    exit_flag = 1;
                    break;
                }
                x[pt_num] = i * a + a / 2.0;
                y[pt_num] = j * a + a / 2.0;
                z[pt_num] = k * a + a / 2.0;
                pt_num++ ;
            }
            if (exit_flag) break;
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


void load_initial_velocity()
{
    double sumvx = 0.0, sumvy = 0.0, sumvz = 0.0;
    double sumv2 = 0.0, sfactor;

    // Initize the seed 
    srand(time(NULL));

    for(int i=0; i<N; i++)
    {
        vx[i] = 2.0 * rand_01() - 1.0;
        vy[i] = 2.0 * rand_01() - 1.0;
        vz[i] = 2.0 * rand_01() - 1.0;

        sumvx = sumvx + vx[i];
        sumvy = sumvy + vy[i];
        sumvz = sumvz + vz[i];
    }


    for(int i=0; i<N; i++)
    {
        vx[i] = vx[i] - sumvx / (double) N;
        vy[i] = vy[i] - sumvy / (double) N;
        vz[i] = vz[i] - sumvz / (double) N;

        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    sfactor = sqrt(3.0* N * init_temp / sumv2);

    for(int i = 0; i < N; i++)
    {
        vx[i] = vx[i] * sfactor;
        vy[i] = vy[i] * sfactor;
        vz[i] = vz[i] * sfactor;
    }

}



/* ------- Read a final config from file --------- */
void load_final_config()
{
    FILE *fid_load;

    fid_load = fopen("data_files/final_config", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config");
        exit(1);
    }

    for ( int i = 0; i < N; i++ )
    {
        if (fscanf(fid_load, "%lf %lf %lf %lf %lf %lf",
             &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]) != 6)
        {
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
    free(x); free(y); free(z);

    free(vx); free(vy); free(vz);

    free(fx); free(fy); free(fz);

    free(neighbour);

    for (int i = 0; i < 13; i++) free(cell_map[i]);
    
}


// NCM = Need to Call in Main();
