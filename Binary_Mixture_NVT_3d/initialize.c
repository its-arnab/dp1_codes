#include "global.h"
#include "subroutine.h"


/* ---------- Initialize global variables -------------- */
void global_variable_initializer()
{
    // general variables 
    L = pow(vol, 1.0/3.0);
    half_L = 0.50 * L;

    sigma[0][0] = 1.00;         /*~~~ sigma_AA ~~~*/ 
    sigma[0][1] = 0.80;         /*~~~ sigma_AB ~~~*/ 
    sigma[1][0] = 0.80;         /*~~~ sigma_BA ~~~*/ 
    sigma[1][1] = 0.88;         /*~~~ sigma_BB ~~~*/ 
    
    epsilon[0][0] = 1.00;       /*~~~ epsilon_AA ~~~*/
    epsilon[0][1] = 1.50;       /*~~~ epsilon_AB ~~~*/
    epsilon[1][0] = 1.50;       /*~~~ epsilon_BA ~~~*/
    epsilon[1][1] = 0.50;       /*~~~ epsilon_BB ~~~*/

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

    xu = malloc(sizeof(double) * N);
    yu = malloc(sizeof(double) * N);
    zu = malloc(sizeof(double) * N);
    if(xu == NULL || yu == NULL || zu == NULL)
    printf("Allocation request denied for unwrapped positions. \n");

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

    pt_kind = malloc(sizeof(int) * N);
    if (pt_kind == NULL)
    printf("Allocation request denied for pt_kind. \n");

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
    int n = 1;

    while (n*n*n < N) n++;
        
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
                x[pt_num] = i * a + 0.50 * a;
                y[pt_num] = j * a + 0.50 * a;
                z[pt_num] = k * a + 0.50 * a;
                pt_num++ ;
            }
            if (exit_flag) break;
        }
        if(exit_flag) break;
    }

    if(pt_num < N)
    printf("Warning ! All particles are NOT placed. \n");

    // For Unwrap position
    for (int i = 0; i < N; i++)
    {
        xu[i] = x[i];
        yu[i] = y[i];
        zu[i] = z[i];
    }
}


void initialise_kind_random()
{
    // starting with all type A particles
    for (int i = 0; i < N; i++) pt_kind[i] = 0; 

    int j = 0;
    do 
    {
        int i = (int)(drand48() * N);        // random integer from 0 to N-1 (both inclusive)
        if (pt_kind[i] != 1)
        {
            pt_kind[i] = 1;        //  assigning type B particles randomly
            j++;
        }

    } while (j < (N-NA) );

}


void load_initial_velocity()
{
    double sumvx = 0.0, sumvy = 0.0, sumvz = 0.0;
    double sumv2 = 0.0, sfactor;


    for(int i=0; i<N; i++)
    {
        vx[i] = 2.0 * drand48() - 1.0;
        vy[i] = 2.0 * drand48() - 1.0;
        vz[i] = 2.0 * drand48() - 1.0;

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

    fid_load = fopen("data_files/T_0.50000/run_60/final_config", "r");
    if (fid_load == NULL) {
        printf("Error: Unable to open file: final_config");
        exit(1);
    }

    // Read the Nose-Hoover variables
    for (int i = 0; i < M_chain; i++)
    {
        if (fscanf(fid_load, "%lf %lf %lf", &eta[i], &p_eta[i], &Q[i]) != 3)
        {
            fprintf(stderr, "Error reading data from file at line %d\n", i);
            exit(EXIT_FAILURE);
        }
    }

    // Read the particle data
    for ( int i = 0; i < N; i++ )
    {
        if (fscanf(fid_load, "%lf %lf %lf %lf %lf %lf %d",
             &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i], &pt_kind[i]) != 7)
        {
            fprintf(stderr, "Error reading data from file at line %d\n", i);
            exit(EXIT_FAILURE);
        }        
    }

    fclose(fid_load);
}


/* ------ Program Initilizer(NCM) ------ */
void Initialize_Program(int mode)
{
    global_variable_initializer();
    allocate_arrays();

    printf(" Variables & Arrays are ready. \n");
    printf("-------------------------------\n");


    // Initize the seed as we use drand48() through out program
    long seed = time(NULL) ^ (getpid() << 16);
    srand48(seed); // seed for drand48()

    if (mode==1)
    {
        load_initial_position();
        initialise_kind_random();
        load_initial_velocity();

        printf(" Initial Position & Velocity given properly. \n");
        printf("---------------------------------------------\n");

        Initilize_NoseHoover();
        printf(" Nose-Hoover variables are initialized. \n");
        printf("----------------------------------------\n");
    }
    else
    {
        load_final_config();
        printf(" Previous configuration reads successfully. \n");
        printf("---------------------------------------------\n");
    }

    // Generate the time step for saving data
    md_step_generator();
    printf(" Time step for saving data is generated. \n");
    printf("----------------------------------------\n");
    
}



/* ------ free heap memory(NCM) ------ */
void deallocate_arrays()
{
    free(x); free(y); free(z);

    free(xu); free(yu); free(zu);

    free(vx); free(vy); free(vz);

    free(fx); free(fy); free(fz);

    free(pt_kind); free(neighbour);

    for (int i = 0; i < 13; i++) free(cell_map[i]);

    free(save_time);
    
}






// NCM = Need to Call in Main();
