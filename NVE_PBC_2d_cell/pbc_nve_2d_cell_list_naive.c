/*-----------------------------------------------------*
 *        NVE simulation(PBC) with cell list           *
 *-----------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

/* -------- For LJ Potential -------- */
const double epsilon = 1.0;
const double sigma = 1.0;
const double r_c = 2.5;
const double sigma_by_rc = sigma/r_c;


/* -------- System Variable -------- */
const int N = 100;
const int max_step = 8000;
const double rho = 1.20;
const double vol = (double) N / rho;
const double dt = 0.001;
const double init_temp = 1.0;

/* -------- particle variables -------- */
double *x, *y;
double *vx, *vy;
double *fx, *fy;

/* -------- For Verlet list -------- */
#define max_neigh 50
const double skin_depth = 0.3 * sigma;
const double rl = (r_c + skin_depth);
int list_counter;
typedef struct
{
    int list[max_neigh];
    int max_counter;
}NEIB;

NEIB *neighbour;

/* -------- For Cell list -------- */
#define max_part_in_cell 50
const double cell_cut = (r_c + skin_depth);  // size of a cell we want
int n_cellx, n_celly, N_cell;                // numbers of cell
double cell_size;                            // true size of cell
int *cell_map[4];                            // 4 cells I choose for 2d
typedef struct                               // list of particles in a cell with max. numbers
{
    int list[max_part_in_cell];
    int max_counter;
}CLIST;    


/* -------- Req. Later -------- */
double sigma2, sigma_by_r2, sigma_by_r6;
double C0, C2;
double L;
double PE, KE;
const double PI = 3.141592653589;


/* ---------- Initialize global variables -------------- */
void global_variable_initializer()
{
    sigma2 = sigma * sigma;

    C0 = 4.0 * pow(sigma_by_rc,6.0) - 7.0 * pow(sigma_by_rc,12.0);
    C2 = 6.0 * pow(sigma_by_rc,14.0) - 3.0 * pow(sigma_by_rc,8.0);

    L = pow(vol, 0.50); 


    /* For cell list */
    n_cellx = n_celly = (int) floor(L / cell_cut);
    N_cell = n_cellx * n_celly;
    cell_size = L / (double)n_cellx;   /* cell_size > cell_cut */

}


/* ---------- Allocate arrays for storing -------------- */
void allocate_arrays()
{
    x = (double *)malloc(sizeof(double) * N);
    y = (double *)malloc(sizeof(double) * N);
    if(x == NULL || y == NULL) printf("Allocation request denied for positions. \n");


    vx = (double *)malloc(sizeof(double) * N);
    vy = (double *)malloc(sizeof(double) * N);
    if(vx == NULL || vy == NULL) printf("Allocation request denied for velocity. \n");


    fx = (double *)malloc(sizeof(double) * N);
    fy = (double *)malloc(sizeof(double) * N);
    if(fx == NULL || fy == NULL) printf("Allocation request denied for Forces. \n");


    neighbour = (NEIB *)malloc(sizeof(NEIB) * N);
    if(neighbour == NULL) printf("Allocation request denied for neighbour. \n");

    
    for (int i = 0; i < 4; i++)
    {
        cell_map[i] = (int *)malloc(sizeof(int) * N_cell );
        if (cell_map[i] == NULL) printf("Allocation request denied for cell_map \n");
    }
    
}



/* ---------- SC configarion -------------- */
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



/* ---------- Random Velocity -------------- */
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



/* ---------- Read a final config -------------- */
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



/* ----------Check temperature at any instace ------- */
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


/* ---- calculate distance between i and j particles ---- */
double particle_distance(int i, int j)
{
    if (j >= N || j < 0) {
        printf("Invalid index j=%d, debug_point=%d\n", j, 1);
        exit(1);
    }

    double dx, dy, dr2;
    
    dx = x[i] - x[j];
    dy = y[i] - y[j];

    // minimum image condition
    if(dx > 0.50 * L) dx = dx - L;
    if(dx < -0.50 * L) dx = dx + L;

    if(dy > 0.50 * L) dy = dy - L;
    if(dy < -0.50 * L) dy = dy + L;

    dr2 = dx*dx + dy*dy;

    return sqrt(dr2);
}


/* ------ give cell index for (i,j) pair box ------- */
int getcell_index(int i, int j)
{
    // apply pbc box
    if (i < 0) i = i + n_cellx;
    if (j < 0) j = j + n_celly;

    if (i >= n_cellx) i = i - n_cellx;
    if (j >= n_celly) j = j - n_celly;

    int cellno = i + j * n_cellx;
    return cellno;
}


/* ------ Give appropiate neighbours(L shape) cell num ---------- */
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



void construct_cell_list()
{
    int cell_index, L_index, k;         // neigh cell index L type 
    int nbc;                            // neighbour cell (1 to 4)
    int part_to_cell_index[N];          // give i particle in which cell

    double dr;


    CLIST cell_list[N_cell];


    // restart counter from 0
    list_counter = 0;

    // --> Initialize Neighbour list
    for (int i = 0; i < N; i++)
    {
        neighbour[i].max_counter = 0;
        for (int j = 0; j < max_neigh; j++)
        {
            neighbour[i].list[j] = 0;
        }
    }

    // --> Initialize cell list
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
        if (k >= N_cell || k < 0) {
            printf("Invalid index k=%d, debug_point=%d\n", k, 2);
            exit(1);
        }

        cell_list[k].list[cell_list[k].max_counter] = i;
        cell_list[k].max_counter ++;
        part_to_cell_index[i] = k;
    }


    /* ~~~~~ generate neighbour list using cell structure ~~~~~ */
    for (int i = 0; i < N; i++)
    {
        cell_index = part_to_cell_index[i];

            /******** Intra cell neighbour counting ********/
        for (int j = 0; j < cell_list[cell_index].max_counter; j++)
        {
            k = cell_list[cell_index].list[j];

            if (k >= N || k < 0) {
                printf("Invalid index k=%d, debug_point=%d\n", j, 3);
                exit(1);
            }

            if (k > i)
            {
                dr = particle_distance(i, k);
                if (dr < rl)
                {
                    neighbour[i].list[neighbour[i].max_counter] = k;
                    neighbour[k].list[neighbour[k].max_counter] = i;
                    neighbour[i].max_counter ++;
                    neighbour[k].max_counter ++;
                } 
            }   
        }

        /******** Inter cell neighbour counting ********/
        for (nbc = 0; nbc < 4; nbc++)
        {
            L_index = cell_map[nbc][cell_index];
            
            if (L_index >= N_cell || L_index < 0) {
                printf("Invalid index k=%d, debug_point=%d\n", L_index, 4);
                exit(1);
            }

            for (int j = 0; j < cell_list[L_index].max_counter; j++)
            {
                k = cell_list[L_index].list[j];

                dr = particle_distance(i, k);
                if (dr < rl)
                {
                    neighbour[i].list[neighbour[i].max_counter] = k;
                    neighbour[k].list[neighbour[k].max_counter] = i;
                    neighbour[i].max_counter ++;
                    neighbour[k].max_counter ++;
                } 
            }
        }

    }


    int temp;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < neighbour[i].max_counter-1; j++)
        {
            for (int k = j+1; k < neighbour[i].max_counter; k++)
            {
                if (neighbour[i].list[j] > neighbour[i].list[k])
                {
                    temp = neighbour[i].list[j];
                    neighbour[i].list[j] = neighbour[i].list[k];
                    neighbour[i].list[k] = temp;
                }
                
            }
            
        }
    }


    // Sorting the neighbour to match verlet
    /***********************************************************
    FILE *fid;
    fid = fopen("data_compare1.txt","w");

    int temp;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < neighbour[i].max_counter-1; j++)
        {
            for (int k = j+1; k < neighbour[i].max_counter; k++)
            {
                if (neighbour[i].list[j] > neighbour[i].list[k])
                {
                    temp = neighbour[i].list[j];
                    neighbour[i].list[j] = neighbour[i].list[k];
                    neighbour[i].list[k] = temp;
                }
                
            }
            
        }

        for (int j = 0; j < neighbour[i].max_counter; j++)
        {
            fprintf(fid, "%d \t", neighbour[i].list[j]);
        }
        fprintf(fid, "\n");
        
    }
    fclose(fid);
    **************************************************************/

}



void construct_verlet_list()
{
    double dr;

    // restart counter from 0
    list_counter = 0; 

    // Initialize Neighbour list
    for (int i = 0; i < N; i++)
    {
        neighbour[i].max_counter = 0;
        for (int j = 0; j < max_neigh; j++)
        {
            neighbour[i].list[j] = 0;
        }
        
    }

    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            dr = particle_distance(i, j);
            if (dr < rl)
            {
                neighbour[i].list[neighbour[i].max_counter] = j;
                neighbour[i].max_counter ++;
            }
        }
    }  
}



void force_calculation()
{
    double dx, dy, dr2;
    double f0;

    int j, tot_neigh;

    for (int i = 0; i < N; fx[i] = fy[i] = 0.0, i++);
    
    PE = 0.0;

    for (int i = 0; i < N; i++)
    {
        tot_neigh = neighbour[i].max_counter;
        for (int k = 0; k < tot_neigh; k++)
        {
            j = neighbour[i].list[k]; // j: neighbour particle index
            if (i < j)
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
                    printf("Particles came so close. \n");
                    exit(0);
                }

                if ( sqrt(dr2) < r_c )
                {
                    sigma_by_r2 = sigma2/dr2;
                    sigma_by_r6 = sigma_by_r2 * sigma_by_r2 * sigma_by_r2;

                    f0 = 8.0 * epsilon * ( (6.0 * sigma_by_r6) *
                                        (sigma_by_r6 - 0.50)/dr2 - C2/(sigma2) );
                    
                    fx[i] = fx[i] + f0 * dx;
                    fx[j] = fx[j] - f0 * dx;

                    fy[i] = fy[i] + f0 * dy;
                    fy[j] = fy[j] - f0 * dy;

                    PE = PE + 4.0 * epsilon * ( sigma_by_r6 * (sigma_by_r6 - 1.0) 
                                                    + C0 + C2 / sigma_by_r2 );
                }

            }
        }   
    }

    PE = PE / (double) N;   //PE per particle
} 




void integrate()
{
    double dx, dy, dr, max_disp;

    max_disp = 0.0;

    for (int i = 0; i < N; i++)
    {
        // velocity-verlet: step1
        vx[i] = vx[i] + fx[i] * (dt/2.0);
        vy[i] = vy[i] + fy[i] * (dt/2.0);


        dx = vx[i] * dt;
        dy = vy[i] * dt;
        dr = sqrt(dx*dx + dy*dy);
        if (dr > max_disp) max_disp = dr;

        x[i] = x[i] + dx;
        y[i] = y[i] + dy;
        
        // applying PBC
        if ( x[i] >= L )  x[i] = x[i] - L;
        else if ( x[i] < 0.0 ) x[i] = x[i] + L;

        if ( y[i] >= L )  y[i] = y[i] - L;
        else if ( y[i] < 0.0 ) y[i] = y[i] + L;
    
    }

    list_counter ++;
    if ( (2.0 * max_disp * list_counter) > skin_depth )
    {
        // printf("clist constructed after step: %d \n", list_counter);
        // construct_verlet_list();
        construct_cell_list();
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

    free(neighbour);

    for (int i = 0; i < 4; i++)
    {
        free(cell_map[i]);
    }
    
}



int main()
{
    
    double temp;
    FILE *fid1, *fid2;

    // Assign global variables
    global_variable_initializer();

    printf("L = %f \n", L);
    printf("cell size = %f \n", cell_size);

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
    printf("-------------------------------- \n");

    maps();
    // construct_verlet_list();
    construct_cell_list();
    printf("Cell list is ready. \n");
    printf("-------------------------------- \n");

    force_calculation();
    printf("Initial force given to system. \n");
    printf("------------------------------------ \n");


    // printf("Initial force given to system. \n");
    // printf("------------------------------------ \n");

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

    // Change to TRUE to save final_config
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


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * The function construct_verlet_list() is newly updated which reduces    *
 * the cost. Force & Integration function also need to update slightly    *
 * The double loop updated in force subrouitine. Integration module       *
 * contains when the neighbour_list need to be updated automatically      *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* -------------------------------------------------------------------------- 
 Benchmark for this program is set by plotting energy vs time from costly
 force calculation and using verlet list. Both programming gives exactly same 
 result. Please remember to read same configuration for both cases. 
 --------------------------------------------------------------------------- */

 


 /* 

  +---+---+---+
  | 6 | 7 | 8 |
  +---+---+---+
  | 3 | 4 | 5 |         >> CELL STRUCTURE 
  +---+---+---+
  | 0 | 1 | 2 |
  +---+---+---+

*/

/*  
 *******************************************************  
 * Output for MAP array for N_cell = 9                 *  
 * It matched if we do it by hand:                     *  
 *                                                     *  
 *   1   2   0   4   5   3   7   8   6                 *  
 *   4   5   3   7   8   6   1   2   0                 *  
 *   3   4   5   6   7   8   0   1   2                 *  
 *   5   3   4   8   6   7   2   0   1                 *  
 *******************************************************  
 */
