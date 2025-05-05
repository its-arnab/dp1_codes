#include "global.h"
#include "subroutine.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * calculate distance between i and j particles   *
 * Input: Particle index i, j                     *
 * Output: (dr/sig)^2 b/w the particles.          *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double particle_distance(int i, int j)
{
    double dx, dy, dz, dr2;
    
    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];

    // minimum image condition
    if(dx > half_L) dx = dx - L;
    else if(dx < - half_L) dx = dx + L;

    if(dy > half_L) dy = dy - L;
    else if(dy < - half_L) dy = dy + L;

    if(dz > half_L) dz = dz - L;
    else if(dz < - half_L) dz = dz + L;

    dr2 = dx*dx + dy*dy + dz*dz;

    return dr2;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Input: xI, yI, zI (cell index in x, y, z direction)          *
 * Output: cellID (cell ID as a number)                         *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int getcell_index(int xI, int yI, int zI)
{

    if (xI < 0) xI = xI + n_cellx;
    if (yI < 0) yI = yI + n_celly;
    if (zI < 0) zI = zI + n_cellz;

    if (xI >= n_cellx) xI = xI - n_cellx;
    if (yI >= n_celly) yI = yI - n_celly;
    if (zI >= n_cellz) zI = zI - n_cellz;
  
    int cellID =  xI + yI * n_cellx + zI * n_cellx * n_celly;
  
    return cellID;

}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Give appropiate neighbours(L shape) cell num for cell list *
 * For a sample output for n_cellx = 3, see comment at end of *
 * the program.                                               * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void maps()
{
    int imap;
    int xI, yI, zI;

    int shift[13][3] = {
        {-1, -1,  1}, {-1,  0,  1}, {-1,  1,  1},       /* Upper z plane */       
        { 0, -1,  1}, { 0,  0,  1}, { 0,  1,  1},       /* Upper z plane */
        { 1, -1,  1}, { 1,  0,  1}, { 1,  1,  1},       /* Upper z plane */
        {-1, -1,  0}, {-1,  0,  0}, {-1,  1,  0},       /* In same plane */
        { 0,  1,  0}                                    /* In same plane */
    };

    for (int i = 0; i < n_cellx; i++)
    {
        for (int j = 0; j < n_celly; j++)
        {
            for (int k = 0; k < n_cellz; k++)
            {
                imap = getcell_index(i, j, k);
                for (int count = 0; count < 13; count++)
                {
                    xI = i + shift[count][0];
                    yI = j + shift[count][1];
                    zI = k + shift[count][2];
                    
                    cell_map[count][imap] = getcell_index(xI, yI, zI);

                }
            }
        }
        
    }
    
}


/* --------- construct neighbour list using cell list --------- */
void construct_cell_list()
{
    int A, B;
    int cell_index, L_index, k;         // neigh cell index L type 
    int nbc;                            // neighbour cell (1 to 4)
    int par_index1, par_index2;         // particle index

    double dr2, r_cutoff2;

    NEIB cell_list[N_cell];


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
        for (int j = 0; j < max_neigh; j++)
        {
            cell_list[i].list[j] = 0;
        }
    }


    /* ~~~~~~~~~ distribute particles to cell ~~~~~~~~~ */
    for (int i = 0; i < N; i++)
    {
        cell_index = floor( x[i]/cell_size ) + floor ( y[i]/cell_size ) *  n_cellx
                    + floor(z[i]/cell_size) * n_cellx * n_celly;
        /* cell_index(k) = i + j * n_cellx + k * n_cellx * n_celly */

        k = cell_index;

        if (k<0 || k >= N_cell)
        {
            printf("Some mistake occurs during calculating cell index. \n");
            exit(1);
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

                A = pt_kind[par_index1];
                B = pt_kind[par_index2];
                sigma2 = sigma[A][B] * sigma[A][B];

                r_cutoff2 = r_c2 * sigma2 + skin_depth2;

                if (dr2 < r_cutoff2)
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

            for (nbc = 0; nbc < 13; nbc++)       // 13 neighbour box for 3d
            {
                L_index = cell_map[nbc][i];     // select index of a neighbour box

                for (int m = 0; m < cell_list[L_index].max_counter; m++)
                { 
                    // select another particles of that neighbour box
                    par_index2 = cell_list[L_index].list[m]; 

                    dr2 = particle_distance(par_index1, par_index2);

                    A = pt_kind[par_index1];
                    B = pt_kind[par_index2];
                    sigma2 = sigma[A][B] * sigma[A][B];
                    r_cutoff2 = r_c2 * sigma2 + skin_depth2;

                    if (dr2 < r_cutoff2)
                    {
                        neighbour[par_index1].list[neighbour[par_index1].max_counter] = par_index2;
                        neighbour[par_index1].max_counter ++;
                    } 
                }
            }
            
        }  
        
    }
    
}


/* ------ Verlet list method to costruct neighbour list ------- */
void construct_verlet_list()
{
    int A, B;
    double dr2, r_cutoff2;

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
            dr2 = particle_distance(i, j);

            A = pt_kind[i];
            B = pt_kind[j];
            sigma2 = sigma[A][B] * sigma[A][B];
            r_cutoff2 = r_c2 * sigma2 + skin_depth2;

            if (dr2 < r_cutoff2)
            {
                neighbour[i].list[neighbour[i].max_counter] = j;
                neighbour[i].max_counter ++;
            }
        }
    }  
}



/*~~~~~ Calculate force of the system ~~~~~*/
void force_calculation()
{
    double dx, dy, dz, dr2;
    double eps, f0;
    int typ_A, typ_B;
    int j, tot_neigh;

    for (int i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0.0;
    
    PE = 0.0;
    virial = 0.0;

    // Change to N in cell list, (N-1) in verlet 
    for (int i = 0; i < N; i++)                 
    {                               
        tot_neigh = neighbour[i].max_counter;
        
        for (int k = 0; k < tot_neigh; k++)
        {
            j = neighbour[i].list[k];   // j: neighbour particle index

            typ_A = pt_kind[i];
            typ_B = pt_kind[j];

            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];

            // minimum image condition
            if(dx > half_L) dx = dx - L;
            else if(dx < - half_L) dx = dx + L;

            if(dy > half_L) dy = dy - L;
            else if(dy < - half_L) dy = dy + L;

            if(dz > half_L) dz = dz - L;
            else if(dz < - half_L) dz = dz + L;

            dr2 = dx*dx + dy*dy + dz*dz;

            sigma2 = sigma[typ_A][typ_B] * sigma[typ_A][typ_B];

            if ( dr2 < sigma2 / 100.0 ) 
            {
                printf("Particles came so close. \n");
                exit(0);
            }

            if ( dr2 < (r_c * r_c * sigma2) )
            {
                sigma_by_r2 = sigma2/dr2;
                sigma_by_r6 = sigma_by_r2 * sigma_by_r2 * sigma_by_r2;

                eps = epsilon[typ_A][typ_B];

                f0 = 8.0 * eps * ( (6.0 * sigma_by_r6) *
                        (sigma_by_r6 - 0.50)/dr2 - C2/(sigma2) );
                
                fx[i] = fx[i] + f0 * dx;
                fx[j] = fx[j] - f0 * dx;

                fy[i] = fy[i] + f0 * dy;
                fy[j] = fy[j] - f0 * dy;

                fz[i] = fz[i] + f0 * dz;
                fz[j] = fz[j] - f0 * dz;

                PE += 4.0 * eps * ( sigma_by_r6 * (sigma_by_r6 - 1.0) 
                        + C0 + C2 / sigma_by_r2 );
                virial += f0 * dr2;
                
            }
        }   
    }

    PE = PE / (double) N;   //PE per particle
}



/*~~~~~~~~~ Initilize Cell list and Calculate force(NCM) ~~~~~~~~~~~*/
void Initilize_force_using_cell()
{
    maps();
    construct_cell_list();
    // construct_verlet_list();
    printf(" Cell List constructed succesfully. \n");
    printf("------------------------------------\n");

    force_calculation();
    printf(" Initial Force calculated succesfully. \n");
    printf("-------------------------------------\n");

    Temperature();
    printf(" Initial Temperature of system = %.6f\n", T_inst);
    printf("-------------------------------------\n");

}





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
    Cell List Neighbor Index Table (27 cells, 13 neighbors each)

    Each row corresponds to a cell (0 to 26).
    Each column represents the 13 neighboring cell indices.

    ------------------------------------------------------------
    |  C  |  N1  N2  N3  N4  N5  N6  N7  N8  N9  N10 N11 N12 N13
    ------------------------------------------------------------
    |  0  |  17  11  14  15   9  12  16  10  13   8   2   5   3
    |  1  |  26  20  23  24  18  21  25  19  22  17  11  14  12
    |  2  |   8   2   5   6   0   3   7   1   4  26  20  23  21
    |  3  |  11  14  17   9  12  15  10  13  16   2   5   8   6
    |  4  |  20  23  26  18  21  24  19  22  25  11  14  17  15
    |  5  |   2   5   8   0   3   6   1   4   7  20  23  26  24
    |  6  |  14  17  11  12  15   9  13  16  10   5   8   2   0
    |  7  |  23  26  20  21  24  18  22  25  19  14  17  11   9
    |  8  |   5   8   2   3   6   0   4   7   1  23  26  20  18
    |  9  |  15   9  12  16  10  13  17  11  14   6   0   3   4
    | 10  |  24  18  21  25  19  22  26  20  23  15   9  12  13
    | 11  |   6   0   3   7   1   4   8   2   5  24  18  21  22
    | 12  |   9  12  15  10  13  16  11  14  17   0   3   6   7
    | 13  |  18  21  24  19  22  25  20  23  26   9  12  15  16
    | 14  |   0   3   6   1   4   7   2   5   8  18  21  24  25
    | 15  |  12  15   9  13  16  10  14  17  11   3   6   0   1
    | 16  |  21  24  18  22  25  19  23  26  20  12  15   9  10
    | 17  |   3   6   0   4   7   1   5   8   2  21  24  18  19
    | 18  |  16  10  13  17  11  14  15   9  12   7   1   4   5
    | 19  |  25  19  22  26  20  23  24  18  21  16  10  13  14
    | 20  |   7   1   4   8   2   5   6   0   3  25  19  22  23
    | 21  |  10  13  16  11  14  17   9  12  15   1   4   7   8
    | 22  |  19  22  25  20  23  26  18  21  24  10  13  16  17
    | 23  |   1   4   7   2   5   8   0   3   6  19  22  25  26
    | 24  |  13  16  10  14  17  11  12  15   9   4   7   1   2
    | 25  |  22  25  19  23  26  20  21  24  18  13  16  10  11
    | 26  |   4   7   1   5   8   2   3   6   0  22  25  19  20
    ------------------------------------------------------------
*/
