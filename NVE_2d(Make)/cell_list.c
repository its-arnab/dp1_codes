#include "global.h"
#include "subroutine.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * calculate distance between i and j particles   *
 * Input: Particle index i, j                     *
 * Output: Distance^2 b/w the particles.          *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double particle_distance(int i, int j)
{
    double dx, dy, dr2;
    
    dx = x[i] - x[j];
    dy = y[i] - y[j];

    // minimum image condition
    if(dx > 0.50 * L) dx = dx - L;
    else if(dx < -0.50 * L) dx = dx + L;

    if(dy > 0.50 * L) dy = dy - L;
    else if(dy < -0.50 * L) dy = dy + L;

    dr2 = dx*dx + dy*dy;

    return dr2;
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

    double dr, dr2;

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
        cell_index = floor( x[i]/cell_size ) + floor ( y[i]/cell_size ) *  n_cellx;
        /* cell_index(k) = i + j * n_cellx */

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


/* ------ Verlet list method to costruct neighbour list ------- */
void construct_verlet_list()
{
    double dr2;

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
            if (dr2 < rl * rl)
            {
                neighbour[i].list[neighbour[i].max_counter] = j;
                neighbour[i].max_counter ++;
            }
        }
    }  
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
    2D Cell List Neighbor Index Table (9 cells, 4 neighbors each)

    Each column represents a cell (0 to 8).
    Each row (from top to bottom) gives one neighbor for that cell.
    
    ---------------------------------------------------------------------------------
          |  Cell 0  Cell 1  Cell 2  Cell 3  Cell 4  Cell 5  Cell 6  Cell 7  Cell 8  
    ---------------------------------------------------------------------------------
    Row 1 |    1       2       0       4       5       3       7       8       6    
    Row 2 |    4       5       3       7       8       6       1       2       0    
    Row 3 |    3       4       5       6       7       8       0       1       2    
    Row 4 |    5       3       4       8       6       7       2       0       1    
    ---------------------------------------------------------------------------------
*/

