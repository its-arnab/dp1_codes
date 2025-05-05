/*-----------------------------------------------------*
 *        NVE simulation(PBC) with cell list           *
 *-----------------------------------------------------*/


#include "global.h"
#include "subroutine.h"


int main()
{
    
    FILE *fid1, *fid2;

    // mode 0 to re-run, 1 to new run
    Initilize_Program(1); 
    

    Initilize_force_using_cell();


    fid1 = fopen("data_files/animation1.xyz", "w");
    if( fid1 == NULL ) return(1);

    fid2 = fopen("data_files/T_PE_KE_E_data", "w");
    if( fid2 == NULL ) return(1);


    printf("simulation started ... \n");
    

    for(int i = 0; i < max_step; i++)
    {
        // writing positions
        // fprintf(fid1, "%d \n", N);
        // fprintf(fid1, "%s \n", " ");

        // for( int j = 0; j < N; j++ )
        // {
        //     fprintf(fid1, "%s \t %f \t %f \t %f \n",
        //             "Ar", x[j], y[j], 0.0);
        // }

        // writing T, PE, KE, E
        if (i % 100 == 0)
        {
            Temperature(); 
            fprintf(fid2, "%d \t %f \t %f \t %f \t %f \n",
                    i, T_inst, KE, PE, (PE+KE));
        }
        integrate();
        
    }

    // Change to TRUE to save final_config
    if (0)
    {
        FILE *fid3;
        fid3 = fopen("data_files/final_config", "w");
        if( fid3 == NULL  ) return(1);
        for( int i = 0; i < N; i++ )
        {
            fprintf(fid3, "%.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\n", 
                    x[i], y[i], z[i], vx[i], vy[i], vz[i]);
        }
        fclose(fid3); 
    }


    deallocate_arrays();

    fclose(fid1);
    fclose(fid2);

    printf("\n ===========================================\n");
    printf("    ✅ SUCCESS! Everything ran smoothly.     \n");
    printf(" ===========================================\n\n");

    return 0;

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Updates in the MD Simulation Code:                                      *
 *                                                                         *
 * 1. The function `construct_cell_list()` has been optimized, reducing    *
 *    computational cost.                                                  *
 * 2. The neighbour list created taking loop cell wise i.e. we take i-th   *
 *    cell, select a particle from that cell and calculate all its         *
 *    neighbour(inter & intra cell)                                        *
 * 3. The integration module also has an automatic update mechanism        *
 *    for the neighbor list when necessary.                                *
 *                                                                         *
 * Benchmarking:                                                           *
 *                                                                         *
 * - Energy vs. time plots are used to compare computational costs of      *
 *   different force calculation methods: brute-force, Verlet list, and    *
 *   cell list. Initially, all methods yield identical results, but        *
 *   differences arise over time due to floating-point round-off errors.   *
 * - Since floating-point operations are not associative, minor            *
 *   discrepancies in results are expected.                                *
 * - Ensure the same initial configuration is used for all comparisons.    *
 *                                                                         *
 * Implementation Notes:                                                   *
 *                                                                         *
 * 1. In the force subroutine, include `i = 1:N` for cell list comparison, *
 *   as opposed to `(N-1)` in the Verlet list, where the last particle     *
 *   has no neighbors.                                                     *
 * 2. The cell list method accounts for intra-cell neighbors, ensuring     *
 *   that even the N-th particle is properly considered.                   *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */



