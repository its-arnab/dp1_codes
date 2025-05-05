/*-----------------------------------------------------*
 *        NVE simulation(PBC) with cell list           *
 *-----------------------------------------------------*/


#include "global.h"
#include "subroutine.h"



int main(int argc, char *argv[])
{
    // Check if the correct number of arguments is provided
    if (argc != 3) {
        fprintf(stderr, "Usage: %s [new|rerun] run_id\n", argv[0]);
        return EXIT_FAILURE;
    }

    int is_fresh_run = strcmp(argv[1], "new") == 0;
    const char *run_id = argv[2];

    char base_path[128], movie_path[200];
    snprintf(base_path, sizeof(base_path), "data_files/T_%.5f/%s", T_res, run_id);
    snprintf(movie_path, sizeof(movie_path), "%s/movie_data", base_path);

    if (is_fresh_run)
    {
        if (directory_exists(base_path))
        {
            printf("[Info] Cleaning directory for new run. \n");
            clean_directory(base_path);
            create_path(movie_path);
            printf("[Info] Creating movie path. \n");
        }
        else
        {
            printf("[Info] Creating directory for new run: %s\n", base_path);
            create_path(base_path);
            create_path(movie_path);
        }
    }
    else
    {
        if (!directory_exists(base_path))
        {
            fprintf(stderr, "[Error] Cannot rerun. Directory does not exist: %s\n", base_path);
            return EXIT_FAILURE;
        }
        printf("[Info] Appending to existing directory: %s\n", base_path);
    }


    // Open file in appropriate mode
    char file_path[256];
    snprintf(file_path, sizeof(file_path), "%s/T_P_KE_PE_E_data", base_path);

    FILE *fid1 = fopen(file_path, is_fresh_run ? "w" : "a");
    if (!fid1) {
        perror("fopen failed");
        return EXIT_FAILURE;
    }


    // mode 0 to re-run, 1 to new run
    Initialize_Program(is_fresh_run); 


    Initilize_force_using_cell();

    for (int fow_run = 0; fow_run < (int) 1e4; fow_run++)
    {
        integrate_nvt();
        if (fow_run % 10000 == 0)
        {
            printf("completed step: %d\n", fow_run);
            fflush(stdout); 
            save_final_config(base_path);
        }

        double sumpx = 0.0;
        double sumpy = 0.0;
        double sumpz = 0.0;
        for (int i = 0; i < N; i++)
        {
            sumpx += vx[i];
            sumpy += vy[i];
            sumpz += vz[i];
        }
        
        if (abs(sumpx) > 1e-5 || abs(sumpy) > 1e-5 || abs(sumpz) > 1e-5)
        {
            printf("completed step: %d\n", fow_run);
            fflush(stdout); 
            break;
        }
        
    }

    
    printf("simulation started ... \n");

    int index = 0;

    for(int step = 0; step <= max_step; step++)
    {
        // writing positions
        if (step == save_time[index])
        {
            save_data_infile(step, movie_path);
            index++;
        }
        
        if (step % 1000 == 0)
        save_final_config(base_path);

        // writing T, PE, KE, E
        if (step % 10 == 0)
        {
            compute_H_eff();

            pressure = rho * T_inst + virial / (3.0 * vol);
        
            fprintf(fid1, "%.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\n",
                            T_inst, pressure, KE, PE, (PE+KE), H_eff);
        }

        integrate_nvt();
        
    }

    // Change to TRUE to save final_config
    if (is_fresh_run)
    {
        FILE *fid3;
        char final_config_path[256];
        snprintf(final_config_path, sizeof(final_config_path), "%s/final_config", base_path);
        fid3 = fopen(final_config_path, "w");
        if( fid3 == NULL  ) return(1);

        for (int i = 0; i < M_chain; i++)
        {
            fprintf(fid3, "%.12lf \t %.12lf \t %.12lf \n", eta[i], p_eta[i], Q[i]);
        }

        for( int i = 0; i < N; i++ )
        {
            fprintf(fid3, "%.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %d\n", 
                    x[i], y[i], z[i], vx[i], vy[i], vz[i], pt_kind[i]);
        }
        fclose(fid3); 
    }


    deallocate_arrays();

    fclose(fid1);

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



