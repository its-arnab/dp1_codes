#include "global.h"
#include "subroutine.h"


/*~~~~~ Calculate force of the system ~~~~~*/
void force_calculation()
{
    double dx, dy, dz, dr2;
    double f0;

    int j, tot_neigh;

    for (int i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0.0;
    
    PE = 0.0;

    // Change to N in cell list, (N-1) in verlet 
    for (int i = 0; i < N; i++)                 
    {                               
        tot_neigh = neighbour[i].max_counter;
        
        for (int k = 0; k < tot_neigh; k++)
        {
            j = neighbour[i].list[k];   // j: neighbour particle index

            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];

            // minimum image condition
            if(dx > 0.50 * L) dx = dx - L;
            else if(dx < -0.50 * L) dx = dx + L;

            if(dy > 0.50 * L) dy = dy - L;
            else if(dy < -0.50 * L) dy = dy + L;

            if(dz > 0.50 * L) dz = dz - L;
            else if(dz < -0.50 * L) dz = dz + L;

            dr2 = dx*dx + dy*dy + dz*dz;

            if ( dr2 < sigma2 / 100.0 ) 
            {
                printf("Particles came so close. %d and %d \n", i, j);
                exit(0);
            }

            if ( dr2 < (r_c * r_c) )
            {
                sigma_by_r2 = sigma2/dr2;
                sigma_by_r6 = sigma_by_r2 * sigma_by_r2 * sigma_by_r2;

                f0 = 8.0 * epsilon * ( (6.0 * sigma_by_r6) *
                                     (sigma_by_r6 - 0.50)/dr2 - C2/(sigma2) );
                
                fx[i] = fx[i] + f0 * dx;
                fx[j] = fx[j] - f0 * dx;

                fy[i] = fy[i] + f0 * dy;
                fy[j] = fy[j] - f0 * dy;

                fz[i] = fz[i] + f0 * dz;
                fz[j] = fz[j] - f0 * dz;

                PE += 4.0 * epsilon * ( sigma_by_r6 * (sigma_by_r6 - 1.0) 
                                                + C0 + C2 / sigma_by_r2 );
                
            }
         }   
     }

    PE = PE / (double) N;   //PE per particle
} 



/* ----------Check temperature at any instace ------- */
void Temperature()
{
    double T, sumv2 = 0;

    for(int i=0; i<N; i++)
    {
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    
    T_inst = sumv2 / (3.0 * N);

    KE = 0.5 * sumv2 / (double) N;  // KE per particles

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


/* ------ Time evolution using Velocity Verlet ------ */
void integrate()
{
    double dx, dy, dz;
    double dr2, max_disp2;

    max_disp2 = 0.0;

    for (int i = 0; i < N; i++)
    {
        // velocity-verlet: step1
        vx[i] = vx[i] + fx[i] * half_dt;
        vy[i] = vy[i] + fy[i] * half_dt;
        vz[i] = vz[i] + fz[i] * half_dt;

        dx = vx[i] * dt;
        dy = vy[i] * dt;
        dz = vz[i] * dt;

        dr2 = dx*dx + dy*dy + dz*dz;
        if (dr2 > max_disp2 ) max_disp2 = dr2;

        x[i] = x[i] + dx;
        y[i] = y[i] + dy;
        z[i] = z[i] + dz;
        
        // applying PBC
        if ( x[i] >= L )  x[i] = x[i] - L;
        else if ( x[i] < 0.0 ) x[i] = x[i] + L;

        if ( y[i] >= L )  y[i] = y[i] - L;
        else if ( y[i] < 0.0 ) y[i] = y[i] + L;

        if ( z[i] >= L )  z[i] = z[i] - L;
        else if ( z[i] < 0.0 ) z[i] = z[i] + L;
    }

    list_counter ++;
    if ( (2.0 * sqrt(max_disp2) * list_counter) > skin_depth )
    {
        // printf("cell list constructed after step: %d \n", list_counter);
        // construct_verlet_list();
        construct_cell_list();
    }

    // Calculate new force
    force_calculation();

    for( int i = 0; i < N; i++)
    {
       vx[i] = vx[i] + fx[i] * half_dt;
       vy[i] = vy[i] + fy[i] * half_dt;
       vz[i] = vz[i] + fz[i] * half_dt;
    }

}
