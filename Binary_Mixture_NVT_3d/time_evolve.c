#include "global.h"
#include "subroutine.h"



/* ----------Check temperature at any instace ------- */
void Temperature()
{
    double sumv2 = 0;

    for(int i=0; i<N; i++)
    {
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    
    T_inst = sumv2 / (3.0 * N);

    KE = 0.5 * sumv2 / (double) N;  // KE per particles

}



/*----- Apply Periodic Boundary Condition -----*/
void Apply_PBC()
{
    for (int i = 0; i < N; i++)
    {
        if ( x[i] >= L )  x[i] = x[i] - L;
        else if ( x[i] < 0.0 ) x[i] = x[i] + L;

        if ( y[i] >= L )  y[i] = y[i] - L;
        else if ( y[i] < 0.0 ) y[i] = y[i] + L;

        if ( z[i] >= L )  z[i] = z[i] - L;
        else if ( z[i] < 0.0 ) z[i] = z[i] + L;
    }
    
}



/* ------ Time evolution using Velocity Verlet ------ */
void integrate_nve()
{
    double dx, dy, dz;
    double dr2;

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

        xu[i] = xu[i] + dx;
        yu[i] = yu[i] + dy;
        zu[i] = zu[i] + dz;
    }

    Apply_PBC();

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

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * ~~~~ M Chain Nose-Hoover Thermostat ~~~~ *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void Initilize_NoseHoover()
{
    Q[0] = DOF * T_res * tau * tau;
    for (int i = 1; i < M_chain; i++)
    {
        Q[i] = T_res * tau * tau;
    }

    for (int i = 0; i < M_chain; i++) eta[i] = 0.0;

    for (int i = 0; i < M_chain; i++)
    {
        p_eta[i] = sqrt(Q[i] * T_res) * random_normal();    
    }
}


/*---- Position Update -----*/
void U1_propagator()
{
    double dx, dy, dz;
    double dr2;

    for (int i = 0; i < N; i++)
    {
        dx = vx[i] * dt;
        dy = vy[i] * dt;
        dz = vz[i] * dt;

        x[i] = x[i] + dx;
        y[i] = y[i] + dy;
        z[i] = z[i] + dz;

        xu[i] = xu[i] + dx;
        yu[i] = yu[i] + dy;
        zu[i] = zu[i] + dz;

        dr2 = dx*dx + dy*dy + dz*dz;
        if (dr2 > max_disp2 ) max_disp2 = dr2;
    }

    Apply_PBC();
    
}


/*---- Half-Kick of Velocity ----*/
void U2_propagator()
{
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] + fx[i] * half_dt;
        vy[i] = vy[i] + fy[i] * half_dt;
        vz[i] = vz[i] + fz[i] * half_dt;
    }
    
}


/*---- Velocity Adjuster & eta Updater ----*/
void U3_propagator()
{
    double ratio = p_eta[0] / Q[0];

    for (int i = 0; i < N; i++)
    {
        vx[i] *= exp(- half_dt * ratio);
        vy[i] *= exp(- half_dt * ratio);
        vz[i] *= exp(- half_dt * ratio);
    }

    /* Update extended variable eta */
    for (int i = 0; i < M_chain; i++)
    {
        eta[i] += half_dt * p_eta[i] / Q[i];
    }
    
}


/*----- Run [0, M-1] in fwd direction -----*/
/*----- Run [M-1, 0] in Bwd direction -----*/
void U4_propagator(int J_start, int J_end)
{
    double Gj = 0.0, c;
    int j = J_start;
    int dj = 1;
    if (J_start > J_end) dj = - 1;

    for (j = J_start; (dj > 0 ? j <= J_end : j >= J_end); j = j + dj)
    {
        if (j == 0)
        {
            for (int i = 0; i < N; i++)
            {
                Gj += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
            }
            Gj = Gj - DOF * T_res;
        }
        else
        {
            Gj = (p_eta[j-1] * p_eta[j-1]) / Q[j-1] - T_res;
        }

        if (j == M_chain - 1)
        {
            p_eta[j] = p_eta[j] + 0.25 * dt * Gj;
        }
        else
        {
            c = (p_eta[j+1] / Q[j+1]) * 0.25 * dt;
            p_eta[j] = p_eta[j] * exp(-c);
            p_eta[j] = p_eta[j] + Gj * 0.25 * dt * ((-expm1(-c))/c);
        }
    }
    
}


/*----- Time evolution using NH scheme -----*/
void integrate_nvt()
{
    U4_propagator(M_chain - 1, 0);
    U3_propagator();
    U4_propagator(0, M_chain - 1);

    U2_propagator();
    U1_propagator();

    list_counter ++;
    if ( (2.0 * sqrt(max_disp2) * list_counter) > skin_depth )
    {
        // printf("cell list constructed after step: %d \n", list_counter);
        // construct_verlet_list();
        construct_cell_list();
    }

    force_calculation();

    U2_propagator();

    U4_propagator(M_chain - 1, 0);
    U3_propagator();
    U4_propagator(0, M_chain - 1);
}


/*----- compute H_eff ------*/
void compute_H_eff()
{
    Temperature();

    H_eff = (KE + PE) + eta[0] * T_res * 3.0;
    
    double tt = 0.0;

    for (int i = 1; i < M_chain; i++)
    tt = tt + T_res * eta[i];

    for (int i = 0; i < M_chain; i++)
    tt += 0.5 * p_eta[i] * p_eta[i] / Q[i];

    H_eff = H_eff + tt / (double) N;
}

