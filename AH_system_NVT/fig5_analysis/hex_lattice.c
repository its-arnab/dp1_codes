#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N  1000

double *x, *y;

/* ---- Make Centered Hexagonal Lattice ---- */
void load_triangular_lattice()
{  
    /* ------------------------------------
       H_n = nth centered hexagonal number
       n = number required for getting H_n
    * ----------------------------------- */
    int n, H_n;
    H_n = N;
    n = ceil(0.5 + sqrt(12 * H_n - 3) / 6.0);  /* from wiki */

    printf("n = %d\n", n);

    int pt_num = 0, exit_flag = 0;
    double x_pos, y_pos;
    double sqrt3_by2 = sqrt(3) / 2.0;
    double a = 1.20;

    for (int i = 0; i < n; i++) {
        y_pos = sqrt3_by2 * a * i;
        for (int j = 0; j < (2 * n - 1 - i); j++) {
            x_pos = (-(2 * n - i - 2) * a) / 2.0 + j * a;

            if (pt_num >= N) {  // Prevent overflow
                exit_flag = 1;
                break;
            }

            x[pt_num] = x_pos;
            y[pt_num] = y_pos;
            pt_num++;

            if (y_pos != 0 && pt_num < N) { // Fix condition and prevent overflow
                x[pt_num] = x_pos;
                y[pt_num] = -y_pos;
                pt_num++;
            }
        }
        if (exit_flag) break;
    }
    
    if (pt_num < N) {
        printf("Warning! Only %d particles placed out of %d.\n", pt_num, N);
    }
}

int main() { 
    x = (double*) malloc(sizeof(double) * N);
    y = (double*) malloc(sizeof(double) * N);

    if (!x || !y) {
        printf("Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    load_triangular_lattice();

    FILE *fid = fopen("output.txt", "w");
    if (fid == NULL) {
        printf("Error opening file.\n");
        free(x);
        free(y);
        return EXIT_FAILURE;
    }

    double L = sqrt(N/0.2);
    for (int i = 0; i < N; i++)
    {
        x[i] += 0.5 * L;
        y[i] += 0.5 * L;
    }
    

    // Print lattice points
    for (int i = 0; i < N; i++) {
        fprintf(fid, "%.12lf %.12lf\n", x[i], y[i]);
    }
    fclose(fid);

    free(x);
    free(y);

    return EXIT_SUCCESS;
}
