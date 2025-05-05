#include "global.h"


void md_step_generator()
{
    int base = 1;
    int blocksize = 80000;
    int numblocks = max_step / blocksize;
    int configperblock = 17;  // up to 2^16 = 65536 fits in a block of 80000

    int totconfig = numblocks * configperblock;
    save_time = malloc(totconfig * sizeof(int));

    if (save_time == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numblocks; i++)
    {
        for (int j = 0; j < configperblock; j++)
        {
            // compute 2^j using left shift
            int time = base * (1 << j) + i * blocksize;

            // safety check: ensure time does not exceed max_step
            if (time >= max_step) break;

            save_time[i * configperblock + j] = time;
        }
    }
}




/*==========================*/
/* Random Number Generator  */
/*==========================*/


/*---- Normal Random Number Generator ----*/
double random_normal()
{
    return sqrt(-2 * log(drand48())) * cos(2.0 * 3.14159265358979 * drand48());
}
