#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "MPC_simplified_qpoases.h"
#define N 2500// Number of data points

int main(int argc, char const *argv[])
{
    FILE *fp;
    double *input = (double*)calloc(N * 160, sizeof(double));

    double *z_op, *x0, *pd, *Q11, *Q22, *R, *umax;
    double z_opt[152];
    int exit_flag[1];
    double u0[1];

    int iter_count[N];

    clock_t start_time, end_time;
    double t_time;

    fp = fopen("data.bin", "rb");
    fread(input, sizeof(double), N * 160, fp);
    fclose(fp);

    
    for (int i = 0; i < N; i++)
    {
        start_time = clock();
        z_op = input + i*160 + 1;
        x0 = z_op + 152;
        pd = x0 + 2;
        Q11 = pd + 1;
        Q22 = Q11 + 1;
        R = Q22 + 1;
        umax = R + 1;

        MPC_simplified_call(
            z_op, 
            x0,
            pd,
            Q11,
            Q22,
            R,
            umax, 
            z_opt, 
            exit_flag, 
            u0,
            iter_count + i);
        end_time = clock();
        t_time = (double)(end_time - start_time) / (double)CLOCKS_PER_SEC;
        printf("%f\n", t_time * 1e3);
    }
    

    

    // printf("Avg compute time: %f us\n", t_time / N * 1e6);

    fp = fopen("data_iter.csv", "w");
    for (int i = 0; i < N; i++)
        fprintf(fp, "%d\n", iter_count[i]);
    fclose(fp);

    return 0;
}
