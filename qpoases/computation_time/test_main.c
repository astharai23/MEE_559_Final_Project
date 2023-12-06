#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
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

    struct timeval st, et;
    double t_time;

    fp = fopen("data.bin", "rb");
    fread(input, sizeof(double), N * 160, fp);
    fclose(fp);

    // start_time = clock();
    for (int i = 0; i < N; i++)
    {
        gettimeofday(&st,NULL);
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
            u0);
        gettimeofday(&et,NULL);
        t_time = (double)(((et.tv_sec - st.tv_sec) * 1000000) + (et.tv_usec - st.tv_usec));
        printf("%f\n", t_time);
    }
    // end_time = clock();

    // t_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // printf("Avg compute time: %f us\n", t_time / N * 1e6);

    return 0;
}
