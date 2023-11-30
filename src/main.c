#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char **argv)
{
    const char *filename = argv[1];

    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    print_CSRMatrix(&A);

    // testing spmv_csr
    double *x = (double *)malloc(A.num_cols * sizeof(double));
    double *y = (double *)malloc(A.num_rows * sizeof(double));

    for (int i = 0; i < A.num_cols; i++)
    {
        x[i] = 1.0;
    }

    spmv_csr(&A, x, y);

    printf("y = \n");
    for (int i = 0; i < A.num_rows; i++)
    {
        printf("%lf\n", y[i]);
    }
    printf("\n");

    // clean up x and y
    free(x);
    free(y);

    // testing solver
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    double *x_sol = (double *)malloc(A.num_cols * sizeof(double));

    for (int i = 0; i < A.num_rows; i++)
    {
        b[i] = 1.0;
    }

    solver_iter_jacobi(&A, b, x_sol, 10000);

    printf("x_sol = \n");
    for (int i = 0; i < A.num_cols; i++)
    {
        printf("%lf\n", x_sol[i]);
    }
    printf("\n");

    // print the residual
    double residual = compute_residual(&A, b, x_sol);
    printf("residual = %lf\n", residual);

    // clean up b and x_sol
    free(b);
    free(x_sol);

    // clean up csr matrix
    free(A.csr_data);
    free(A.col_ind);
    free(A.row_ptr);

    return 0;
}