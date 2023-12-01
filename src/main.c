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

    // set up x and b for jacobi method
    double *x = (double *)malloc(A.num_cols * sizeof(double));
    double *b = (double *)malloc(A.num_rows * sizeof(double));

    // initialize x and b
    for (int i = 0; i < A.num_cols; i++)
    {
        x[i] = 0.0;
    }

    for (int i = 0; i < A.num_rows; i++)
    {
        b[i] = 1.0;
    }

    // solve the linear system
    solver_iter_jacobi(&A, b, x, 100);

    // print the solution
    printf("Solution:\n");
    for (int i = 0; i < A.num_cols; i++)
    {
        printf("%f\n", x[i]);
    }

    // compute the residual
    double residual = compute_residual(&A, b, x);
    printf("Residual: %f\n", residual);

    // free memory
    free(A.csr_data);
    free(A.col_ind);
    free(A.row_ptr);
    free(x);
    free(b);

    return 0;
}