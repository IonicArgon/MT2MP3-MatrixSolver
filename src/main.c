#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

#define ITERATIONS 5

int main(int argc, char **argv)
{
    const char *filename = argv[1];

    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);
    print_CSRMatrix(A);
    
    // set up b and x for jacobi method
    double *b = (double *)malloc(sizeof(double) * A->num_rows);
    double *x = (double *)malloc(sizeof(double) * A->num_rows);
    for (int i = 0; i < A->num_rows; i++)
    {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // solve the linear system
    printf("Solving linear system using Jacobi Method...\n");
    solver_iter_jacobi(A, b, x, ITERATIONS, true);

    // compute the residual
    double residual = compute_residual(A, b, x);

    // print result and residual
    printf("Results of Jacobi Method:\n");
    printf("\t- residual: %f\n", residual);
    printf("\t- x: ");
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("%f ", x[i]);
    }
    printf("\n");

    // clean up and reset for gauss-seidel method
    free_CSRMatrix(A);
    free(b);
    free(x);

    A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);

    // set up b and x for gauss-seidel method
    b = (double *)malloc(sizeof(double) * A->num_rows);
    x = (double *)malloc(sizeof(double) * A->num_rows);
    for (int i = 0; i < A->num_rows; i++)
    {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // solve the linear system
    printf("Solving linear system using Gauss-Seidel Method...\n");
    solver_iter_gauss_seidel(A, b, x, ITERATIONS, true);

    // compute the residual
    residual = compute_residual(A, b, x);

    // print result and residual
    printf("Results of Gauss-Seidel Method:\n");
    printf("\t- residual: %f\n", residual);
    printf("\t- x: ");
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("%f ", x[i]);
    }
    printf("\n");

    // clean up
    free_CSRMatrix(A);
    free(b);
    free(x);
    
    return 0;
}