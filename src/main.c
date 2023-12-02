#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

void export_solution(const char *filename, double *x, int n)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%f\n", x[i]);
    }

    fclose(fp);
}

int main(int argc, char **argv)
{
    // check for correct number of arguments
    if (argc != 2)
    {
        printf("Usage: ./main <filename>\n");
        return 1;
    }

    const char *filename = argv[1];

    // ask for number of iterations and method
    int num_iterations;
    char method;
    int preconditioner;

#ifdef USER_INPUT
    printf("Enter number of iterations: ");
    scanf("%d", &num_iterations);
    printf("Enter method (j for jacobi, g for gauss-seidel): ");
    scanf(" %c", &method);
    printf("Enter preconditioner (0 for none, 1 for Jacobi Preconditioning): ");
    scanf("%d", &preconditioner);
#else

#if defined(ITERATIONS)
    num_iterations = ITERATIONS;
#else
    num_iterations = 100;
#endif

#if defined(JACOBI)
    method = 'j';
#elif defined(GAUSS_SEIDEL)
    method = 'g';
#else
    method = 'j';
#endif

#if defined(PRECONDITIONING)
    preconditioner = 1;
#else
    preconditioner = 0;
#endif

#endif

    // validate method
    if (method != 'j' && method != 'g')
    {
        printf("Invalid method\n");
        return 1;
    }

    // validate preconditioner
    if (preconditioner != 0 && preconditioner != 1)
    {
        printf("Invalid preconditioner\n");
        return 1;
    }

    // read matrix from file
    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);
    // raw_print_CSRMatrix(A);

    // set up b and x for jacobi method
    double *b = (double *)malloc(sizeof(double) * A->num_rows);
    double *x = (double *)malloc(sizeof(double) * A->num_rows);
    for (int i = 0; i < A->num_rows; i++)
    {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // solve the linear system
    if (method == 'j')
    {
        printf("Solving with Jacobi Method...\n");
        solver_iter_jacobi(A, b, x, num_iterations, (bool)preconditioner);
    }
    else
    {
        printf("Solving with Gauss-Seidel Method...\n");
        solver_iter_gauss_seidel(A, b, x, num_iterations, (bool)preconditioner);
    }

    // compute the residual
    double residual = compute_residual(A, b, x);

    // print result and residual
    if (method == 'j')
    {
        printf("Jacobi Method:\n");
    }
    else
    {
        printf("Gauss-Seidel Method:\n");
    }

    printf("\t- residual: %f\n", residual);

    // export solution to file
    export_solution("solution.txt", x, A->num_rows);

    // clean up
    free_CSRMatrix(A);
    free(b);
    free(x);

    return 0;
}