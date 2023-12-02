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

    // print out program information
    printf("Running with the following parameters:\n");
    printf("\t- filename: %s\n", filename);
    printf("\t- num_iterations: %d\n", num_iterations);
    printf("\t- method: ");
    if (method == 'j')
    {
        printf("Jacobi method\n");
    }
    else
    {
        printf("Gauss-Seidel method\n");
    }
    printf("\t- preconditioner: ");
    if (preconditioner == 0)
    {
        printf("None\n");
    }
    else
    {
        printf("Jacobi preconditioning\n\n");
    }

    // read matrix from file
    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);

#if defined(PRINT)
    printf("Note: You should recompile this without the PRINT flag for large matrices.\n");

    if (PRINT == 2)
    {
        print_CSRMatrix(A);
    }
    else if (PRINT == 1)
    {
        raw_print_CSRMatrix(A);
    }
    else
    {
        raw_print_CSRMatrix(A);
    }

#endif

    // set up b and x for jacobi method
    double *b = (double *)malloc(sizeof(double) * A->num_rows);
    double *x = (double *)malloc(sizeof(double) * A->num_rows);
    for (int i = 0; i < A->num_rows; i++)
    {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // solve the linear system
    printf("Solving linear system...\n");
    if (method == 'j')
    {
        solver_iter_jacobi(A, b, x, num_iterations, (bool)preconditioner);
    }
    else
    {
        solver_iter_gauss_seidel(A, b, x, num_iterations, (bool)preconditioner);
    }

    // compute the residual
    double residual = compute_residual(A, b, x);

    // print result and residual
    printf("Result:\n");
    printf("\t- Residual: %e\n", residual);

    // export solution to file
#if defined(PRINT)
    printf("\t- Solution of x: \n");
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("\t\t %f\n", x[i]);
    }
    printf("\n");
#else
    printf("Exporting solution to file...\n");
    export_solution("solution.txt", x, A->num_rows);
#endif

    // clean up
    free_CSRMatrix(A);
    free(b);
    free(x);
    free(A);

    return 0;
}