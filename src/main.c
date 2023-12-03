#include "functions.h"

int main(int argc, char const *argv[])
{
    // check if the user provided the correct number of arguments
    if (argc != 2)
    {
        printf("Usage: %s <matrix market file>\n", argv[0]);
        return 1;
    }

    char method;
    int max_iter;
    int precondition;
    double threshold;

#if defined (USER_INPUT)
    // get the method from the user
    printf("Please select a method:\n");
    printf("\t1. Jacobi (j)\n");
    printf("\t2. Gauss-Seidel (g)\n");
    scanf("%c", &method);

    // get the maximum number of iterations from the user
    printf("Please enter the maximum number of iterations:\n");
    scanf("%d", &max_iter);

    // get the preconditioner from the user
    printf("Please select a preconditioner:\n");
    printf("\t1. None (0)\n");
    printf("\t2. Jacobi (1)\n");
    scanf("%d", &precondition);

    // check if the user entered a valid method
    if (method != 'j' && method != 'g')
    {
        printf("Invalid method\n");
        return 1;
    }

    // check if the user entered a valid preconditioner
    if (precondition != 0 && precondition != 1)
    {
        printf("Invalid preconditioner\n");
        return 1;
    }

    // check if the user entered a valid maximum number of iterations
    if (max_iter < 1)
    {
        printf("Invalid maximum number of iterations\n");
        return 1;
    }
#else
    // set the method based on defined constants at compile time
    #if defined (JACOBI)
        method = 'j';
    #elif defined (GAUSS_SEIDEL)
        method = 'g';
    #endif

    // set the maximum number of iterations based on defined constants at compile time
    #if defined (MAX_ITER)
        max_iter = MAX_ITER;
    #endif

    // set the threshold based on defined constants at compile time
    #if defined (THRESHOLD)
        threshold = THRESHOLD;
    #endif

    // set the preconditioner based on defined constants at compile time
    #if defined (PRECONDITIONING)
        precondition = 1;
    #else
        precondition = 0;
    #endif
#endif

    // print program information
    printf("Method: ");
    if (method == 'j')
    {
        printf("Jacobi\n");
    }
    else if (method == 'g')
    {
        printf("Gauss-Seidel\n");
    }

    printf("Number of iterations: %d\n", max_iter);
    printf("Threshold: %e\n", threshold);

    printf("Preconditioner: ");
    if (precondition == 0)
    {
        printf("None\n");
    }
    else if (precondition == 1)
    {
        printf("Jacobi\n");
    }
    printf("\n");


    // create a CSRMatrix
    const char *filename = argv[1];
    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);

    // create a vector x and b to prepare for solving
    double *x = (double *)malloc(A->num_cols * sizeof(double));
    double *b = (double *)malloc(A->num_cols * sizeof(double));

    // initialize x and b
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
        b[i] = 0.0;
    }

    // print the matrix based on compile time constants
    #if defined (PRINT)
    if (PRINT == 1)
    {
        CSR_raw_print(A);
    }
    else if (PRINT == 2)
    {
        CSR_pretty_print(A);
    }
    #endif

    // check if the matrix is triangular
    char triangular = CSR_triangular_test(A);
    if (triangular == 'N')
    {
        printf("The matrix is not triangular\n");
    }
    else if (triangular == 'L')
    {
        printf("The matrix is lower triangular\n");
    }
    else if (triangular == 'U')
    {
        printf("The matrix is upper triangular\n");
    }
    printf("\n");

    // show that the matrix multiplication works with an array of ones
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 1.0;
    }

    // compute the matrix-vector product
    spmv_csr(A, x, b);

    // print the result if PRINT is defined
    #if defined (PRINT)
    printf("\nMatrix-vector product:\n");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("\t%f\n", b[i]);
    }
    printf("\n");
    
    // otherwise, write the result to a file
    #else
    FILE *fp = fopen("smvp_output.txt", "w");
    for (int i = 0; i < A->num_cols; i++)
    {
        fprintf(fp, "%f\n", b[i]);
    }
    fclose(fp);
    #endif

    // reset x and b
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
        b[i] = 1.0;
    }

    // solve the system
    printf("Solving the system...\n");
    if (method == 'j')
    {
        solver_iter_jacobi(A, b, x, max_iter, threshold, precondition);
    }
    else if (method == 'g')
    {
        printf("Gauss-Seidel is not implemented yet\n");
        //solver_iter_gauss_seidel(A, b, x, max_iter, threshold, precondition);
    }

    // compute the residual
    double residual = compute_residual(A, b, x);
    printf("\nResidual: %e\n", residual);

    // print the solution if PRINT is defined
    #if defined (PRINT)
    printf("\nSolution:\n");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("\t%f\n", x[i]);
    }

    // otherwise, write the solution to a file
    #else
    fp = fopen("solution.txt", "w");
    for (int i = 0; i < A->num_cols; i++)
    {
        fprintf(fp, "%f\n", x[i]);
    }
    fclose(fp);
    #endif

    // clean up
    CSR_free(A);
    free(x);
    free(b);
    free(A);
}