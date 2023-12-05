#include "functions.h"
#include <sys/time.h>

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
    int diagonal_check;
    int precondition;
    double threshold;
    double omega;

#if defined (USER_INPUT)
    // get the method from the user
    printf("Please select a method:\n");
    printf("\t1. Jacobi (j)\n");
    printf("\t2. Gauss-Seidel (g)\n");
    scanf("%c", &method);

    // get the maximum number of iterations from the user
    printf("Please enter the maximum number of iterations:\n");
    scanf("%d", &max_iter);

    // ask if the user wants a diagonal check
    printf("Perform diagonal check?\n");
    printf("\t1. None (0)\n");
    printf("\t2. Yes (1)\n");
    scanf("%d", &diagonal_check);

    // get the threshold from the user
    printf("Please enter the threshold:\n");
    scanf("%lf", &threshold);

    // check if the user entered a valid method
    if (method != 'j' && method != 'g')
    {
        printf("Invalid method\n");
        return 1;
    }

    // check if the diagonal check is valid
    if (diagonal_check != 0 && diagonal_check != 1)
    {
        printf("Invalid choice for diagonal check\n");
        return 1;
    }

    // check if the user entered a valid threshold
    if (threshold < 0.0)
    {
        printf("Invalid threshold\n");
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
    #elif defined (SOR)
        method = 's';

        // set the omega based on defined constants at compile time
        #if defined (OMEGA)
            omega = OMEGA;
        #endif

    #endif

    // set the maximum number of iterations based on defined constants at compile time
    #if defined (MAX_ITER)
        max_iter = MAX_ITER;
    #endif

    // set the threshold based on defined constants at compile time
    #if defined (THRESHOLD)
        threshold = THRESHOLD;
    #endif

    // set the diagonal check based on defined constants at compile time
    #if defined (DIAGONAL_CHECK)
        diagonal_check = 1;
    #else
        diagonal_check = 0;
    #endif
#endif

    // print program information
    printf("Method: ");
    if (method == 'j')
    {
        printf("Jacobi\n");
    }
    else if (method == 's')
    {
        printf("Successive Over-Relaxation\n");
        printf("Omega: %f\n", omega);
    }

    printf("Number of iterations: %d\n", max_iter);
    printf("Threshold: %e\n", threshold);

    printf("Diagonal check: ");
    if (diagonal_check == 0)
    {
        printf("No\n");
    }
    else if (diagonal_check == 1)
    {
        printf("Yes\n");
    }
    printf("\n");

    // get start (wall clock time)
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // create a CSRMatrix
    const char *filename = argv[1];
    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    CSRMatrix *AT = NULL;
    ReadMMtoCSR(filename, A);
    CSR_raw_print(A, false);

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

    // if the matrix is lower triangular, create the transpose
    if (triangular == 'L')
    {
        AT = (CSRMatrix *)malloc(sizeof(CSRMatrix));
        AT->num_cols = A->num_cols;
        AT->num_rows = A->num_rows;
        AT->num_non_zeros = A->num_non_zeros;
        AT->csr_data = (double *)malloc(AT->num_non_zeros * sizeof(double));
        AT->col_ind = (int *)malloc(AT->num_non_zeros * sizeof(int));
        AT->row_ptr = (int *)malloc((AT->num_rows + 1) * sizeof(int));
        memcpy(AT->csr_data, A->csr_data, AT->num_non_zeros * sizeof(double));
        memcpy(AT->col_ind, A->col_ind, AT->num_non_zeros * sizeof(int));
        memcpy(AT->row_ptr, A->row_ptr, (AT->num_rows + 1) * sizeof(int));

        CSR_transpose(AT);
    }

    // check if the matrix is strictly diagonally dominant
    bool sdd = CSR_strictly_diagonally_dominant(A, AT);
    if (sdd)
    {
        printf("The matrix is strictly diagonally dominant\n");
    }
    else
    {
        printf("The matrix is not strictly diagonally dominant\n");
    }
    printf("\n");

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
        CSR_raw_print(A, false);
    }
    else if (PRINT == 2)
    {
        CSR_pretty_print(A);
    }
    #endif

    // show that the matrix multiplication works with an array of ones
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 1.0;
    }

    // compute the matrix-vector product
    printf("Computing the matrix-vector product...\n");
    spmv_csr(A, AT, x, b);
    printf("\n");

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
        fprintf(fp, "%e\n", b[i]);
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
    printf("Solving the system...\n\n");
    if (method == 'j')
    {
        solver_iter_jacobi(A, AT, b, x, max_iter, threshold, diagonal_check);
    }
    else if (method == 's')
    {
        solver_iter_SOR(A, AT, b, x, max_iter, threshold, omega, diagonal_check);
    }

    // compute the residual
    double residual = compute_residual(A, AT, b, x);
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
        fprintf(fp, "%e\n", x[i]);
    }
    fclose(fp);
    #endif

    // clean up
    CSR_free(A);
    free(x);
    free(b);
    free(A);

    if (AT != NULL)
    {
        CSR_free(AT);
        free(AT);
    }

    // get end (wall clock time)
    gettimeofday(&end, NULL);

    // compute the elapsed time
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("\nElapsed time: %f seconds\n", elapsed_time);
}