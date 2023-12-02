#include "functions.h"

int main(int argc, char const *argv[])
{
    // check if the user provided the correct number of arguments
    if (argc != 2)
    {
        printf("Usage: %s <matrix market file>\n", argv[0]);
        return 1;
    }

    // create a CSRMatrix
    const char *filename = argv[1];
    CSRMatrix *A = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ReadMMtoCSR(filename, A);

    // print the matrix
    CSR_pretty_print(A);

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
    
    // create a vector x and b to prepare for solving
    double *x = (double *)malloc(A->num_cols * sizeof(double));
    double *b = (double *)malloc(A->num_cols * sizeof(double));

    // initialize x and b
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
        b[i] = 1.0;
    }

    // solve the system
    solver_iter_jacobi(A, b, x, 1000, true);

    // compute the residual
    double residual = compute_residual(A, b, x);
    printf("\nResidual: %e\n", residual);

    // print the solution
    printf("\nSolution:\n");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("\t%f\n", x[i]);
    }

    // clean up
    CSR_free(A);
    free(x);
    free(b);
    free(A);
}