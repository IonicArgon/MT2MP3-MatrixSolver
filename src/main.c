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
    
    // create a vector x filled with 1s
    double *x = (double *)malloc(A->num_cols * sizeof(double));
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 1.0;
    }

    // create a vector y to store the result of the matrix-vector product
    double *y = (double *)malloc(A->num_rows * sizeof(double));

    // compute the matrix-vector product
    spmv_csr(A, x, y);

    // print the result
    printf("Result of the matrix-vector product:\n");
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("\t%f\n", y[i]);
    }

    // clean up
    CSR_free(A);
    free(x);
    free(y);
    free(A);
}