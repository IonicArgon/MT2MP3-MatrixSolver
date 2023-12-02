#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// ###########################################################
// Do not change this part
typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);
void spmv_csr(CSRMatrix *A, const double *x, double *y);

// ###########################################################

/* <Here you can add the declaration of functions you need.>
<The actual implementation must be in functions.c>
Here what "potentially" you need:
1. A function called "solver" receiving const CSRMatrix A, double *b, double *x 
and solving the linear system 
2. A function called "compute_residual" to compute the residual like r=Ax-b.
This shows how much x values you found are accurate, but 
printing the whole vector r might not be a good idea. So
3. A function called compute_norm to compute the norm of vector residual
*/

// extra stuff
bool fuzzy_equals(double a, double b, double epsilon);
double compute_norm(double *x, int n);
double compute_residual(CSRMatrix *A, const double *b, const double *x);

// solving stuff
void preconditioner_jacobi_gauss(CSRMatrix *A, double *diagonal);
void solver_iter_jacobi(CSRMatrix *A, const double *b, double *x, const int max_iter, bool precondition);

// matrix specific functions
void CSR_raw_print(const CSRMatrix *A);
void CSR_pretty_print(const CSRMatrix *A);
void CSR_free(CSRMatrix *A);
char CSR_triangular_test(const CSRMatrix *A);
void CSR_transpose(CSRMatrix *A);
void CSR_row_swap(CSRMatrix *A, int row1, int row2);

#endif // FUNCTIONS_H