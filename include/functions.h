#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdbool.h>

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
void spmv_csr(const CSRMatrix *A, const double *x, double *y);

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

// project requirements
void print_CSRMatrix(const CSRMatrix *A);
void raw_print_CSRMatrix(const CSRMatrix *A);
double compute_residual(const CSRMatrix *A, double *b, double *x);
double compute_norm(double *r, int n);

// extras
bool fuzzy_equals(double a, double b, double epsilon);

// a bunch of csr matrix functions
void CSR_row_swap(CSRMatrix *A, int row1, int row2);
void free_CSRMatrix(CSRMatrix *A);

// a bunch of solvers
// jacobi method solver
void solver_iter_jacobi(CSRMatrix *A, double *b, double *x, int max_iter, bool row_swap);
// gauss-seidel method solver
void solver_iter_gauss_seidel(CSRMatrix *A, double *b, double *x, int max_iter, bool row_swap);

#endif // FUNCTIONS_H