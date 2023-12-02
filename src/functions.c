#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

// --- bare minimum project requirements ---

typedef struct {
    int row;
    int col;
    double val;
} Element;

int compare(const void *a, const void *b)
{
    int row_a = ((Element *)a)->row;
    int row_b = ((Element *)b)->row;
    return (row_a - row_b);
}

// this function reads a matrix market file and stores the data in a CSRMatrix
void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    // open file
    FILE *file = fopen(filename, "r");

    // check if file exists
    if (file == NULL)
    {
        printf("Error: file does not exist\n");
        return;
    }

    // ignore the header (the lines start with '%')
    char line[256];
    do
    {
        char *throwaway = fgets(line, sizeof(line), file);
    } while (line[0] == '%');

    // the first line after that is the matrix size
    // rows, cols, number of non-zero elements
    sscanf(line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

    // allocate memory for the arrays
    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = (int *)malloc((matrix->num_rows + 1) * sizeof(int));

    // check if memory allocation was successful
    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // fill the arrays with zeros
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        matrix->csr_data[i] = 0.0;
        matrix->col_ind[i] = 0;
    }

    for (int i = 0; i < matrix->num_rows + 1; i++)
    {
        matrix->row_ptr[i] = 0;
    }

    // create an array of elements to store the data temporarily
    Element *elements = (Element *)malloc(matrix->num_non_zeros * sizeof(Element));

    // check if memory allocation was successful
    if (elements == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }
    
    // fill temporary array with zeros
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        elements[i].row = 0;
        elements[i].col = 0;
        elements[i].val = 0.0;
    }

    // read the data from the file into the temporary array
    int index = 0;
    while (fgets(line, sizeof(line), file))
    {
        int row, col;
        double val;
        sscanf(line, "%d %d %lf", &row, &col, &val);
        elements[index].row = row - 1;
        elements[index].col = col - 1;
        elements[index].val = val;
        index++;
    }

    // sort the temporary array by row
    qsort(elements, matrix->num_non_zeros, sizeof(Element), compare);

    // write csr data
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        matrix->csr_data[i] = elements[i].val;
        matrix->col_ind[i] = elements[i].col;
    }

    // row pointer data is the number of non-zero elements in each row plus the previous row pointer
    int row = 0;
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        if (elements[i].row == row)
        {
            matrix->row_ptr[row + 1]++;
        }
        else
        {
            row++;
            matrix->row_ptr[row + 1]++;
        }
    }

    for (int i = 0; i < matrix->num_rows; i++)
    {
        matrix->row_ptr[i + 1] += matrix->row_ptr[i];
    }

    // free the temporary array
    free(elements);

    // close the file
    fclose(file);
}

// this function computes the sparse matrix vector multiplication
void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{
    // y = A * x
    for (int i = 0; i < A->num_rows; i++)
    {
        double sum = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            sum += A->csr_data[j] * x[A->col_ind[j]];
        }
        y[i] = sum;
    }
}

// this function computes the residual of Ax-b
double compute_residual(const CSRMatrix *A, double *b, double *x)
{
    double *r = (double *)malloc(A->num_rows * sizeof(double));
    spmv_csr(A, x, r);

    for (int i = 0; i < A->num_rows; i++)
    {
        r[i] = r[i] - b[i];
    }

    double norm = compute_norm(r, A->num_rows);

    free(r);

    return norm;
}

// this function computes the norm of a vector
double compute_norm(double *r, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += r[i] * r[i];
    }
    norm = sqrt(norm);

    return norm;
}

// --- end of bare minimum project requirements ---

// --- here below is just extras ---

// this function checks if two doubles are equal within a certain epsilon
bool fuzzy_equals(double a, double b, double epsilon)
{
    if (fabs(a - b) < epsilon)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// --- end of extras ---

// --- here below is just csr matrix functions ---

// this function swaps two rows in a CSRMatrix
void CSR_row_swap(CSRMatrix *A, int row1, int row2)
{
    // create temporary arrays
    double *temp_csr_data = (double *)malloc(A->num_non_zeros * sizeof(double));
    int *temp_col_ind = (int *)malloc(A->num_non_zeros * sizeof(int));
    int *temp_row_ptr = (int *)malloc((A->num_rows + 1) * sizeof(int));

    // fill temporary arrays with zeros
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        temp_csr_data[i] = 0.0;
        temp_col_ind[i] = 0;
    }

    for (int i = 0; i < A->num_rows + 1; i++)
    {
        temp_row_ptr[i] = 0;
    }

    // copy data from A into temporary arrays until we hit row1
    int temp_csr_index = 0;
    int temp_row_ptr_index = 0;
    while (temp_row_ptr_index < row1)
    {
        for (int i = A->row_ptr[temp_row_ptr_index]; i < A->row_ptr[temp_row_ptr_index + 1]; i++)
        {
            temp_csr_data[temp_csr_index] = A->csr_data[i];
            temp_col_ind[temp_csr_index] = A->col_ind[i];
            temp_csr_index++;
        }
        temp_row_ptr[temp_row_ptr_index + 1] = temp_csr_index;
        temp_row_ptr_index++;
    }

    // copy the data from row2 into the spot where row1 was
    for (int i = A->row_ptr[row2]; i < A->row_ptr[row2 + 1]; i++)
    {
        temp_csr_data[temp_csr_index] = A->csr_data[i];
        temp_col_ind[temp_csr_index] = A->col_ind[i];
        temp_csr_index++;
    }

    temp_row_ptr[temp_row_ptr_index + 1] = temp_csr_index;
    temp_row_ptr_index++;

    // copy the data from A into temporary arrays until we hit row2
    while (temp_row_ptr_index < row2)
    {
        for (int i = A->row_ptr[temp_row_ptr_index]; i < A->row_ptr[temp_row_ptr_index + 1]; i++)
        {
            temp_csr_data[temp_csr_index] = A->csr_data[i];
            temp_col_ind[temp_csr_index] = A->col_ind[i];
            temp_csr_index++;
        }
        temp_row_ptr[temp_row_ptr_index + 1] = temp_csr_index;
        temp_row_ptr_index++;
    }

    // copy the data from row1 into the spot where row2 was
    for (int i = A->row_ptr[row1]; i < A->row_ptr[row1 + 1]; i++)
    {
        temp_csr_data[temp_csr_index] = A->csr_data[i];
        temp_col_ind[temp_csr_index] = A->col_ind[i];
        temp_csr_index++;
    }

    temp_row_ptr[temp_row_ptr_index + 1] = temp_csr_index;
    temp_row_ptr_index++;

    // copy the data from A into temporary arrays until we hit the end
    while (temp_row_ptr_index < A->num_rows)
    {
        for (int i = A->row_ptr[temp_row_ptr_index]; i < A->row_ptr[temp_row_ptr_index + 1]; i++)
        {
            temp_csr_data[temp_csr_index] = A->csr_data[i];
            temp_col_ind[temp_csr_index] = A->col_ind[i];
            temp_csr_index++;
        }
        temp_row_ptr[temp_row_ptr_index + 1] = temp_csr_index;
        temp_row_ptr_index++;
    }

    // copy the data from temporary arrays into A
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        A->csr_data[i] = temp_csr_data[i];
        A->col_ind[i] = temp_col_ind[i];
    }

    for (int i = 0; i < A->num_rows + 1; i++)
    {
        A->row_ptr[i] = temp_row_ptr[i];
    }

    // free temporary arrays
    free(temp_csr_data);
    free(temp_col_ind);
    free(temp_row_ptr);
}

// this function frees the memory allocated for a CSRMatrix
void free_CSRMatrix(CSRMatrix *A)
{
    free(A->csr_data);
    free(A->col_ind);
    free(A->row_ptr);
    A->csr_data = NULL;
    A->col_ind = NULL;
    A->row_ptr = NULL;
    A->num_rows = 0;
    A->num_cols = 0;
    A->num_non_zeros = 0;
}

// --- end of csr matrix functions ---

// --- here below is just solvers and parts of solvers ---

// --- preconditioners

// Jacobi/Gauss-Seidel preconditioner
// ensures that the matrix at the very least contains all non-zero elements
// in the diagonal
void preconditioner_jacobi_gauss(CSRMatrix *A, double *diag)
{
    // check for zeros in the diagonal and swap rows if needed
    int swapped_rows[A->num_rows];
    for (int i = 0; i < A->num_rows; i++)
    {
        swapped_rows[i] = 0;
    }

    for (int i = 0; i < A->num_rows; i++)
    {
        if (diag[i] == 0.0)
        {
            // find a row that has a value in csr_data in the same column as the zero
            printf("Warning: zero found at row %d and column %d\n", i, i);
            printf("Swapping rows to fill the zero in...\n");
            printf("Note: this does not guarantee that the matrix will be solvable!\n\n");
            int row_to_swap = -1;
            for (int j = 0; j < A->num_rows; j++)
            {
                for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                {
                    // if the row we're looking at has already been swapped, skip it
                    if (swapped_rows[j])
                    {
                        continue;
                    }

                    if (A->col_ind[k] == i && fuzzy_equals(A->csr_data[k], 0.0, 0.0000000000001) == false)
                    {
                        row_to_swap = j;
                        break;
                    }
                }
                if (row_to_swap != -1)
                {
                    break;
                }
            }

            // if we found a row to swap, then swap it
            if (row_to_swap != -1 && !swapped_rows[row_to_swap] && !swapped_rows[i])
            {
                printf("Swapping row %d with row %d\n\n", i, row_to_swap);
                if (i > row_to_swap)
                {
                    int temp = i;
                    i = row_to_swap;
                    row_to_swap = temp;
                }
                swapped_rows[i] = 1;

                CSR_row_swap(A, i, row_to_swap);

                // update the diagonal array
                for (int j = 0; j < A->num_rows; j++)
                {
                    diag[j] = 0.0;
                }

                for (int j = 0; j < A->num_rows; j++)
                {
                    for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                    {
                        if (A->col_ind[k] == j)
                        {
                            diag[j] = A->csr_data[k];
                        }
                    }
                }
            }
            else
            {
                return;
            }
        }
    }
}

// --- solvers

// Jacobi method solver
void solver_iter_jacobi(CSRMatrix *A, double *b, double *x, int max_iter, bool row_swap)
{
    // create a new array for diagonal elements, initialize it to zero
    double *diag = (double *)malloc(A->num_rows * sizeof(double));
    for (int i = 0; i < A->num_rows; i++)
    {
        diag[i] = 0.0;
    }

    // then copy over any diagonal elements from A
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] == i)
            {
                diag[i] = A->csr_data[j];
            }
        }
    }

    if (row_swap)
    {
        printf("Warning: row swapping is enabled. The original matrix will not be preserved.\n\n");
        preconditioner_jacobi_gauss(A, diag);
    }

    // initialize x to zero
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
    }

    // if we still have zeros in the diagonal, then we can't solve the system
    for (int i = 0; i < A->num_rows; i++)
    {
        if (diag[i] == 0.0)
        {
            printf("Error: Matrix is singular\n");
            return;
        }
    }

    // iterate and compute x
    // residual calculation is done in the main function
    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int i = 0; i < A->num_rows; i++)
        {
            double sum = 0.0;
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                if (A->col_ind[j] != i)
                {
                    sum += A->csr_data[j] * x[A->col_ind[j]];
                }
            }
            x[i] = (b[i] - sum) / diag[i];
        }
    }

    // free the diagonal array
    free(diag);
}

// Gauss-Seidel method solver
void solver_iter_gauss_seidel(CSRMatrix *A, double *b, double *x, int max_iter, bool row_swap)
{
    // create a new array for diagonal elements, initialize it to zero
    double *diag = (double *)malloc(A->num_rows * sizeof(double));
    for (int i = 0; i < A->num_rows; i++)
    {
        diag[i] = 0.0;
    }

    // then copy over any diagonal elements from A
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] == i)
            {
                diag[i] = A->csr_data[j];
            }
        }
    }

    if (row_swap)
    {
        printf("Warning: row swapping is enabled. The original matrix will not be preserved.\n\n");
        preconditioner_jacobi_gauss(A, diag);
    }

    // initialize x to zero
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
    }

    // if we still have zeros in the diagonal, then we can't solve the system
    for (int i = 0; i < A->num_rows; i++)
    {
        if (diag[i] == 0.0)
        {
            printf("Error: Matrix is singular\n");
            return;
        }
    }

    // iterate and compute x
    // residual calculation is done in the main function
    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int i = 0; i < A->num_rows; i++)
        {
            double s1 = 0.0;
            double s2 = 0.0;

            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                if (A->col_ind[j] < i)
                {
                    s1 += A->csr_data[j] * x[A->col_ind[j]];
                }
                else if (A->col_ind[j] > i)
                {
                    s2 += A->csr_data[j] * x[A->col_ind[j]];
                }
            }

            x[i] = (b[i] - s1 - s2) / diag[i];
        }
    }

    // free the diagonal array
    free(diag);
}

// --- end of solvers ---

// --- here below is just extra project requirements ---

// pretty print the matrix
void print_CSRMatrix(const CSRMatrix *A)
{
    printf("Pretty print of CSR Matrix:\n");
    printf("\t- num_rows: %d\n", A->num_rows);
    printf("\t- num_cols: %d\n", A->num_cols);
    printf("\t- num_non_zeros: %d\n\n", A->num_non_zeros);

    // print the whole matrix, row by row, including the zeros
    // something like:
    //    0 1 2 3 4 5 6 7 8 9
    // 0 [0 0 0 0 0 0 0 0 0 0]
    // 1 [0 0 0 0 0 0 0 0 0 0]
    // ...

    // calculate the spacing needed for each element
    int max_element_width = 0;
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        int element_width = snprintf(NULL, 0, "%lf", A->csr_data[i]);
        if (element_width > max_element_width)
        {
            max_element_width = element_width;
        }
    }

    // calculate the spacing on the left for each row based on the
    // number of rows in the matrix
    int max_row_width = snprintf(NULL, 0, "%d", A->num_rows);
    int row_spacing = max_row_width + 3;

    // print the column row first
    printf("%*s", row_spacing - 1, "");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("%*d", max_element_width + 1, i);
    }
    printf("\n");

    // print the matrix
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("%*d [", max_row_width, i);
        for (int j = 0; j < A->num_cols; j++)
        {
            int found = 0;
            for (int k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++)
            {
                if (A->col_ind[k] == j)
                {
                    printf("%*lf", max_element_width + 1, A->csr_data[k]);
                    found = 1;
                    break;
                }
            }
            if (!found)
            {
                printf("%*d", max_element_width + 1, 0);
            }
        }
        printf("]\n");
    }
    printf("\n");
}

// print the matrix in raw format
void raw_print_CSRMatrix(const CSRMatrix *A)
{
    printf("Raw print of CSR Matrix:\n");
    printf("\t- num_rows: %d\n", A->num_rows);
    printf("\t- num_cols: %d\n", A->num_cols);
    printf("\t- num_non_zeros: %d\n\n", A->num_non_zeros);

    printf("csr_data:\n");
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        printf("%lf ", A->csr_data[i]);
    }
    printf("\n\n");

    printf("col_ind:\n");
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        printf("%d ", A->col_ind[i]);
    }
    printf("\n\n");

    printf("row_ptr:\n");
    for (int i = 0; i < A->num_rows + 1; i++)
    {
        printf("%d ", A->row_ptr[i]);
    }
    printf("\n\n");
}

// --- end of extra project requirements ---