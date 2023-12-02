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
    matrix->csr_data = (double *)malloc(2 * matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = (int *)malloc(2 * matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = (int *)malloc((matrix->num_rows + 1) * sizeof(int));

    // check if memory allocation was successful
    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // fill the arrays with zeros
    for (int i = 0; i < 2 * matrix->num_non_zeros; i++)
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

// this function performs a sparse matrix vector multiplication
// since CSRMatrix only stores a triangular matrix (we store half of
// the symmetric matrix), we calculate our result like so:
// 1. y1 = A*x (ignoring the diagonal)
// 2. y2 = A^T*x (ignoring the diagonal)
// 3. y = y1 + y2 + (diagonal elements of A)*x
void spmv_csr(CSRMatrix *A, const double *x, double *y)
{
    // calculate y1
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] != i)
            {
                y[i] += A->csr_data[j] * x[A->col_ind[j]];
            }
        }
    }

    // calculate y2
    CSR_transpose(A);
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] != i)
            {
                y[i] += A->csr_data[j] * x[A->col_ind[j]];
            }
        }
    }
    CSR_transpose(A);

    // calculate (diagonal elements of A)*x
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] == i)
            {
                y[i] += A->csr_data[j] * x[A->col_ind[j]];
            }
        }
    }
}

// --- end of bare minimum project requirements ---


// --- matrix specific functions ---

// this function raw prints a CSRMatrix
void CSR_raw_print(const CSRMatrix *A)
{
    printf("Raw print of CSRMatrix:\n");
    printf("\t- Number of rows: %d\n", A->num_rows);
    printf("\t- Number of columns: %d\n", A->num_cols);
    printf("\t- Number of non-zero elements: %d\n", A->num_non_zeros);

    printf("\t- csr_data:\n\t\t");
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        printf("%lf ", A->csr_data[i]);
    }
    printf("\n");

    printf("\t- col_ind:\n\t\t");
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        printf("%d ", A->col_ind[i]);
    }
    printf("\n");

    printf("\t- row_ptr:\n\t\t");
    for (int i = 0; i < A->num_rows + 1; i++)
    {
        printf("%d ", A->row_ptr[i]);
    }
    printf("\n");
}

// this function pretty prints a CSRMatrix
// note that this only prints one triangle of the matrix
void CSR_pretty_print(const CSRMatrix *A)
{
    printf("Pretty print of CSRMatrix:\n");
    printf("\t- Number of rows: %d\n", A->num_rows);
    printf("\t- Number of columns: %d\n", A->num_cols);
    printf("\t- Number of non-zero elements: %d\n", A->num_non_zeros);

    // calculate max element width
    int max_width = 0;
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        int element_width = snprintf(NULL, 0, "%lf", A->csr_data[i]);
        if (element_width > max_width)
        {
            max_width = element_width;
        }
    }

    // calculate spacing on left side of matrix based on rows
    int max_row_width = snprintf(NULL, 0, "%d", A->num_rows);
    int left_spacing = max_row_width + 3;

    // print column indices
    printf("%*s", left_spacing, "");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("%*d", max_width + 1, i);
    }
    printf("\n");

    // print top border
    printf("%*s", left_spacing, "");
    for (int i = 0; i < A->num_cols; i++)
    {
        for (int j = 0; j < max_width + 1; j++)
        {
            printf("-");
        }
    }
    printf("\n");

    // print matrix including zeros
    int index = 0;
    for (int i = 0; i < A->num_rows; i++)
    {
        printf("%*d |", max_row_width, i);
        for (int j = 0; j < A->num_cols; j++)
        {
            if (j == A->col_ind[index])
            {
                printf("%*lf", max_width + 1, A->csr_data[index]);
                index++;
            }
            else
            {
                printf("%*s", max_width + 1, "");
            }
        }
        printf("\n");
    }
    printf("\n");
}


// this function frees the memory allocated for a CSRMatrix
void CSR_free(CSRMatrix *A)
{
    free(A->csr_data);
    free(A->col_ind);
    free(A->row_ptr);
}

// this function checks the triangularity of a matrix
// returns 'U' if upper triangular, 'L' if lower triangular, 'N' if neither
char CSR_triangular_test(const CSRMatrix *A)
{
    // check if matrix is square
    if (A->num_rows != A->num_cols)
    {
        printf("Error: matrix is not square\n");
        return 'N';
    }

    char triangular = 'N';

    // check if matrix is upper triangular
    // for each row, check if there is an element before the diagonal element
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] < i)
            {
                triangular = 'N';
            }
            else if (A->col_ind[j] == i)
            {
                triangular = 'U';
                break;
            }
        }
    }

    // check if matrix is lower triangular
    // for each row, check if there is an element after the diagonal element
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] > i)
            {
                triangular = 'N';
            }
            else if (A->col_ind[j] == i)
            {
                triangular = 'L';
                break;
            }
        }
    }

    return triangular;
}

// this function transposes a CSRMatrix
void CSR_transpose(CSRMatrix *A)
{
    // create temporary array of elements
    Element *elements = (Element *)malloc(A->num_non_zeros * sizeof(Element));

    // check if memory allocation was successful
    if (elements == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // fill temporary array with zeros
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        elements[i].row = 0;
        elements[i].col = 0;
        elements[i].val = 0.0;
    }

    // iterate through the matrix and store the elements in the temporary array
    // the row and column of each element are swapped
    int index = 0;
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            elements[index].row = A->col_ind[j];
            elements[index].col = i;
            elements[index].val = A->csr_data[j];
            index++;
        }
    }

    // sort the temporary array by row
    qsort(elements, A->num_non_zeros, sizeof(Element), compare);

    // write csr data
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        A->csr_data[i] = elements[i].val;
        A->col_ind[i] = elements[i].col;
    }

    // Initialize row_ptr to 0
    for (int i = 0; i <= A->num_rows; i++) {
        A->row_ptr[i] = 0;
    }

    // Count the number of non-zero elements in each row
    for (int i = 0; i < A->num_non_zeros; i++) {
        A->row_ptr[elements[i].row + 1]++;
    }

    // Calculate the cumulative sum for row_ptr
    for (int i = 0; i < A->num_rows; i++) {
        A->row_ptr[i + 1] += A->row_ptr[i];
    }

    // free the temporary array
    free(elements);
}

// --- end of matrix specific functions ---