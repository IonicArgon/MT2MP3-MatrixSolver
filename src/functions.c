#include "functions.h"

// --- bare minimum project requirements ---

typedef struct
{
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

// this function performs a sparse matrix vector multiplication
// since CSRMatrix only stores a triangular matrix (we store half of
// the symmetric matrix), we calculate our result like so:
// 1. y1 = A*x (ignoring the diagonal)
// 2. y2 = A^T*x (ignoring the diagonal)
// 3. y = y1 + y2 + (diagonal elements of A)*x
void spmv_csr(CSRMatrix *A, CSRMatrix *AT, const double *x, double *y)
{
    // if AT is NULL, we do a regular matrix-vector multiplication
    if (AT == NULL)
    {
        for (int i = 0; i < A->num_rows; i++)
        {
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                y[i] += A->csr_data[j] * x[A->col_ind[j]];
            }
        }
    }

    // if AT is not NULL, then y = (A_strictly_lower_triangular)*x + AT*x
    else if (AT != NULL)
    {
        // calculate for the strictly lower triangular part first
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

        // calculate for the transpose
        for (int i = 0; i < AT->num_rows; i++)
        {
            for (int j = AT->row_ptr[i]; j < AT->row_ptr[i + 1]; j++)
            {
                y[i] += AT->csr_data[j] * x[AT->col_ind[j]];
            }
        }
    }
}

// --- end of bare minimum project requirements ---

// --- extra stuff ---

// this function checks if two doubles are equal within a certain tolerance
bool fuzzy_equals(double a, double b, double epsilon)
{
    return fabs(a - b) < epsilon;
}

// this function computes the norm of a vector
double compute_norm(double *x, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
}

// this function computes the residual of a matrix
double compute_residual(CSRMatrix *A, CSRMatrix *AT, const double *b, const double *x)
{
    // create a vector to store the result of A*x
    double *Ax = (double *)malloc(A->num_cols * sizeof(double));

    // check if memory allocation was successful
    if (Ax == NULL)
    {
        printf("Error: memory allocation failed\n");
        return 0.0;
    }

    // initialize Ax to 0
    for (int i = 0; i < A->num_cols; i++)
    {
        Ax[i] = 0.0;
    }

    // calculate Ax
    spmv_csr(A, AT, x, Ax);

    // calculate r = Ax - b
    double *r = (double *)malloc(A->num_cols * sizeof(double));

    // check if memory allocation was successful
    if (r == NULL)
    {
        printf("Error: memory allocation failed\n");
        return 0.0;
    }

    for (int i = 0; i < A->num_cols; i++)
    {
        r[i] = Ax[i] - b[i];
    }

    // calculate the norm of r
    double norm = compute_norm(r, A->num_cols);

    // free memory
    free(Ax);
    free(r);

    return norm;
}

// --- end of extra stuff ---

// --- solving stuff ---

// jacobi preconditioner will check the diagonal elements of A and
// perform row swaps if necessary to make sure the diagonal elements
// are non-zero
void diagonal_checker(CSRMatrix *A, double *diagonal)
{
    // check for zeros in the diagonal and swap rows if needed
    int swapped_rows[A->num_rows];
    for (int i = 0; i < A->num_rows; i++)
    {
        swapped_rows[i] = 0;
    }

    for (int i = 0; i < A->num_rows; i++)
    {
        if (diagonal[i] == 0.0)
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
                    diagonal[j] = 0.0;
                }

                for (int j = 0; j < A->num_rows; j++)
                {
                    for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                    {
                        if (A->col_ind[k] == j)
                        {
                            diagonal[j] = A->csr_data[k];
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

// jacobi method solver
void solver_iter_jacobi(CSRMatrix *A, CSRMatrix *AT, double *b, double *x, const int max_iter, double threshold, bool diagonal_check)
{
    // create a vector to store the diagonal elements of A
    double *diagonal = (double *)malloc(A->num_rows * sizeof(double));

    // check if memory allocation was successful
    if (diagonal == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // fill the vector with the diagonal elements of A
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] == i)
            {
                diagonal[i] = A->csr_data[j];
            }
        }
    }

    // check if we need to do a diagonal check
    if (diagonal_check == true)
    {
        printf("Warning: Checking diagonal entries for zeros. The original matrix might not be preserved.\n");
        diagonal_checker(A, diagonal);
        printf("Diagonal check complete.\n");
    }

    // check if the diagonal elements of A are all non-zero
    for (int i = 0; i < A->num_rows; i++)
    {
        if (diagonal[i] == 0.0)
        {
            printf("Error: Matrix is singular and therefore cannot be solved.\n");
            return;
        }
    }

    // initialize x to 0
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
    }

    // make a copy of x to store the previous iteration
    double *x_prev = (double *)malloc(A->num_cols * sizeof(double));

    // check if memory allocation was successful
    if (x_prev == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // if we're not given the transpose, we can just do a normal jacobi iteration
    if (AT == NULL)
    {
        for (int i = 0; i < max_iter; i++)
        {
            // calculate x_prev = x
            for (int j = 0; j < A->num_cols; j++)
            {
                x_prev[j] = x[j];
            }

            // jacobi
            for (int j = 0; j < A->num_rows; j++)
            {
                double sum = 0.0;
                for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                {
                    if (A->col_ind[k] != j)
                    {
                        sum += A->csr_data[k] * x_prev[A->col_ind[k]];
                    }
                }
                x[j] = (b[j] - sum) / diagonal[j];
            }

            // check if we've converged
            double residual = compute_residual(A, AT, b, x);
            if (residual < threshold)
            {
                printf("\nConverged after %d iterations.\n", i + 1);
                break;
            }

            printf("\rIteration: %d", i + 1);
            printf(" | Residual: %e", residual);
            fflush(stdout);
        }
    }

    // otherwise, then we have to consider the tranpose as well
    else if (AT != NULL)
    {
        for (int i = 0; i < max_iter; i++)
        {
            // calculate x_prev = x
            for (int j = 0; j < A->num_cols; j++)
            {
                x_prev[j] = x[j];
            }

            // jacobi
            for (int j = 0; j < A->num_rows; j++)
            {
                double sum = 0.0;
                for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                {
                    if (A->col_ind[k] != j)
                    {
                        sum += A->csr_data[k] * x_prev[A->col_ind[k]];
                    }
                }

                for (int k = AT->row_ptr[j]; k < AT->row_ptr[j + 1]; k++)
                {
                    if (AT->col_ind[k] != j)
                    {
                        sum += AT->csr_data[k] * x_prev[AT->col_ind[k]];
                    }
                }

                x[j] = (b[j] - sum) / diagonal[j];
            }

            // check if we've converged
            double residual = compute_residual(A, AT, b, x);
            if (residual < threshold)
            {
                printf("\nConverged after %d iterations.\n", i + 1);
                break;
            }

            printf("\rIteration: %d", i + 1);
            printf(" | Residual: %e", residual);
            fflush(stdout);
        }
    }

    // free the diagonal vector
    free(diagonal);

    printf("\nSolver complete.\n");
    printf("\n");
}

// successive over-relaxation
void solver_iter_SOR(CSRMatrix *A, CSRMatrix *AT, double *b, double *x, const int max_iter, double threshold, double omega, bool diagonal_check)
{
    // extract our diagonals
    double *diagonal = (double *)malloc(A->num_rows * sizeof(double));

    // check if memory allocation was successful
    if (diagonal == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // fill the vector with the diagonal elements of A
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] == i)
            {
                diagonal[i] = A->csr_data[j];
            }
        }
    }

    // check if we need to do a diagonal check
    if (diagonal_check == true)
    {
        printf("Warning: Checking diagonal entries for zeros. The original matrix might not be preserved.\n");
        diagonal_checker(A, diagonal);
        printf("Diagonal check complete.\n");
    }

    // check if the diagonal elements of A are all non-zero
    for (int i = 0; i < A->num_rows; i++)
    {
        if (diagonal[i] == 0.0)
        {
            printf("Error: Matrix is singular and therefore cannot be solved.\n");
            return;
        }
    }

    // initialize x to 0
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
    }

    // make a copy of x to store the previous iteration
    double *x_prev = (double *)malloc(A->num_cols * sizeof(double));

    // check if memory allocation was successful
    if (x_prev == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // if we're not given the transpose, we can just do a normal SOR
    if (AT == NULL)
    {
        for (int i = 0; i < max_iter; i++)
        {
            // calculate x_prev = x
            for (int j = 0; j < A->num_cols; j++)
            {
                x_prev[j] = x[j];
            }

            // SOR
            for (int j = 0; j < A->num_rows; j++)
            {
                double s1 = 0.0;
                double s2 = 0.0;
                for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                {
                    if (A->col_ind[k] < j)
                    {
                        s1 += A->csr_data[k] * x[A->col_ind[k]];
                    }
                    else if (A->col_ind[k] > j)
                    {
                        s2 += A->csr_data[k] * x_prev[A->col_ind[k]];
                    }
                }
                x[j] = (1 - omega) * x_prev[j] + (omega / diagonal[j]) * (b[j] - s1 - s2);
            }

            // check if we've converged
            double residual = compute_residual(A, AT, b, x);
            if (residual < threshold)
            {
                printf("\nConverged after %d iterations.\n", i + 1);
                break;
            }

            printf("\rIteration: %d", i + 1);
            printf(" | Residual: %e", residual);
            fflush(stdout);
        }
    }

    // otherwise if we're given the transpose, we must consider it symmetric
    else if (AT != NULL)
    {
        for (int i = 0; i < max_iter; i++)
        {
            // calculate x_prev = x
            for (int j = 0; j < A->num_cols; j++)
            {
                x_prev[j] = x[j];
            }

            // SOR
            for (int j = 0; j < A->num_rows; j++)
            {
                double s1 = 0.0;
                double s2 = 0.0;
                for (int k = A->row_ptr[j]; k < A->row_ptr[j + 1]; k++)
                {
                    if (A->col_ind[k] < j)
                    {
                        s1 += A->csr_data[k] * x[A->col_ind[k]];
                    }
                    else if (A->col_ind[k] > j)
                    {
                        s2 += A->csr_data[k] * x_prev[A->col_ind[k]];
                    }
                }

                for (int k = AT->row_ptr[j]; k < AT->row_ptr[j + 1]; k++)
                {
                    if (AT->col_ind[k] < j)
                    {
                        s1 += AT->csr_data[k] * x[AT->col_ind[k]];
                    }
                    else if (AT->col_ind[k] > j)
                    {
                        s2 += AT->csr_data[k] * x_prev[AT->col_ind[k]];
                    }
                }

                x[j] = (1 - omega) * x_prev[j] + (omega / diagonal[j]) * (b[j] - s1 - s2);
            }

            // check if we've converged
            double residual = compute_residual(A, AT, b, x);
            if (residual < threshold)
            {
                printf("\nConverged after %d iterations.\n", i + 1);
                break;
            }

            printf("\rIteration: %d", i + 1);
            printf(" | Residual: %e", residual);
            fflush(stdout);
        }
    }
}

// -- end of solving stuff --

// --- matrix specific functions ---

// this function raw prints a CSRMatrix
void CSR_raw_print(const CSRMatrix *A, bool print_values)
{
    printf("Raw print of CSRMatrix:\n");
    printf("\t- Number of rows: %d\n", A->num_rows);
    printf("\t- Number of columns: %d\n", A->num_cols);
    printf("\t- Number of non-zero elements: %d\n", A->num_non_zeros);

    // if print_values is false, we don't print the values of col_ind and row_ptr
    if (print_values == false)
    {
        return;
    }

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
    printf("%*s", left_spacing, "+-");
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
            int found = 0;
            for (int k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++)
            {
                if (A->col_ind[k] == j)
                {
                    printf("%*lf", max_width + 1, A->csr_data[k]);
                    found = 1;
                    break;
                }
            }
            if (found == 0)
            {
                printf("%*d", max_width + 1, 0);
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
    bool upper = false;
    bool lower = false;
    char triangular = 'N';

    // check if matrix is upper triangular
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] < i)
            {
                upper = false;
                break;
            }
            else
            {
                upper = true;
            }
        }
        if (upper == false)
        {
            break;
        }
    }

    // check if matrix is lower triangular
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            if (A->col_ind[j] > i)
            {
                lower = false;
                break;
            }
            else
            {
                lower = true;
            }
        }
        if (lower == false)
        {
            break;
        }
    }

    // set triangularity
    if (upper == true)
    {
        triangular = 'U';
    }
    else if (lower == true)
    {
        triangular = 'L';
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
    for (int i = 0; i <= A->num_rows; i++)
    {
        A->row_ptr[i] = 0;
    }

    // Count the number of non-zero elements in each row
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        A->row_ptr[elements[i].row + 1]++;
    }

    // Calculate the cumulative sum for row_ptr
    for (int i = 0; i < A->num_rows; i++)
    {
        A->row_ptr[i + 1] += A->row_ptr[i];
    }

    // free the temporary array
    free(elements);
}

// this function swaps two rows of a CSRMatrix
void CSR_row_swap(CSRMatrix *A, int row1, int row2)
{
    // create temporary arrays
    double *temp_csr_data = (double *)malloc(A->num_non_zeros * sizeof(double));
    int *temp_col_ind = (int *)malloc(A->num_non_zeros * sizeof(int));
    int *temp_row_ptr = (int *)malloc((A->num_rows + 1) * sizeof(int));

    // check if memory allocation was successful
    if (temp_csr_data == NULL || temp_col_ind == NULL || temp_row_ptr == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

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

// check if a matrix is strictly diagonally dominant
bool CSR_strictly_diagonally_dominant(const CSRMatrix *A, const CSRMatrix *AT)
{
    // if no AT, we do a regular check
    if (AT == NULL)
    {
        for (int i = 0; i < A->num_rows; i++)
        {
            double sum = 0.0;
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                if (A->col_ind[j] != i)
                {
                    sum += fabs(A->csr_data[j]);
                }
            }
            if (fabs(A->csr_data[A->row_ptr[i] + i]) <= sum)
            {
                return false;
            }
        }
    }

    // if we're an AT, we have a symmetric matrix
    else if (AT != NULL)
    {
        double sum = 0.0;
        for (int i = 0; i < A->num_rows; i++)
        {
            // sum the lower triangular part
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                if (A->col_ind[j] != i)
                {
                    sum += fabs(A->csr_data[j]);
                }
            }

            // sum the upper triangular part
            for (int j = AT->row_ptr[i]; j < AT->row_ptr[i + 1]; j++)
            {
                if (AT->col_ind[j] != i)
                {
                    sum += fabs(AT->csr_data[j]);
                }
            }

            // now check
            if (fabs(A->csr_data[A->row_ptr[i] + i]) <= sum)
            {
                return false;
            }
        }
    }

    return true;
}

// --- end of matrix specific functions ---