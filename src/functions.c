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
void spmv_csr(CSRMatrix *A, const double *x, double *y)
{
    // first check if the matrix is triangular
    char triangular = CSR_triangular_test(A);

    // if the matrix is not triangular or is upper triangular, we can just do a normal matrix-vector multiplication
    if (triangular == 'N' || triangular == 'U')
    {
        for (int i = 0; i < A->num_rows; i++)
        {
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
            {
                y[i] += A->csr_data[j] * x[A->col_ind[j]];
            }
        }
        return;
    }

    // lower triangular matrices are defined as symmetric, skew, or hermitian
    // matrices in MTX format so we calculate using A*x + A^T*x + (diagonal elements of A)*x

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
double compute_residual(CSRMatrix *A, const double *b, const double *x)
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
    spmv_csr(A, x, Ax);

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

// jacobi method solver
void solver_iter_jacobi(CSRMatrix *A, const double *b, double *x, const int max_iter, double threshold, bool precondition)
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

    // check if we need to use a preconditioner
    if (precondition == true)
    {
        printf("Warning: preconditioning is enabled. The original matrix will not be preserved.\n");
        preconditioner_jacobi_gauss(A, diagonal);
        printf("Preconditioning complete.\n");
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

    // check the kind of matrix we're dealing with
    char triangular = CSR_triangular_test(A);

    // if it's non triangular or upper triangular, we can just do a normal jacobi iteration
    if (triangular == 'N' || triangular == 'U')
    {
        for (int iter = 0; iter < max_iter; iter++)
        {
            for (int i = 0; i < A->num_rows; i++)
            {
                double s = 0.0;

                for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
                {
                    if (A->col_ind[j] != i)
                    {
                        s += A->csr_data[j] * x[A->col_ind[j]];
                    }
                }

                x[i] = (b[i] - s) / diagonal[i];
            }

            // check if we've reached the threshold
            double residual = compute_residual(A, b, x);
            if (compute_residual(A, b, x) < threshold)
            {
                printf("\nResidual reached threshold. Stopping iterations.");
                break;
            }

            printf("\rIteration: %d", iter);
            printf(" | Residual: %.12lf", residual);
            fflush(stdout);
        }
    }

    // if it's lower triangular, the matrix is symmetric, skew, or hermitian
    // we have iterate but take into account both halfs of the matrix (reflection across x axis)
    else if (triangular == 'L')
    {
        // precompute A^T
        CSRMatrix *A_transpose = (CSRMatrix *)malloc(sizeof(CSRMatrix));

        // check if memory allocation was successful
        if (A_transpose == NULL)
        {
            printf("Error: memory allocation failed\n");
            return;
        }

        // copy A over to A^T
        A_transpose->num_rows = A->num_rows;
        A_transpose->num_cols = A->num_cols;
        A_transpose->num_non_zeros = A->num_non_zeros;

        A_transpose->csr_data = (double *)malloc(A_transpose->num_non_zeros * sizeof(double));
        A_transpose->col_ind = (int *)malloc(A_transpose->num_non_zeros * sizeof(int));
        A_transpose->row_ptr = (int *)malloc((A_transpose->num_rows + 1) * sizeof(int));

        // check if memory allocation was successful
        if (A_transpose->csr_data == NULL || A_transpose->col_ind == NULL || A_transpose->row_ptr == NULL)
        {
            printf("Error: memory allocation failed\n");
            return;
        }

        // now copy over the data
        memcpy(A_transpose->csr_data, A->csr_data, A_transpose->num_non_zeros * sizeof(double));
        memcpy(A_transpose->col_ind, A->col_ind, A_transpose->num_non_zeros * sizeof(int));
        memcpy(A_transpose->row_ptr, A->row_ptr, (A_transpose->num_rows + 1) * sizeof(int));

        // transpose A^T
        CSR_transpose(A_transpose);

        for (int iter = 0; iter < max_iter; iter++)
        {
            for (int i = 0; i < A->num_rows; i++)
            {
                double s = 0.0;

                for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
                {
                    if (A->col_ind[j] != i)
                    {
                        s += A->csr_data[j] * x[A->col_ind[j]];
                    }
                }

                // same for transpose
                for (int j = A_transpose->row_ptr[i]; j < A_transpose->row_ptr[i + 1]; j++)
                {
                    if (A_transpose->col_ind[j] != i)
                    {
                        s += A_transpose->csr_data[j] * x[A_transpose->col_ind[j]];
                    }
                }

                x[i] = (b[i] - s) / diagonal[i];
            }

            // check if we've reached the threshold
            double residual = compute_residual(A, b, x);
            if (compute_residual(A, b, x) < threshold)
            {
                printf("\nResidual reached threshold. Stopping iterations.");
                break;
            }

            printf("\rIteration: %d", iter);
            printf(" | Residual: %.12lf", residual);
            fflush(stdout);
        }

        // free A^T
        CSR_free(A_transpose);
        free(A_transpose);
    }

    // free the diagonal vector
    free(diagonal);

    printf("\nSolver complete.\n");
    printf("\n");
}

/*
// gauss-seidel method solver
void solver_iter_gauss_seidel(CSRMatrix *A, const double *b, double *x, const int max_iter, double threshold, bool precondition)
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

    // check if we need to use a preconditioner
    if (precondition == true)
    {
        printf("Warning: preconditioning is enabled. The original matrix will not be preserved.\n");
        preconditioner_jacobi_gauss(A, diagonal);
    }

    // initialize x to 0
    for (int i = 0; i < A->num_cols; i++)
    {
        x[i] = 0.0;
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

    // check if the matrix is triangular
    char triangular = CSR_triangular_test(A);

    // iterate and compute x
    // if it's non triangular or upper triangular, we can just do a normal gauss-seidel iteration
    if (triangular == 'N' || triangular == 'U')
    {
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

                x[i] = (b[i] - s1 - s2) / diagonal[i];
            }
        }
    }

    // if it's lower triangular, the matrix is symmetric, skew, or hermitian
    // we have iterate but take into account both halfs of the matrix (reflection across x axis)
    else if (triangular == 'L')
    {
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

                CSR_transpose(A);
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
                CSR_transpose(A);

                x[i] = (b[i] - s1 - s2) / diagonal[i];
            }
        }
    }

    // free the diagonal vector
    free(diagonal);
}
*/

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
bool CSR_strictly_diagonally_dominant(const CSRMatrix *A)
{
    // check if the matrix is triangular
    char triangular = CSR_triangular_test(A);

    // if it's not triangular or is upper triangular, we can just do a normal check
    if (triangular == 'N' || triangular == 'U')
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

    // if it's lower triangular, the matrix is symmetric, skew, or hermitian
    // so we also need the transpose
    if (triangular == 'L')
    {
        // precompute A^T
        CSRMatrix *A_transpose = (CSRMatrix *)malloc(sizeof(CSRMatrix));

        // check if memory allocation was successful
        if (A_transpose == NULL)
        {
            printf("Error: memory allocation failed\n");
            return false;
        }

        // copy A over to A^T
        A_transpose->num_rows = A->num_rows;
        A_transpose->num_cols = A->num_cols;
        A_transpose->num_non_zeros = A->num_non_zeros;

        A_transpose->csr_data = (double *)malloc(A_transpose->num_non_zeros * sizeof(double));
        A_transpose->col_ind = (int *)malloc(A_transpose->num_non_zeros * sizeof(int));
        A_transpose->row_ptr = (int *)malloc((A_transpose->num_rows + 1) * sizeof(int));

        // check if memory allocation was successful
        if (A_transpose->csr_data == NULL || A_transpose->col_ind == NULL || A_transpose->row_ptr == NULL)
        {
            printf("Error: memory allocation failed\n");
            return false;
        }

        // now copy over the data
        memcpy(A_transpose->csr_data, A->csr_data, A_transpose->num_non_zeros * sizeof(double));
        memcpy(A_transpose->col_ind, A->col_ind, A_transpose->num_non_zeros * sizeof(int));
        memcpy(A_transpose->row_ptr, A->row_ptr, (A_transpose->num_rows + 1) * sizeof(int));

        // transpose A^T
        CSR_transpose(A_transpose);

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
            for (int j = A_transpose->row_ptr[i]; j < A_transpose->row_ptr[i + 1]; j++)
            {
                if (A_transpose->col_ind[j] != i)
                {
                    sum += fabs(A_transpose->csr_data[j]);
                }
            }

            // now check
            if (fabs(A->csr_data[A->row_ptr[i] + i]) <= sum)
            {
                return false;
            }
        }

        // free A^T
        CSR_free(A_transpose);
        free(A_transpose);
    }

    return true;
}

// --- end of matrix specific functions ---