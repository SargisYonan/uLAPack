/**
 * @file ulapack.c
 * @brief Definitions of uLAPack matrix manipulation functions.
 *
 * @author Sargis Yonan
 * @date July 8, 2017
 *
 * @version 1.0.1
 **/

/*
 * CLib math required for sqrt() function.
 */
#include <math.h>

#include "ulapack.h"

/**
 * @name _ulapack_is_valid_memory
 * Static subroutine to determine if a passed in object has been initialized.
 *
 * @param[in] mat A pointer to an object.
 *
 * @return ulapack_error if the object has not yet been initialized.
 *         ulapack_success if the object passed in has been initialized.
 */
static MatrixError_t _ulapack_is_valid_memory(const Matrix_t *mat) {
    if (!mat) {
        return ulapack_error;
    }

    return ulapack_success;
}

MatrixError_t ulapack_init(
    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t *matrix,
    #else
        Matrix_t **matrix, 
    #endif
    const uint64_t n_rows, const uint64_t n_cols) {

    /*
     * Case when using static memory allocation.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        /*
         * Check to make sure a legitimate matrix object block was given.
         */
        if (!matrix) {
            return ulapack_invalid_argument;
        }

        /*
         * Return an error if a dimension given exceeds the preset maximum 
         * dimensions.
         */
        if (n_rows > ULAPACK_MAX_MATRIX_N_ROWS || 
            n_cols > ULAPACK_MAX_MATRIX_N_COLS) {

            matrix->n_rows = 0;
            matrix->n_cols = 0;

            return ulapack_invalid_argument;
        }

        /*
         * Set the matrix dimensions.
         */
        matrix->n_rows = n_rows;
        matrix->n_cols = n_cols;

        /*
         * If specified, initialize all memory to zeros.
         */
        #ifdef ULAPACK_INITIALIZE_MEMORY
            ulapack_set(matrix, 0);
        #endif

        return ulapack_success;
    /*
     * Case when using dynamic memory allocation.
     */
    #else
        /*
         * Declare loop counters here in case C99 is not used.
         */
        uint64_t row_itor = 0;
        uint64_t col_itor = 0;
        uint64_t free_rows_itor = 0;

        /*
         * Check to make sure the matrix object pointer given is not allocated.
         */
        if (*matrix) {
            return ulapack_invalid_argument;
        }

        /*
         * Allocate memory for overall object and check for successful
         * allocation.
         */
        *matrix = ULAPACK_ALLOC(sizeof(Matrix_t));
        if (!matrix) {
            return ulapack_oom;
        }
        /*
         * Allocate memory for the matrix rows.
         */
        (*matrix)->entry = ULAPACK_ALLOC(sizeof(MatrixEntry_t) * n_rows);
        
        /*
         * Check for successful alloc.
         */
        if (!((*matrix)->entry))
        {
            /*
             * Free the last unusable memory allocated.
             */
            ULAPACK_FREE(*matrix);
            *matrix = NULL;

            return ulapack_oom;
        }

        /*
         * For each row in the matrix, allocate enough columns.
         */
        for(row_itor = 0; row_itor < n_rows; row_itor++) {
            /*
             * Allocate memory for each row.
             */
            (*matrix)->entry[row_itor] = ULAPACK_ALLOC(
                sizeof(MatrixEntry_t) * n_cols);

            /*
             * Check for error in alloc. If an error occurred, free all memory
             * allocated up to this point, and then return with error.
             */
            if (!((*matrix)->entry[row_itor])) {
                for(free_rows_itor = 0; 
                    free_rows_itor < row_itor; 
                    free_rows_itor++) {
                    ULAPACK_FREE((*matrix)->entry[row_itor]);
                }

                ULAPACK_FREE(*matrix);
                *matrix = NULL;

                return ulapack_oom;
            }

            /*
             * If specified, initialize all allocated memory to zeros.
             */
            #ifdef ULAPACK_INITIALIZE_MEMORY
                for (col_itor = 0; col_itor < n_cols; col_itor++) {
                    (*matrix)->entry[row_itor][col_itor] = 0;
                }
            #endif
        }

        /*
         * Set the matrix dimensions.
         */
        (*matrix)->n_rows = n_rows;
        (*matrix)->n_cols = n_cols;

        return ulapack_success;

    #endif
}

#ifdef ULAPACK_USE_DYNAMIC_ALLOC
MatrixError_t ulapack_destroy(Matrix_t *obj) {
    uint64_t row_itor = 0;

    if (obj) {
        if (obj->entry) {
            /*
             * First free and NULL out all dynamically allocated matrix rows.
             */
            for (row_itor = 0; row_itor < obj->n_rows; row_itor++) {
                ULAPACK_FREE(obj->entry[row_itor]);
                obj->entry[row_itor] = NULL;
            }
            
            /*
             * Free and NULL out the matrix entry pointer.
             */
            if (obj->entry) {
                ULAPACK_FREE(obj->entry);
                obj->entry = NULL;
            }
        }

        /*
         * Free and NULLify the actual object.
         */
        ULAPACK_FREE(obj);
        obj = NULL;

        return ulapack_success;
    }

    return ulapack_uninit_obj;
}
#endif

#ifdef ULAPACK_USE_PRINT
MatrixError_t ulapack_print(const Matrix_t * const matrix,
                            FILE *stream) {

    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        
        fprintf(stream, "[ ");
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            fprintf(stream, ULAPACK_PRINT_DELIMITER " ", 
                matrix->entry[row_itor][col_itor]);
        }

        fprintf(stream, "]\n");
    }

    fprintf(stream, "\n");

    return ulapack_success;
}
#endif

MatrixError_t ulapack_edit_entry(Matrix_t * const matrix, 
                                 const uint64_t row, const uint64_t col,
                                 const MatrixEntry_t value) {
    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify that the coordinates to edit are within the initialized dimensions
     * of the matrix object.
     */
    if (row > matrix->n_rows || col > matrix->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * Edit the value of the matrix at (row, col).
     */
    matrix->entry[row][col] = value;

    return ulapack_success;
}

MatrixError_t ulapack_get_entry(const Matrix_t * const matrix, 
                                  const uint64_t row, const uint64_t col,
                                  MatrixEntry_t * const value) {
    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify that the coordinates to edit are within the initialized dimensions
     * of the matrix object.
     */
    if (row > matrix->n_rows || col > matrix->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * The value pointer should not be NULL because it should be written to.
     */
    if (!value) {
        return ulapack_invalid_argument;
    }

    /*
     * Store the value of the matrix at (row, col).
     */
    *value = matrix->entry[row][col];
    return ulapack_success;
}

MatrixError_t ulapack_size(const Matrix_t * const matrix, 
                           uint64_t * const rows, uint64_t * const cols) {
    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    if (!rows || !cols) {
        return ulapack_invalid_argument;
    }

    *rows = matrix->n_rows;
    *cols = matrix->n_cols;

    return ulapack_success;
}

MatrixError_t ulapack_set(Matrix_t * const matrix,
                          const MatrixEntry_t value) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Set the elements of the matrix.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            matrix->entry[row_itor][col_itor] = value;
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_is_equal(const Matrix_t * const A, 
                               const Matrix_t * const B,
                               MatrixError_t * const is_equal) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(A) || !_ulapack_is_valid_memory(B)) {
        return ulapack_uninit_obj;
    }

    /*
     * Check if the dimensions of the operands are identical.
     */
    if (A->n_rows != B->n_rows || A->n_cols != B->n_cols) {
        *is_equal = ulapack_error;
        return ulapack_success;
    }

    /*
     * The is_equal pointer should not be NULL because it should be written to.
     */
    if (!is_equal) {
        return ulapack_invalid_argument;
    }

    /*
     * Compare the entries of the operands. Check if A == B.
     */
    for (row_itor = 0; row_itor < A->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < A->n_cols; col_itor++) {
            /*
             * If a mismatch occurs, set is_equal to ulapack_error, and return
             * ulapack_success to indicate a successful comparison.
             */
            if (A->entry[row_itor][col_itor] != B->entry[row_itor][col_itor]) {
                *is_equal = ulapack_error;
                return ulapack_success;
            }
        }
    }

    /*
     * The matrices are assumed to be identical at this point.
     */
    *is_equal = ulapack_success;
    return ulapack_success;
}

/**
 * @name ulapack_is_vector
 * Check if the input operand is a vector: 1xN or Nx1.
 *
 * @param[in] vector An initialized vector object operand.
 * @param[out] is_vector Pass back if the input is a vector. ulapack_success if
 *             the input is a vector, and ulapack_error if the input is not a
 *             vector.
 *
 * @return ULAPack success code ulapack_success is returned if the input vector
 *         is a initialized and checked for vector dimensions successful.
 */
MatrixError_t ulapack_is_vector(const Matrix_t * const vector, 
                                MatrixError_t * const is_vector) {
    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(vector)) {
        return ulapack_uninit_obj;
    }

    if ((vector->n_rows == 1 && vector->n_cols >= 1) ||
        (vector->n_cols == 1 && vector->n_rows >= 1)) {
            *is_vector = ulapack_success;
    } else {
        *is_vector = ulapack_error;
    }

    return ulapack_success;
}

MatrixError_t ulapack_eye(Matrix_t * const matrix) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Check to ensure that the dimensions of the operands are identical.
     */
    if (matrix->n_rows != matrix->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * Set the diagonal elements to ones, and everything else to zeros.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            
            if (row_itor != col_itor){
                matrix->entry[row_itor][col_itor] = 0;
            } else {
                matrix->entry[row_itor][col_itor] = 1;
            }

        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_sum(const Matrix_t * const matrix,
                          MatrixEntry_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!result) {
        return ulapack_invalid_argument;
    }

    /*
     * Reset result.
     */
    *result = 0;

    /*
     * Add the elements of the matrix.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            *result += matrix->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_add(const Matrix_t * const A, 
                          const Matrix_t * const B,
                          Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(A) || !_ulapack_is_valid_memory(B)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify that the result matrix object has been initialized.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_uninit_obj;
    }

    /*
     * Check to ensure that the dimensions of the operands are identical.
     */
    if (A->n_rows != B->n_rows || A->n_cols != B->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = A->n_rows;
        result->n_cols = A->n_cols;

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of A and B.
     */
    #else
        if (A->n_rows != result->n_rows || A->n_cols != result->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Compute the sum of the operands: A + B, and store the result into the
     * result object.
     */
    for (row_itor = 0; row_itor < A->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < A->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            A->entry[row_itor][col_itor] + B->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_subtract(const Matrix_t * const A, 
                               const Matrix_t * const B,
                               Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(A) || !_ulapack_is_valid_memory(B)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify that the result matrix object has been initialized.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_uninit_obj;
    }

    /*
     * Check to ensure that the dimensions of the operands are identical.
     */
    if (A->n_rows != B->n_rows || A->n_cols != B->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = A->n_rows;
        result->n_cols = A->n_cols;

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of A and B.
     */
    #else
        if (A->n_rows != result->n_rows || A->n_cols != result->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Compute the difference of the operands: A - B, and store the result into 
     * the result object.
     */
    for (row_itor = 0; row_itor < A->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < A->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            A->entry[row_itor][col_itor] - B->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_scale(Matrix_t * const matrix, 
                            const MatrixEntry_t scalar,
                            Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify that the result matrix object has been initialized.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_uninit_obj;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = matrix->n_rows;
        result->n_cols = matrix->n_cols;

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of A and B.
     */
    #else
        if (matrix->n_rows != result->n_rows || 
            matrix->n_cols != result->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Scale the elements of the matrix by the value specified.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            matrix->entry[row_itor][col_itor] * scalar;
        }
    }

    return ulapack_success;  
}

MatrixError_t ulapack_norm(const Matrix_t * const matrix, 
                           MatrixEntry_t * const norm) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The norm pointer should not be NULL because it should be written to.
     */
    if (!norm) {
        return ulapack_invalid_argument;
    }

    /*
     * Reset norm variable.
     */
    *norm = 0;

    /*
     * Scale the elements of the matrix by the value specified.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            *norm += matrix->entry[row_itor][col_itor] * 
                     matrix->entry[row_itor][col_itor];
        }
    }

    *norm = sqrt(*norm);

    return ulapack_success;
}

MatrixError_t ulapack_trace(const Matrix_t * const matrix, 
                            MatrixEntry_t * const trace) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The trace pointer should not be NULL because it should be written to.
     */
    if (!trace) {
        return ulapack_invalid_argument;
    }

    /*
     * Reset trace variable.
     */
    *trace = 0;

    /*
     * Scale the elements of the matrix by the value specified.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            *trace += matrix->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_dot(const Matrix_t * const vector_a,
                          const Matrix_t * const vector_b,
                          MatrixEntry_t * const dot) {
    /*
     * Vector element iterator.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t elem_itor = 0;

    /*
     * is_vector return value.
     */
    MatrixError_t isvect = ulapack_error;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(vector_a) || 
        !_ulapack_is_valid_memory(vector_b)) {
        return ulapack_uninit_obj;
    }

    /*
     * Check if input a is a vector. Only need to check one because the
     * dimensions must be equal in the next check.
     */
    ulapack_is_vector(vector_a, &isvect);
    if (!isvect) {
        return ulapack_invalid_argument;
    }

    /*
     * Check that dimensions match.
     */
    if ((vector_a->n_rows != vector_b->n_rows) || 
        (vector_a->n_cols != vector_b->n_cols)) {
        return ulapack_invalid_argument;
    }

    /*
     * Reset the result variable.
     */
    *dot = 0;

    /*
     * Check if row vector.
     */
    if (vector_a->n_cols > vector_a->n_rows) {
        /*
         * Perform the dot product operation on the columns if a row vector.
         */
        for (elem_itor = 0; elem_itor < vector_a->n_cols; elem_itor++) {
            *dot += vector_a->entry[0][elem_itor] * 
                    vector_b->entry[0][elem_itor];
        }
    } else {
        /*
         * Assume the vector is a column vector, and operate on the rows.
         */
        for (elem_itor = 0; elem_itor < vector_a->n_rows; elem_itor++) {
            *dot += vector_a->entry[elem_itor][0] * 
                    vector_b->entry[elem_itor][0];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_product(const Matrix_t * const A,
                              const Matrix_t * const B,
                              Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;
    uint64_t swap_itor = 0;

    /*
     * Verify that valid operand objects have been passed in.
     */
    if (!_ulapack_is_valid_memory(A) || !_ulapack_is_valid_memory(B)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_invalid_argument;
    }

    /*
     * Check for legal dimensions.
     */
    if (A->n_cols != B->n_rows) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = A->n_rows;
        result->n_cols = B->n_cols;

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the product of A and B.
     */
    #else
        if (A->n_rows != result->n_rows || 
            B->n_cols != result->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Compute and store product in the result matrix.
     */
    for(row_itor = 0; row_itor < A->n_rows; row_itor++) {
        for(col_itor = 0; col_itor < B->n_cols; col_itor++) {
            for(swap_itor = 0; swap_itor < A->n_cols; swap_itor++) {
                result->entry[row_itor][col_itor] += 
                A->entry[row_itor][swap_itor] * B->entry[swap_itor][col_itor];
            }
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_transpose(const Matrix_t * const matrix,
                                Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * transpose of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = matrix->n_cols;
        result->n_cols = matrix->n_rows;

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the transpose of the matrix.
     */
    #else
        if (matrix->n_rows != result->n_cols || 
            matrix->n_cols != result->n_rows) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Compute and store transpose in the result matrix.
     */
    for(row_itor = 0; row_itor < result->n_rows; row_itor++) {
        for(col_itor = 0; col_itor < result->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            matrix->entry[col_itor][row_itor];
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_det(const Matrix_t * const matrix, 
                          MatrixEntry_t * const det) {
    /*
     * Matrix element iterator.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t elem_itor = 0;

    /*
     * Make a copy of the input matrix to put it into upper triangular form.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t matrix_upper_tri;
        Matrix_t matrix_lower_tri;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *matrix_upper_tri = NULL;
        Matrix_t *matrix_lower_tri = NULL;
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * Verify if the matrix is square.
     */
    if (matrix->n_rows != matrix->n_cols) {
        return ulapack_invalid_argument;
    }

    /*
     * Reset determinant result.
     */
    *det = 0;

   /*
    * Switch on the square dimension of the matrix.
    */
    switch (matrix->n_rows) {
        /*
         * No matrix left.
         */
        case 0:
            break;
        /*
         * Scalar case. Determinant is itself.
         */
        case 1:
            *det = matrix->entry[0][0];
            break;
        /*
         * For the 2x2 case, just use the closed form solution for the matrix
         * determinant.
         */
        case 2:
            *det = (matrix->entry[0][0] * matrix->entry[1][1]) - 
                   (matrix->entry[1][0] * matrix->entry[0][1]);
            break;

        case 3:
            for(elem_itor = 0; elem_itor < 3; elem_itor++) {
                *det += (matrix->entry[0][elem_itor] * 
                    (matrix->entry[1][(elem_itor + 1) % 3] * 
                        matrix->entry[2][(elem_itor + 2) % 3] - 
                        matrix->entry[1][(elem_itor + 2) % 3] * 
                    matrix->entry[2][(elem_itor + 1) % 3]));
            }
            break;
        /*
         * For all cases where the dimensions of the matrix are 4x4 or more.
         */
        default:
            *det = 1;

            /*
             * Put the new matrix into upper triangular form.
             */
            #ifdef ULAPACK_USE_STATIC_ALLOC
                /*
                 * No need to compile in this code since the matrices won't
                 * exceed the closed form switch cases.
                 */
                #if ULAPACK_MAX_MATRIX_N_ROWS <= 3
                    return ulapack_invalid_argument;
                #else

                ulapack_init(&matrix_upper_tri, matrix->n_rows, matrix->n_cols);
                ulapack_init(&matrix_lower_tri, matrix->n_rows, matrix->n_cols);

                if (ulapack_lu(matrix, &matrix_upper_tri, &matrix_lower_tri) != 
                    ulapack_success) {
                    return ulapack_error;
                }

                for (elem_itor = 0; 
                    elem_itor < matrix_upper_tri.n_rows && 
                    elem_itor < matrix_upper_tri.n_cols; 
                    elem_itor++) {

                    *det *= 
                    isnormal(matrix_upper_tri.entry[elem_itor][elem_itor]) ? 
                            matrix_upper_tri.entry[elem_itor][elem_itor] : 0;
                }
                #endif
            #endif

            #ifdef ULAPACK_USE_DYNAMIC_ALLOC
                if (ulapack_init(&matrix_upper_tri, matrix->n_rows, 
                                                    matrix->n_cols) !=                  
                    ulapack_success) {

                    return ulapack_error;
                }

                if (ulapack_init(&matrix_lower_tri, 
                    matrix->n_rows, matrix->n_cols) != ulapack_success) {
                    return ulapack_error;
                }

                if (ulapack_lu(matrix, matrix_upper_tri, matrix_lower_tri) != 
                    ulapack_success) {
                    ulapack_destroy(matrix_upper_tri);
                    ulapack_destroy(matrix_lower_tri);

                    return ulapack_error;
                }

                for (elem_itor = 0; 
                    elem_itor < matrix_upper_tri->n_rows && 
                    elem_itor < matrix_upper_tri->n_cols; 
                    elem_itor++) {

                    *det *= 
                    isnormal(matrix_upper_tri->entry[elem_itor][elem_itor]) ? 
                            matrix_upper_tri->entry[elem_itor][elem_itor] : 0;                
                }

                /*
                 * Free the copy matrix memory.
                 */
                ulapack_destroy(matrix_upper_tri);
                ulapack_destroy(matrix_lower_tri);

            #endif

            break;
    }

    return ulapack_success;

}

MatrixError_t ulapack_lu(const Matrix_t * const matrix, 
                            Matrix_t * const upper_matrix,
                            Matrix_t * const lower_matrix) {

    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    uint64_t sub_row = 0;
    uint64_t sub_col = 0;

    /*
     * A copy of the input matrix.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t matrix_copy_mem;
        Matrix_t *matrix_copy = &matrix_copy_mem;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *matrix_copy = NULL;
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(upper_matrix)) {
        return ulapack_invalid_argument;
    }
    if (!_ulapack_is_valid_memory(lower_matrix)) {
        return ulapack_invalid_argument;
    }
    /*
     * The dimensions of the result pointer are changed to equal that of the
     * transpose of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        upper_matrix->n_rows = matrix->n_rows;
        upper_matrix->n_cols = matrix->n_cols;

        lower_matrix->n_rows = matrix->n_rows;
        lower_matrix->n_cols = matrix->n_cols;
    #endif
    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the transpose of the matrix.
     */
    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (matrix->n_rows != upper_matrix->n_rows || 
            matrix->n_cols != upper_matrix->n_cols ||
            matrix->n_rows != lower_matrix->n_rows ||
            matrix->n_cols != lower_matrix->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    #ifdef ULAPACK_USE_STATIC_ALLOC
        ulapack_init(matrix_copy, matrix->n_rows, matrix->n_cols);
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_init(&matrix_copy, matrix->n_rows, matrix->n_cols);
    #endif

    ulapack_copy(matrix, matrix_copy);

    /*
     * Using Cormen's algorithm from Introduction to Algorithms.
     */

    /*
     * Initialize memory to zero.
     */
    ulapack_set(upper_matrix, 0);
    ulapack_set(lower_matrix, 0);

    /*
     * Initialize lu components.
     */
    for (row_itor = 0; row_itor < upper_matrix->n_rows; row_itor++) {
        for(col_itor = 0; col_itor < upper_matrix->n_cols; col_itor++) {
            if (row_itor > col_itor) {
                upper_matrix->entry[row_itor][col_itor] = 0;
            } else if (row_itor < col_itor) {
                lower_matrix->entry[row_itor][col_itor] = 0;
            } else if (row_itor == col_itor) {
                lower_matrix->entry[row_itor][col_itor] = 1;
            }
        }
    }

    /*
     * Decompose.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        upper_matrix->entry[row_itor][row_itor] = 
        matrix_copy->entry[row_itor][row_itor];

        for (col_itor = row_itor + 1; col_itor < matrix->n_cols; col_itor++) {
            upper_matrix->entry[row_itor][col_itor] = 
            matrix_copy->entry[row_itor][col_itor];

            lower_matrix->entry[col_itor][row_itor] = 
                matrix_copy->entry[col_itor][row_itor] / 
                upper_matrix->entry[row_itor][row_itor];
        }
        for(sub_row = row_itor + 1; sub_row < matrix->n_rows; sub_row++){
            for(sub_col = row_itor + 1; sub_col < matrix->n_cols; sub_col++) {
                matrix_copy->entry[sub_row][sub_col] = 
                matrix_copy->entry[sub_row][sub_col] - 
                    (lower_matrix->entry[sub_row][row_itor] * 
                    upper_matrix->entry[row_itor][sub_col]);
            }
        }
    }

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_destroy(matrix_copy);
    #endif

    return ulapack_success;
}

/**
 * @name _ulapack_swap_rows
 * A private function to swap two rows in a matrix.
 *
 * @note The input matrix must be square.
 *
 * @param[in/out] matrix The matrix to manipulate.
 * @param row_a One of the operand rows.
 * @param row_b The other operand row.
 */
#ifndef ULAPACK_INVERSE_VIA_LU
static void _ulapack_swap_rows(Matrix_t * const matrix, 
                               const uint64_t row_a, 
                               const uint64_t row_b) {
    MatrixEntry_t swap_element = 0;

    uint64_t elem_itor = 0;

    for (elem_itor = 0; elem_itor < matrix->n_rows; elem_itor++) {
        swap_element = matrix->entry[row_a][elem_itor];
        matrix->entry[row_a][elem_itor] = matrix->entry[row_b][elem_itor];
        matrix->entry[row_b][elem_itor] = swap_element;
    }
}
#endif

/**
 * @name ulapack_backsub
 * Perform a backwards substitution between two matrices.
 *
 * @note Both matrices must be square and equal in size.
 *
 * @param[in/out] matrix_a The first matrix.
 * @param[in/out] matrix_b The second matrix.
 */
static void _ulapack_backsub(Matrix_t * const matrix_a, 
                             Matrix_t * const matrix_b) {
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;
    uint64_t sub_itor = 0;

    MatrixEntry_t pivot = 0;
    MatrixEntry_t ratio = 0;

    for(int row_itor = matrix_a->n_rows - 1; 
        (SIndex_t)row_itor >= 0; 
        row_itor--) {

        pivot = matrix_a->entry[row_itor][row_itor];
        for (col_itor = row_itor - 1; (SIndex_t)col_itor >= 0; col_itor--) {
            ratio = -1 * matrix_a->entry[col_itor][row_itor] / pivot;

            for(sub_itor = 0; sub_itor < matrix_a->n_rows; sub_itor++) {

                matrix_a->entry[col_itor][sub_itor] += ratio * 
                matrix_a->entry[row_itor][sub_itor];

                matrix_b->entry[col_itor][sub_itor] += ratio * 
                matrix_b->entry[row_itor][sub_itor];
            }
        }
    }

    for (row_itor = 0; row_itor < matrix_a->n_rows; row_itor++) {
        
        for (col_itor = 0; col_itor < matrix_a->n_rows; col_itor++) {
            matrix_b->entry[row_itor][col_itor] /= 
            matrix_a->entry[row_itor][row_itor];
        }

        matrix_a->entry[row_itor][row_itor] = 1;
    }
}

/**
 * @name ulapack_forwardsub
 * Perform a forward substitution between two matrices.
 *
 * @note Both matrices must be square and equal in size.
 *
 * @param[in/out] matrix_a The first matrix (must be upper triangular).
 * @param[in/out] matrix_b The second matrix.
 */
#ifdef ULAPACK_INVERSE_VIA_LU
static void _ulapack_forwardsub(Matrix_t * const matrix_a, 
                                Matrix_t * const matrix_b) {
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;
    uint64_t sub_itor = 0;

    MatrixEntry_t pivot = 0;
    MatrixEntry_t ratio = 0;

    for (row_itor = 0; row_itor < matrix_a->n_rows; row_itor++) {
        pivot = matrix_a->entry[row_itor][row_itor];

        for (col_itor = row_itor + 1; col_itor < matrix_a->n_rows; col_itor++) {
            ratio = -1 * matrix_a->entry[col_itor][row_itor]/pivot;

            for (sub_itor = 0; sub_itor < matrix_a->n_rows; sub_itor++) {
                matrix_a->entry[col_itor][sub_itor] += 
                    ratio * matrix_a->entry[row_itor][sub_itor];

                matrix_b->entry[col_itor][sub_itor] += 
                    ratio * matrix_b->entry[row_itor][sub_itor];
            }
        }
    }
}
#endif

MatrixError_t ulapack_inverse(const Matrix_t * const matrix,
                              Matrix_t * const inverse) {

    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    #ifndef ULAPACK_INVERSE_VIA_LU
        uint64_t row_itor = 0;
        uint64_t col_itor = 0;

        uint64_t pivot_index = 0;

        MatrixEntry_t pivot_element = 0;
        MatrixEntry_t ratio = 0;
    #endif

    /*
     * A copy of the input matrix.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t matrix_copy_mem;
        Matrix_t *matrix_copy = &matrix_copy_mem;

        #ifdef ULAPACK_INVERSE_VIA_LU
            Matrix_t upper_matrix_mem;
            Matrix_t lower_matrix_mem;

            Matrix_t *upper_matrix = &upper_matrix_mem;
            Matrix_t *lower_matrix = &lower_matrix_mem;
        #endif   
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *matrix_copy = NULL;
        
        #ifdef ULAPACK_INVERSE_VIA_LU
            Matrix_t *upper_matrix = NULL;
            Matrix_t *lower_matrix = NULL;
        #endif
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(inverse)) {
        return ulapack_uninit_obj;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * transpose of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        inverse->n_rows = matrix->n_rows;
        inverse->n_cols = matrix->n_cols;
    #endif

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the transpose of the matrix.
     */
    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (matrix->n_rows != inverse->n_rows || 
            matrix->n_cols != inverse->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    #ifdef ULAPACK_USE_STATIC_ALLOC
        ulapack_init(matrix_copy, matrix->n_rows, matrix->n_cols);

        #ifdef ULAPACK_INVERSE_VIA_LU
            ulapack_init(upper_matrix, matrix->n_rows, matrix->n_cols);
            ulapack_init(lower_matrix, matrix->n_rows, matrix->n_cols);
        #endif
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_init(&matrix_copy, matrix->n_rows, matrix->n_cols);

        #ifdef ULAPACK_INVERSE_VIA_LU
            ulapack_init(&upper_matrix, matrix->n_rows, matrix->n_cols);
            ulapack_init(&lower_matrix, matrix->n_rows, matrix->n_cols);
        #endif
    #endif

    ulapack_copy(matrix, matrix_copy);

    /*
     * Start by setting the inverse result to the identity matrix.
     */
    ulapack_eye(inverse);

    #ifdef ULAPACK_INVERSE_VIA_LU

        ulapack_lu(matrix, upper_matrix, lower_matrix);
        _ulapack_forwardsub(lower_matrix, inverse);
        _ulapack_backsub(upper_matrix, inverse);

        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(upper_matrix);
            ulapack_destroy(lower_matrix);
        #endif

    /*
     * Method not using LU decomposition.
     */
    #else
        for (row_itor = 0; row_itor < matrix_copy->n_rows; row_itor++) {
            
            pivot_index = row_itor;
            while (pivot_element < MINIMUM_THRESHOLD_TOLERANCE && 
                pivot_index < matrix_copy->n_rows) {

                pivot_element = matrix_copy->entry[pivot_index][row_itor];
                pivot_index++;
            }

            if (pivot_element < MINIMUM_THRESHOLD_TOLERANCE) {
                /*
                 * The matrix is singular.
                 */
                return ulapack_error;
            }

            if (row_itor != pivot_index - 1) {
                _ulapack_swap_rows(matrix_copy, row_itor, pivot_index);
                _ulapack_swap_rows(inverse, row_itor, pivot_index);
            }

            for (col_itor = row_itor + 1; 
                col_itor < matrix_copy->n_rows; 
                col_itor++) {

                ratio = -1 * matrix_copy->entry[col_itor][row_itor] / 
                             pivot_element;

                for(pivot_index = 0; 
                    pivot_index < matrix_copy->n_rows; 
                    pivot_index++) {
                    
                    matrix_copy->entry[col_itor][pivot_index] += 
                    ratio * matrix_copy->entry[row_itor][pivot_index];

                    inverse->entry[col_itor][pivot_index] += 
                    ratio * inverse->entry[row_itor][pivot_index];
                }   
            }
        }

        _ulapack_backsub(matrix_copy, inverse);
    
    #endif // ULAPACK_INVERSE_VIA_LU

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_destroy(matrix_copy);
    #endif

    return ulapack_success;
}

MatrixError_t ulapack_pinverse(const Matrix_t * const matrix,
                               Matrix_t * const pinverse) {
    MatrixError_t ret_code;

    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t matrix_transpose_mem;
        Matrix_t matrix_transpose_term_product_mem;
        Matrix_t matrix_transpose_term_inverse_mem;

        Matrix_t *matrix_transpose = &matrix_transpose_mem;

        Matrix_t *matrix_transpose_term_product = 
        &matrix_transpose_term_product_mem;
        
        Matrix_t *matrix_transpose_term_inverse = 
        &matrix_transpose_term_inverse_mem;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *matrix_transpose = NULL;
        Matrix_t *matrix_transpose_term_product = NULL;
        Matrix_t *matrix_transpose_term_inverse = NULL;
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(pinverse)) {
        return ulapack_uninit_obj;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * pseudo inverse of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        pinverse->n_rows = matrix->n_cols;
        pinverse->n_cols = matrix->n_rows;
    #endif
        
    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the pseudo inverse of the matrix.
     */
    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (matrix->n_rows != pinverse->n_cols || 
            matrix->n_cols != pinverse->n_rows) {
            return ulapack_invalid_argument;
        }
    #endif

    #ifdef ULAPACK_USE_STATIC_ALLOC
        /*
         * Memory for matrix^T.
         */
        if (ulapack_init(matrix_transpose, 
            matrix->n_cols, matrix->n_rows) != ulapack_success) {
            return ulapack_error;
        }

        /*
         * Memory for (matrix^T * matrix).
         */
        if (ulapack_init(matrix_transpose_term_product, 
            matrix->n_cols, matrix->n_cols) != ulapack_success) {
            return ulapack_error;
        }

        /*
         * Memory for (matrix^T * matrix)^-1.
         */
        if (ulapack_init(matrix_transpose_term_inverse, 
            matrix->n_cols, matrix->n_cols) != ulapack_success) {
            return ulapack_error;
        }
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        /*
         * Memory for matrix^T.
         */
        ret_code = ulapack_init(&matrix_transpose, 
            matrix->n_cols, matrix->n_rows);

        if (ret_code != ulapack_success) {
            return ret_code;
        }

        /*
         * Memory for (matrix^T * matrix).
         */
        ret_code = ulapack_init(&matrix_transpose_term_product, 
            matrix->n_cols, matrix->n_cols);

        if (ret_code != ulapack_success) {
            ulapack_destroy(matrix_transpose);
            return ret_code;
        }

        /*
         * Memory for (matrix^T * matrix)^-1.
         */
        ret_code = ulapack_init(&matrix_transpose_term_inverse, 
            matrix->n_cols, matrix->n_cols);

        if (ret_code != ulapack_success) {
            ulapack_destroy(matrix_transpose);
            ulapack_destroy(matrix_transpose_term_product);
            return ret_code;
        }
    #endif

    ulapack_set(pinverse, 0);
    ulapack_set(matrix_transpose, 0);
    ulapack_set(matrix_transpose_term_product, 0);
    ulapack_set(matrix_transpose_term_inverse, 0);

    /*
     * A^T.
     */
    ret_code = ulapack_transpose(matrix, matrix_transpose);

    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(matrix_transpose);
            ulapack_destroy(matrix_transpose_term_product);
            ulapack_destroy(matrix_transpose_term_inverse);
        #endif

        return ret_code;
    }

    /*
     * (A^T * A).
     */
    ret_code = ulapack_product(matrix_transpose, matrix, 
        matrix_transpose_term_product);

    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(matrix_transpose);
            ulapack_destroy(matrix_transpose_term_product);
            ulapack_destroy(matrix_transpose_term_inverse);
        #endif

        return ret_code;
    }

    /*
     * (A^T * A)^-1.
     */
    ret_code = ulapack_inverse(matrix_transpose_term_product, 
        matrix_transpose_term_inverse);

    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(matrix_transpose);
            ulapack_destroy(matrix_transpose_term_product);
            ulapack_destroy(matrix_transpose_term_inverse);
        #endif

        return ret_code;
    }

    /*
     * A^* = (A^T * A)^-1 * A^T.
     */
    ret_code = ulapack_product(matrix_transpose_term_inverse, matrix_transpose, 
        pinverse);

    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(matrix_transpose);
            ulapack_destroy(matrix_transpose_term_product);
            ulapack_destroy(matrix_transpose_term_inverse);
        #endif

        return ret_code;
    }


    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_destroy(matrix_transpose);
        ulapack_destroy(matrix_transpose_term_product);
        ulapack_destroy(matrix_transpose_term_inverse);
    #endif

    return ulapack_success;
}

MatrixError_t ulapack_least_squares(const Matrix_t * const A,
                                    const Matrix_t * const y,
                                    Matrix_t * const x) {

    MatrixError_t ret_code;

    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t pinv_mem;
        Matrix_t *pinv = &pinv_mem;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *pinv = NULL;
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(A) || !_ulapack_is_valid_memory(y)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(x)) {
        return ulapack_uninit_obj;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * pseudo inverse of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        x->n_rows = A->n_cols;
        x->n_cols = 1;
    #endif
        
    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the expected dimensions.
     */
    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (x->n_rows != A->n_cols || 
            x->n_cols != 1) {
            return ulapack_invalid_argument;
        }
    #endif

    #ifdef ULAPACK_USE_STATIC_ALLOC
        /*
         * Memory for A^*.
         */
        if (ulapack_init(pinv, A->n_cols, A->n_rows) != ulapack_success) {
            return ulapack_error;
        }
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        /*
         * Memory for A^*.
         */
        ret_code = ulapack_init(&pinv, A->n_cols, A->n_rows);

        if (ret_code != ulapack_success) {
            return ret_code;
        }
    #endif

    /*
     * A^*.
     */
    ret_code = ulapack_pinverse(A, pinv);

    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(pinv);
        #endif

        return ret_code;
    }

    /*
     * x ~= (A^* * y).
     */
    ret_code = ulapack_product(pinv, y, x);
    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(pinv);
        #endif

        return ret_code;
    }


    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_destroy(pinv);
    #endif

    return ulapack_success;
}

/**
 * @name ulapack_elem_pow
 * Take the power of a matrix element.
 *
 * @param element The element value to take to the power of p.
 * @param p The power to take the element to.
 *
 * @return element^p.
 */
static MatrixEntry_t ulapack_elem_pow(const MatrixEntry_t elem, 
                                      const MatrixEntry_t p) {
    uint64_t pitor = 0;

    MatrixEntry_t power = 0;

    if (p == 0) {
        return 1;
    } else {

        if (elem == 0) {
            return 0;
        }

        power = 1;
        for (pitor = 0; pitor < p; pitor++) {
            power *= elem;
        }        
    }

    return power;
}

MatrixError_t ulapack_vandermonde(const Matrix_t * const x,
                                  const Index_t order, 
                                  Matrix_t * const V) {
    uint64_t col_itor = 0;
    uint64_t v_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(x)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(V)) {
        return ulapack_uninit_obj;
    }

    if (order > x->n_rows) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * Vandermonde of the operands if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        V->n_rows = x->n_rows;
        V->n_cols = order;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (V->n_rows != x->n_rows || 
            V->n_cols != order) {
            return ulapack_invalid_argument;
        }
    #endif

    for (v_itor = 0; v_itor < x->n_rows; v_itor++) {
        for (col_itor = 0; col_itor < order; col_itor++) {
            V->entry[v_itor][order - 1 - col_itor] = 
            ulapack_elem_pow(x->entry[v_itor][0], col_itor);
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_power(const Matrix_t * const matrix, 
                            const uint64_t power,
                            Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }
    
    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_invalid_argument;
    }

    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = matrix->n_rows;
        result->n_cols = matrix->n_cols;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (result->n_rows != matrix->n_rows || 
            result->n_cols != matrix->n_cols) {
            return ulapack_invalid_argument;
        }
    #endif

    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            ulapack_elem_pow(matrix->entry[row_itor][col_itor], power);
        }
    }

    return ulapack_success;
}

MatrixError_t ulapack_polyfit(const Matrix_t * const x, 
                              const Matrix_t * const y,
                              const uint64_t n,
                              Matrix_t * const p) {
    MatrixError_t ret_code;

    #ifdef ULAPACK_USE_STATIC_ALLOC
        Matrix_t V_mem;

        Matrix_t *V = &V_mem;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        Matrix_t *V = NULL;
    #endif

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(x) || !_ulapack_is_valid_memory(y)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(p)) {
        return ulapack_uninit_obj;
    }

    /*
     * x and y must have equal dimensions.
     */
    if (x->n_rows != y->n_rows || x->n_cols != y->n_cols) {
        return ulapack_invalid_argument;
    }

    #ifdef ULAPACK_USE_STATIC_ALLOC
        p->n_rows = n + 1;
        p->n_cols = 1;
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if (p->n_rows != n + 1 || p->n_cols != 1) {
            return ulapack_invalid_argument;
        }
    #endif

    #ifdef ULAPACK_USE_STATIC_ALLOC
        ret_code = ulapack_init(V, x->n_rows, n + 1);
        if (ret_code != ulapack_success) {
            return ret_code;
        }
    #endif

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ret_code = ulapack_init(&V, x->n_rows, n + 1);
        if (ret_code != ulapack_success) {
            return ret_code;
        }

    #endif

    /*
     * Compute the Vandermonde matrix.
     */
    ret_code = ulapack_vandermonde(x, n + 1, V);
    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(V);
        #endif

        return ret_code;
    }

    ret_code = ulapack_least_squares(V, y, p);
    if (ret_code != ulapack_success) {
        #ifdef ULAPACK_USE_DYNAMIC_ALLOC
            ulapack_destroy(V);
        #endif

        return ret_code;
    }

    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        ulapack_destroy(V);
    #endif
            
    return ulapack_success;
}

MatrixError_t ulapack_copy(const Matrix_t * const matrix,
                           Matrix_t * const result) {
    /*
     * Matrix iterators.
     * Declared at top of function for legacy compiler comparability.
     */
    uint64_t row_itor = 0;
    uint64_t col_itor = 0;

    /*
     * Verify that a valid operand object have been passed in.
     */
    if (!_ulapack_is_valid_memory(matrix)) {
        return ulapack_uninit_obj;
    }

    /*
     * The result pointer should not be NULL because it should be written to.
     */
    if (!_ulapack_is_valid_memory(result)) {
        return ulapack_invalid_argument;
    }

    /*
     * The dimensions of the result pointer are changed to equal that of the
     * input matrix if static allocation is used.
     */
    #ifdef ULAPACK_USE_STATIC_ALLOC
        result->n_rows = matrix->n_rows;
        result->n_cols = matrix->n_cols;
    #endif

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the input matrix.
     */
    #ifdef ULAPACK_USE_DYNAMIC_ALLOC
        if ((matrix->n_rows != result->n_rows) || 
            (matrix->n_cols != result->n_cols)) {
            return ulapack_invalid_argument;
        }
    #endif

    /*
     * Copy the matrix elements.
     */
    for(row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for(col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            matrix->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;    
}
