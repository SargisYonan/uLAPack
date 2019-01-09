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
static ulapack_error_t _ulapack_is_valid_memory(const ulapack_matrix_t *mat) {
    if (!mat) {
        return ulapack_error;
    }

    return ulapack_success;
}

ulapack_error_t ulapack_initialize_matrix(
    #ifdef ULAPACK_USE_STATIC_ALLOC
        ulapack_matrix_t *matrix,
    #else
        ulapack_matrix_t **matrix, 
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
            return ulapack_invalid_argument;
        }

        /*
         * Set the matrix dimensions.
         */
        matrix->n_rows = n_rows;
        matrix->n_cols = n_cols;

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
        *matrix = ULAPACK_ALLOC(sizeof(ulapack_matrix_t));
        if (!matrix) {
            return ulapack_oom;
        }
        /*
         * Allocate memory for the matrix rows.
         */
        (*matrix)->entry = ULAPACK_ALLOC(sizeof(ulapack_entry_t) * n_rows);
        
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
                sizeof(ulapack_entry_t) * n_cols);

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

#ifndef ULAPACK_USE_STATIC_ALLOC
ulapack_error_t ulapack_destroy(ulapack_matrix_t *obj) {
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

ulapack_error_t ulapack_edit_entry(ulapack_matrix_t * const matrix, 
                                   const uint64_t row, const uint64_t col,
                                   const uint64_t value) {
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

ulapack_error_t ulapack_get_entry(const ulapack_matrix_t * const matrix, 
                                  const uint64_t row, const uint64_t col,
                                  ulapack_entry_t * const value) {
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

ulapack_error_t ulapack_size(const ulapack_matrix_t * const matrix, 
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

ulapack_error_t ulapack_set(ulapack_matrix_t * const matrix,
                            const ulapack_entry_t value) {
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

ulapack_error_t ulapack_is_equal(const ulapack_matrix_t * const A, 
                                 const ulapack_matrix_t * const B,
                                 ulapack_entry_t * const is_equal) {
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

ulapack_error_t ulapack_eye(ulapack_matrix_t * const matrix) {
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

ulapack_error_t ulapack_sum(const ulapack_matrix_t * const matrix,
                            ulapack_entry_t * const result) {
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
     * Add the elements of the matrix.
     */
    for (row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for (col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            *result += matrix->entry[row_itor][col_itor];
        }
    }

    return ulapack_success;
}

ulapack_error_t ulapack_add(const ulapack_matrix_t * const A, 
                            const ulapack_matrix_t * const B,
                            ulapack_matrix_t * const result) {
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

ulapack_error_t ulapack_subtract(const ulapack_matrix_t * const A, 
                                 const ulapack_matrix_t * const B,
                                 ulapack_matrix_t * const result) {
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

ulapack_error_t ulapack_scale(ulapack_matrix_t * const matrix, 
                              const ulapack_entry_t scalar,
                              ulapack_matrix_t * const result) {
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

ulapack_error_t ulapack_norm(const ulapack_matrix_t * const matrix, 
                             ulapack_entry_t * const norm) {
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

ulapack_error_t ulapack_product(const ulapack_matrix_t * const A,
                                const ulapack_matrix_t * const B,
                                ulapack_matrix_t * const result) {
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

ulapack_error_t ulapack_transpose(const ulapack_matrix_t * const matrix,
                                  ulapack_matrix_t * const result) {
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
    for(row_itor = 0; row_itor < matrix->n_rows; row_itor++) {
        for(col_itor = 0; col_itor < matrix->n_cols; col_itor++) {
            result->entry[row_itor][col_itor] = 
            matrix->entry[col_itor][row_itor];
        }
    }

    return ulapack_success;
}

ulapack_error_t ulapack_copy(const ulapack_matrix_t * const matrix,
                                  ulapack_matrix_t * const result) {
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

    /*
     * If dynamic memory allocation is specified, the dimensions of result 
     * should match that of the dimensions of the input matrix.
     */
    #else
        if (matrix->n_rows != result->n_rows || 
            matrix->n_cols != result->n_cols) {
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