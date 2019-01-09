/**
 * @file ulapack.h
 * @brief Declarations of uLAPack matrix manipulation functions.
 *
 * @author Sargis Yonan
 * @date July 8, 2017
 *
 * @version 1.0.1
 **/

/**
 * @todo new line the return codes in dox for readability.
 */

/*
 * File header guards.
 */
#ifndef __ULAPACK_H__
#define __ULAPACK_H__

/*
 * Including standard integer header for uintX_t types.
 */
#include <stdint.h>

/*
 * uLAPack options header.
 */
#include "ulapack_options.h"

/*
 * uLAPack type definitions.
 */
#include "ulapack_type.h"

/**
 * @name ulapack_initialize_matrix
 *
 * Initializes uLAPack matrix object.
 *
 * @note If static memory allocation is used, then only the n_rows and n_cols
 *       members of the matrix structure are set. Otherwise, memory is allocated
 *       using the specified memory allocator call.
 *
 * @param[in/out] matrix A pointer to a matrix object type.
 * @param n_rows The number of rows the matrix will have.
 * @param n_cols The number of columns the matrix will have.
 *
 * @return ulapack_success is returned upon successful modification and
 *         initialization of the matrix object.
 *
 *         ulapack_oom is returned if dynamic memory allocation is specified for 
 *         use, and the memory required for the object cannot be allocated. 
 *
 *         ulapack_invalid_argument is returned if a NULL pointer was passed in
 *         while using static memory allocation, or a non-NULL pointer was 
 *         passed in using dynamic memory allocation. 
 *
 *         ulapack_invalid_argument is returned when the rows and 
 *         columns passed in exceed the preset maximum values set in the uLAPack
 *         options.
 */
ulapack_error_t ulapack_initialize_matrix(
    #ifdef ULAPACK_USE_STATIC_ALLOC
        ulapack_matrix_t *matrix,
    #else
        ulapack_matrix_t **matrix, 
    #endif
    const uint64_t n_rows, const uint64_t n_cols);

/**
 * @name ulapack_destroy
 * Free a dynamically allocated ULAPack matrix object, and reset its pointers.
 *
 * @note This function can only be called in dynamic memory allocation is
 *       specified by explicitly not defining static allocation.
 * @note The free function defined by the macro ULAPACK_FREE is used as the
 *       memory freeing call.
 *
 * @param[in] obj A pointer to the ULAPack matrix object to destroy.
 *
 * @return ULAPack success code ulapack_success is returned if the matrix object
 *         passed in was valid and freed successfully.
 *         ulapack_uninit_obj is returned if the matrix object is not
 *         initialized.
 */
#ifndef ULAPACK_USE_STATIC_ALLOC
ulapack_error_t ulapack_destroy(ulapack_matrix_t *obj);
#endif

/**
 * @name ulapack_edit_entry
 * Edit a single uLAPack matrix object's element entry.
 *
 * @param[in/out] matrix The matrix object to edit an entry within.
 * @param row The row to edit the matrix within.
 * @param col The column to edit the matrix within.
 * @param value The value to overwrite the matrix element at the coordinates
 *        (row, col).
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         entry edit.
 *         ulapack_uninit_obj is returned if the matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the coordinates to edit are
 *         out of the range of the initialized matrix dimensions.
 */
ulapack_error_t ulapack_edit_entry(ulapack_matrix_t * const matrix, 
                                   const uint64_t row, const uint64_t col,
                                   const uint64_t value);

/**
 * @name ulapack_get_entry
 * Get a single uLAPack matrix object's element entry.
 *
 * @param[in] matrix The matrix object to get an entry from.
 * @param row The row to edit the matrix within.
 * @param col The column to edit the matrix within.
 * @param[out] value The variable to copy the matrix value at coordinates
 *             (row,col) into. ulapack_success is passed back if the matrices 
 *             are equal. ulapack_error is returned if the matrices are not 
 *             equal.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         copy of the matrix element at coordinates (row, col).
 *         ulapack_uninit_obj is returned if the matrix object is not
 *         initialized.
 *         ulapack_invalid_argument is returned if value is NULL.
 *         ulapack_invalid_argument is returned if the coordinates to get from
 *         are out of the range of the initialized matrix dimensions.
 */
ulapack_error_t ulapack_get_entry(const ulapack_matrix_t * const matrix, 
                                  const uint64_t row, const uint64_t col,
                                  ulapack_entry_t * const value);

/**
 * @name ulapack_size
 * Get the dimensions of a matrix object.
 *
 * @param[in] matrix The matrix object to get dimensions of.
 * @param[out] rows The number of rows in the matrix object.
 * @param[out] cols The number of columns in the matrix object.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         copy of matrix dimensions into the references rows and cols.
 *         ulapack_uninit_obj is returned if the matrix object is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if rows or cols is NULL.
 */
ulapack_error_t ulapack_size(const ulapack_matrix_t * const matrix, 
                             uint64_t * const rows, uint64_t * const cols);

/**
 * @name ulapack_set
 * Set every entry in a matrix to a specified value.
 *
 * @param[in/out] matrix An initialized matrix object to modify.
 * @param value The value to set the matrix elements to.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         modification of the matrix object.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 */
ulapack_error_t ulapack_set(ulapack_matrix_t * const matrix,
                            const ulapack_entry_t value);

/**
 * @name ulapack_is_equal
 * Check if the two matrix operands passed in are equal to each other (A == B).
 *
 * @param[in] A An initialized matrix object operand.
 * @param[in] B An initialized matrix object operand.
 * @param[out] is_equal The conditional that stores the equality (A == B). 
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         comparison between the two matrices.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the is_equal pointer is NULL.
 */
ulapack_error_t ulapack_is_equal(const ulapack_matrix_t * const A, 
                                 const ulapack_matrix_t * const B,
                                 ulapack_entry_t * const is_equal);

/**
 * @name ulapack_eye
 * Set a matrix equal to the identity matrix of that dimension.
 *
 * @note The matrix must be square (n_rows == n_cols).
 *
 * @param[in/out] matrix An initialized matrix object.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         modification of the matrix object.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the matrix is not square.
 */
ulapack_error_t ulapack_eye(ulapack_matrix_t * const matrix);

/**
 * @name ulapack_sum
 * Add the elements of a matrix or vector and pass the value back by reference.
 *
 * @param[in] matrix An initialized matrix object.
 * @param[out] result The resultant of the sum of the elements of the matrix 
 *             passed in.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         summation of the matrix object.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 */
ulapack_error_t ulapack_sum(const ulapack_matrix_t * const matrix,
                            ulapack_entry_t * const result);

/**
 * @name ulapack_add
 * Add two matrices or vectors together.
 *
 * @note the dimensions of the two matrices must be identical.
 * @note if static allocation is specified, the dimensions of the result matrix
 *       are modified to equal the dimensions of the operand matrices.
 *
 * @param[in] A An initialized operand matrix.
 * @param[in] B An initialized operand matrix.
 * @param[out] An initialized matrix object to store the result of A + B.
 *
 * @return Upon a successful summation ulapack_success is returned.
 *         ulapack_invalid_argument is returned if size(A) != size(B).
 *         ulapack_uninit_obj is returned if the operand objects passed in are 
 *         not initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 *         ulapack_invalid_argument is returned if dynamic allocation is 
 *         specified, and the dimensions of the result matrix do not equal that
 *         of the operands.
 */
ulapack_error_t ulapack_add(const ulapack_matrix_t * const A, 
                            const ulapack_matrix_t * const B,
                            ulapack_matrix_t * const result);

/**
 * @name ulapack_subtract
 * Take the difference between two matrices or vectors.
 *
 * @note the dimensions of the two matrices must be identical.
 * @note if static allocation is specified, the dimensions of the result matrix
 *       are modified to equal the dimensions of the operand matrices.
 *
 * @param[in] A An initialized operand matrix.
 * @param[in] B An initialized operand matrix.
 * @param[out] result An initialized matrix object to store the result of A - B.
 *
 * @return Upon a successful subtraction ulapack_success is returned.
 *         ulapack_invalid_argument is returned if size(A) != size(B).
 *         ulapack_uninit_obj is returned if the operand objects passed in are 
 *         not initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 *         ulapack_invalid_argument is returned if dynamic allocation is 
 *         specified, and the dimensions of the result matrix do not equal that
 *         of the operands.
 */
ulapack_error_t ulapack_subtract(const ulapack_matrix_t * const A, 
                                 const ulapack_matrix_t * const B,
                                 ulapack_matrix_t * const result);

/**
 * @name ulapack_scale
 * Multiply a matrix/vector by a scalar.
 *
 * @param[in] matrix An initialized operand matrix object to scale.
 * @param scalar A value to scale a matrix by.
 * @param[out] result The matrix to store the result of the scaling into.
 *
 * @note If static memory is specified, the dimensions of the result matrix are
 *       set to that of the operand matrix.
* @note If dynamic memory allocation is specified, the initialized dimensions 
 *       of the result matrix must be equal to that of the input matrix.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         operation.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if dynamic allocation is 
 *         specified, and the dimensions of the result matrix do not equal that
 *         of the operands.
 */
ulapack_error_t ulapack_scale(ulapack_matrix_t * const matrix, 
                              const ulapack_entry_t scalar,
                              ulapack_matrix_t * const result);

/**
 * @name ulapack_norm
 * Take the Frobenius norm of a matrix or norm of a vector.
 *
 * @param[in] matrix An initialized matrix object to take the norm of.
 * @param[out] norm A value to scale a matrix by.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         operation.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the norm pointer is NULL.
 */
ulapack_error_t ulapack_norm(const ulapack_matrix_t * const matrix,
                             ulapack_entry_t * const norm);

/**
 * @name ulapack_product
 * Multiply two matrices/vectors.
 *
 * @param[in] A An initialized matrix object operand.
 * @param[in] B An initialized matrix object operand.
 * @param[out] result An initialized matrix object to store the product of A and
 *             B. result = AB.
 *
 * @note For a matrix, A, with dimensions NxM, and a matrix B, with dimensions
 *       MxK, the product of the two matrices, AB, will have dimensions NxK.
 * @note If dynamic memory allocation is specified, the initialized dimensions 
 *       of the result matrix must be NxK.
 * @note If static memory allocation is specified, the dimensions of the result
 *       matrix are set to NxK.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         operation.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 *         ulapack_invalid_argument is returned if dynamic memory allocation is
 *         specified and the dimensions of the result matrix do not equal the
 *         required dimensions of the matrix/vector multiplication.
 */
ulapack_error_t ulapack_product(const ulapack_matrix_t * const A,
                                const ulapack_matrix_t * const B,
                                ulapack_matrix_t * const result);


/**
 * @name ulapack_transpose
 * Take the transpose of a matrix.
 *
 * @param[in] matrix An initialized matrix object operand.
 * @param[out] result The transpose of the input matrix.
 *
 * @note For a matrix, with dimensions NxM, the transpose result matrix will
 *       have dimensions MxN.
 * @note If static memory allocation is specified, the dimensions of the result
 *       matrix are set to NxK.
 * @note If dynamic memory allocation is used, the dimensions of the result
 *       matrix must be MxN.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         operation.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 *         ulapack_invalid_argument is returned if dynamic memory allocation is
 *         specified and the dimensions of the result matrix do not equal the
 *         required dimensions of the matrix/vector transpose.
 */
ulapack_error_t ulapack_transpose(const ulapack_matrix_t * const matrix,
                                  ulapack_matrix_t * const result);

/**
 * @name ulapack_copy
 * Copy a matrix.
 *
 * @param[in] matrix An initialized matrix object operand.
 * @param[out] result The copy of the input matrix.
 *
 * @note For a matrix, with dimensions NxM, the copy result matrix will
 *       have dimensions NxM.
 * @note If static memory allocation is specified, the dimensions of the result
 *       matrix are set to NxM.
 * @note If dynamic memory allocation is used, the dimensions of the result
 *       matrix must be NxM.
 *
 * @return ULAPack success code ulapack_success is returned upon a successful
 *         operation.
 *         ulapack_uninit_obj is returned if a matrix object passed in is not 
 *         initialized.
 *         ulapack_invalid_argument is returned if the result pointer is NULL.
 *         ulapack_invalid_argument is returned if dynamic memory allocation is
 *         specified and the dimensions of the result matrix do not equal the
 *         required dimensions of the matrix/vector copy.
 */
ulapack_error_t ulapack_copy(const ulapack_matrix_t * const matrix,
                             ulapack_matrix_t * const result);

/*
 * End header guard definition.
 */
#endif
