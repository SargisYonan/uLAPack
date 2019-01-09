/**
 * @file main.c
 * @brief Main test driving code for uLAPack.
 *
 * @author Sargis Yonan
 * @date July 8, 2017
 *
 * @version 1.0.1
 **/

/*
 * printf()
 */
#include <stdio.h>

#include "ulapack.h"

#include <stdbool.h>

static bool ut_iserror(ulapack_error_t err) {
    if (err != ulapack_success) {
        printf("Unit test error: return code %d\n", (int)err);
        return true;
    }

    return false;
}

static void print_matrix(const ulapack_matrix_t *mat) {
    for (uint64_t row_itor = 0; row_itor < mat->n_rows; row_itor++) {
        for (uint64_t col_itor = 0; col_itor < mat->n_cols; col_itor++) {
            printf("%lf ", mat->entry[row_itor][col_itor]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
}

static ulapack_error_t test_initialization(void) {
    ulapack_error_t ret_code;

    #ifdef ULAPACK_USE_STATIC_ALLOC

        printf("Testing static memory allocation.\n");
        ulapack_matrix_t A;
        ulapack_matrix_t B;
        ulapack_matrix_t result;

        /*
         * Test static initialization.
         */
        ret_code = ulapack_initialize_matrix(&A, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix A.\n");
            return ret_code;
        }
        if (A.n_rows != 3 || A.n_cols != 3) {
            printf("Matrix dimensions do not match expected values.\n");
            return ulapack_error;
        }

        ret_code = ulapack_initialize_matrix(&B, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix B.\n");
            return ret_code;
        }
        if (B.n_rows != 3 || B.n_cols != 3) {
            printf("Matrix dimensions do not match expected values.\n");
            return ulapack_error;
        }

        ret_code = ulapack_initialize_matrix(&result, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix result.\n");
            return ret_code;
        }
        if (result.n_rows != 3 || result.n_cols != 3) {
            printf("Matrix dimensions do not match expected values.\n");
            return ulapack_error;
        }

    #else

        printf("Testing dynamic memory allocation.\n");
        ulapack_matrix_t *A = NULL;
        ulapack_matrix_t *B = NULL;
        ulapack_matrix_t *result = NULL;

        /*
         * Test dynamic initialization.
         */
        ret_code = ulapack_initialize_matrix(&A, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix A.\n");
            return ret_code;
        }
        if (A->n_rows != 3 || A->n_cols != 3) {
            printf("Matrix dimensions do not match expected values.\n");
            return ulapack_error;
        }

        ret_code = ulapack_initialize_matrix(&B, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix B.\n");
            return ret_code;
        }
        if (B->n_rows != 3 || B->n_cols != 3) {
             printf("Matrix dimensions do not match expected values.\n");
             return ulapack_error;
        }

        ret_code = ulapack_initialize_matrix(&result, 3u, 3u);
        if (ut_iserror(ret_code)) {
            printf("Cannot initialize matrix result.\n");
            return ret_code;
        }
        if (result->n_rows != 3 || result->n_cols != 3) {
            printf("Matrix dimensions do not match expected values.\n");
            return ulapack_error;
        }

    #endif

    return ulapack_success;
}

static ulapack_error_t test_addition(void) {
    #ifdef ULAPACK_USE_STATIC_ALLOC
        ulapack_matrix_t A;
        ulapack_matrix_t B;
        ulapack_matrix_t result;

        ulapack_initialize_matrix(&A, 3u, 3u);
        ulapack_initialize_matrix(&B, 3u, 3u);
        ulapack_initialize_matrix(&result, 3u, 3u);

        ulapack_eye(&A);
        ulapack_scale(&A, 3, &A);

        ulapack_set(&B, 3.14);
        ulapack_add(&A, &B, &result);

        print_matrix(&result);

    #else
        ulapack_matrix_t *A = NULL;
        ulapack_matrix_t *B = NULL;
        ulapack_matrix_t *result = NULL;

        ulapack_initialize_matrix(&A, 3u, 3u);
        ulapack_initialize_matrix(&B, 3u, 3u);
        ulapack_initialize_matrix(&result, 3u, 3u);

        ulapack_eye(A);
        ulapack_scale(A, 3, A);

        ulapack_set(B, 3.14);
        ulapack_add(A, B, result);

        print_matrix(result);

    #endif

    return ulapack_success;
}

int main(void) {
    if (ut_iserror(test_initialization())) {
        return ulapack_error;
    }

    if (ut_iserror(test_addition())) {
        return ulapack_error;
    }

    return 0;
}


