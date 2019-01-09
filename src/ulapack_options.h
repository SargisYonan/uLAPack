/**
 * @file ulapack_options.h
 * @brief Settings profile for uLAPack operations and usage.
 *
 * @author Sargis Yonan
 * @date July 8, 2017
 *
 * @version 1.0.1
 **/

/*
 * File header guards.
 */
#ifndef __ULAPACK_OPTIONS_H__
#define __ULAPACK_OPTIONS_H__

/*
 * Define the entry data type.
 */
#define ULAPACK_MATRIX_ENTRY_TYPE double
typedef ULAPACK_MATRIX_ENTRY_TYPE ulapack_entry_t;

/*
 * Defined if memory should be initialized to zeros.
 */
#define ULAPACK_INITIALIZE_MEMORY

/*
 * Keep defined for static memory allocation declarations.
 */
// #define ULAPACK_USE_STATIC_ALLOC

#ifndef ULAPACK_USE_STATIC_ALLOC
    /*
     * Specify the memory allocator function to use which takes the form of stdlib's
     * malloc function : void *alloc(size_t bytes).
     */
    #include <stdlib.h>
    #define ULAPACK_ALLOC malloc

    /*
     * Specify the memory freeing function to use which takes the form of stdlib's
     * free function : void free(void* ptr).
     */
    #include <stdlib.h>
    #define ULAPACK_FREE free
#endif

#ifdef ULAPACK_USE_STATIC_ALLOC
    /*
     * Set the largest row and column sizes for uLAPack matrix types.
     *
     * @note only required to be specified when static allocation is used.
     */
    #define ULAPACK_MAX_MATRIX_N_ROWS (10u)
    #define ULAPACK_MAX_MATRIX_N_COLS (10u)

    /*
     * The maximum number of rows and columns must be specified
     * for static matrix allocation.
     */
    #ifndef ULAPACK_MAX_MATRIX_N_ROWS
        #error "ULAPACK_MAX_MATRIX_N_ROWS not defined."
    #endif
    #ifndef ULAPACK_MAX_MATRIX_N_COLS
        #error "ULAPACK_MAX_MATRIX_N_COLS not defined."
    #endif
#endif

#ifndef ULAPACK_USE_STATIC_ALLOC
    #ifndef ULAPACK_ALLOC
        #error "ULAPACK_ALLOC is not defined."
    #endif
#endif

/*
 * End header guard definition.
 */
#endif