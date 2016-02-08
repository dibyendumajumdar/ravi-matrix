/******************************************************************************
* Copyright (C) 2015 Dibyendu Majumdar
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
******************************************************************************/

#ifndef RAVI_MATRIXLIB_H_
#define RAVI_MATRIXLIB_H_

#include <ravi_matrix_conf.h>

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 
* We will cast the Ravi numeric arrays into this format
* exploiting the fact that there is an extra element in Ravi
* arrays
* When not using Ravi this will be our vector or matrix
* user defined type
*/
typedef struct ravi_matrix_t ravi_matrix_t;
struct ravi_matrix_t {
  int32_t m; /* rows */
  int32_t n; /* columns */
  double data[1];
};

enum ravi_matrix_norm_type {
  RAVI_MATRIX_NORM_ONE,
  RAVI_MATRIX_NORM_INFINITY,
  RAVI_MATRIX_NORM_FROBENIUS
};
typedef enum ravi_matrix_norm_type ravi_matrix_norm_type;

typedef struct ravi_matrix_ops_t ravi_matrix_ops_t;
struct ravi_matrix_ops_t {
  // workhorse for matrix multiplication
  // C=alpha*A*B + beta*C
  // Simple wrapper around dgemm
  // rows_A = number of rows in A
  // cols_A = number of columns in A
  // rows_B = number of rows in B
  // cols_B = number of columns in B
  // rows_C = number of rows in C
  // cols_C = number of columns in C
  // To perform C=A*B set alpha to 1.0 and beta to 0.0
  void (*multiply)(int32_t rows_A, int32_t cols_A, double *A, int32_t rows_B, int32_t cols_B, double *B, 
                   int32_t rows_C, int32_t cols_C, double *C, bool transposeA,
                   bool transposeB, double alpha, double beta);

  double (*norm)(int32_t m, int32_t n, double *a, ravi_matrix_norm_type normType);

  // LU Factorisation
  // The matrix 'a' will be updated
  // ipsize must be at least min(m,n)
  // ipiv must be an array of size ipsize
  // Returns 0 on success
  int(*lufactor)(int32_t m, int32_t n, double *a, int32_t ipsize, int *ipiv);

  // A = A + alpha*B
  void(*add)(int32_t rows, int32_t cols, double *A, const double *B, double alpha);

  bool (*estimate_rcond)(const ravi_matrix_t *A, double *rcond);

  // SVD
  // S must be matrix min(m,n) x 1 (vector)
  // U must be matrix m x m
  // V must be matrix n x n
  bool (*svd)(const ravi_matrix_t *A, ravi_matrix_t *S, ravi_matrix_t *U, ravi_matrix_t *V);

  // A = -A
  void (*negate)(ravi_matrix_t *A);

  // A = copy(B)
  void (*copy)(ravi_matrix_t *A, const ravi_matrix_t *B);

  // M must be a square matrix
  // V a column vector with rows same as M
  bool (*solve)(ravi_matrix_t *M, ravi_matrix_t *V);

  // M must rows > columns
  // V a column vector with rows same as M
  bool (*lsq_solve)(ravi_matrix_t *M, ravi_matrix_t *V, double rcond, bool use_svd);

  // A must be square matrix
  bool (*inverse)(ravi_matrix_t *A);

  // transposed must be size nxm where original is sized mxn
  void(*transpose)(int32_t rows, int32_t cols, double *b, const double *a);
};

RAVIMATRIX_API const ravi_matrix_ops_t *ravi_matrix_get_implementation();

#ifdef __cplusplus
}
#endif

#endif
