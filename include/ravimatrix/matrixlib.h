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

#ifndef MATRIXLIB_H_
#define MATRIXLIB_H_

#include <stdbool.h>
#include <stdint.h>

#include <ravimatrix/conf.h>

#ifdef __cplusplus
extern "C" {
#endif

/* We will cast the Ravi numeric arrays into this format
* exploiting the fact that there is an extra element in Ravi
* arrays
* When not using Ravi this will be our vector or matrix
* user defined type
*/
typedef struct matrix_t matrix_t;
struct matrix_t {
  int32_t m; /* rows */
  int32_t n; /* columns */
  double data[1];
};

typedef struct matrix_ops_t matrix_ops_t;
struct matrix_ops_t {
  // workhorse for matrix multiplication
  // C=A*B
  bool (*multiply)(matrix_t *C, matrix_t *A, matrix_t *B, bool transposeA,
                   bool transposeB);

  double (*norm1)(matrix_t *A);

  int (*lufactor)(matrix_t *A);

  bool (*estimate_rcond)(const matrix_t *A, double *rcond);

  // SVD
  // S must be matrix min(m,n) x 1 (vector)
  // U must be matrix m x m
  // V must be matrix n x n
  bool (*svd)(const matrix_t *A, matrix_t *S, matrix_t *U, matrix_t *V);

  // A = -A
  void (*negate)(matrix_t *A);

  // A += B
  bool (*add)(matrix_t *A, const matrix_t *B);

  // A -= B
  bool (*sub)(matrix_t *A, const matrix_t *B);

  // A = copy(B)
  void (*copy)(matrix_t *A, const matrix_t *B);

  // M must be a square matrix
  // V a column vector with rows same as M
  bool (*solve)(matrix_t *M, matrix_t *V);

  // M must rows > columns
  // V a column vector with rows same as M
  bool (*lsq_solve)(matrix_t *M, matrix_t *V, double rcond, bool use_svd);

  // A must be square matrix
  bool (*inverse)(matrix_t *A);

  // transposed must be size nxm where original is sized mxn
  void(*transpose)(matrix_t *transposed, matrix_t *original);
};

API const matrix_ops_t *ravi_matrix_get_implementation();

#ifdef __cplusplus
}
#endif

#endif
