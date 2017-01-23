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

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <ravi_matrix_conf.h>
#include <ravi_linalg.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ravi_vector_t ravi_vector_;
struct ravi_vector_t {
  int32_t n; /* columns */
  double* data;
};

typedef struct ravi_matrix_t ravi_matrix_t;
struct ravi_matrix_t {
  int32_t m; /* rows */
  int32_t n; /* columns */
  double* data;
};

enum ravi_matrix_norm_type {
  RAVI_MATRIX_NORM_ONE,
  RAVI_MATRIX_NORM_INFINITY,
  RAVI_MATRIX_NORM_FROBENIUS
};
typedef enum ravi_matrix_norm_type ravi_matrix_norm_type;

// workhorse for matrix multiplication
// C=alpha*A*B + beta*C
// Simple wrapper around dgemm
// To perform C=A*B set alpha to 1.0 and beta to 0.0
static inline void ravi_matrix_multiply(ravi_matrix_t* A,
                                        ravi_matrix_t* B,
                                        ravi_matrix_t* C,
                                        bool transposeA,
                                        bool transposeB,
                                        double alpha,
                                        double beta) {
  int m = transposeA
              ? A->n
              : A->m; /* If transposing then m = columns(A) else rows(A) */
  int n = transposeB
              ? B->m
              : B->n; /* If transposing then n = rows(B) else columns(B) */
  int k = transposeA
              ? A->m
              : A->n; /* If transposing A then k = rows(A) else columns(A) */
  assert(C->m == m);
  assert(C->n == n);
  int lda = A->m;
  int ldb = B->m;
  int ldc = C->m;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
  cblas_dgemm(CblasColMajor, transposeA ? CblasTrans : CblasNoTrans,
              transposeB ? CblasTrans : CblasNoTrans, m, n, k, alpha, A->data,
              lda, B->data, ldb, beta, C->data, ldc);
#else
  char transa = transposeA ? 'T' : 'N';
  char transb = transposeB ? 'T' : 'N';
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A->data, &lda, B->data, &ldb,
         &beta, C->data, &ldc);
#endif
}

static double ravi_matrix_norm(ravi_matrix_t* A,
                               ravi_matrix_norm_type normType) {
  char normt = normType == RAVI_MATRIX_NORM_ONE
                   ? '1'
                   : (normType == RAVI_MATRIX_NORM_INFINITY ? 'i' : 'f');
  int lda = A->m;
  double* work = (double*)alloca(sizeof(double) * A->m);
  return dlange_(&normt, &A->m, &A->n, A->data, &lda, work);
}

// LU Factorisation
// The matrix 'A' will be updated
// ipsize must be at least min(m,n)
// ipiv must be an array of size ipsize
// Returns 0 on success
static inline int ravi_matrix_lufactor(ravi_matrix_t* A,
                                       int32_t ipsize,
                                       int* ipiv) {
  int info = 0;
  int size = A->m < A->n ? A->m : A->n;
  if (ipsize < size) {
    info = -5;
    goto Lerror;
  }
  for (int i = 0; i < ipsize; i++)
    ipiv[i] = 0;
  int lda = (1 < A->m ? A->m : 1);
  dgetrf_(&A->m, &A->n, A->data, &lda, ipiv, &info);
Lerror:
  return info;
}

// A = A*scalar
static inline void ravi_matrix_scalar_multiply(ravi_matrix_t* A,
                                               double scalar) {
  int n = A->m * A->n;
  int inca = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
  cblas_dscal(n, scalar, A->data, inca);
#else
  dscal_(&n, &scalar, A->data, &inca);
#endif
}

// A = A + alpha*B
static inline void ravi_matrix_add(ravi_matrix_t* A,
                                   ravi_matrix_t* B,
                                   double alpha) {
  int n = A->m * A->n;
  int incx = 1;
  int incy = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
  cblas_daxpy(n, alpha, B->data, incx, A->data, incy);
#elif 0
  const double* e = &A->data[A->m * A->n];
  double* p = A->data;
  const double* q = B->data;
  for (; p != e; p++, q++)
    *p += alpha * (*q);
#else
  daxpy_(&n, &alpha, B->data, &incx, A->data, &incy);
#endif
}

// Copy B into A
// A = copy(B)
static inline void ravi_matrix_copy(ravi_matrix_t* A, const ravi_matrix_t* B) {
  assert(A->m == B->m && A->n == B->n);
  const double* e = &A->data[A->m * A->n];
  double* p = &A->data[0];
  const double* q = &B->data[0];
  for (; p != e; p++, q++)
    *p = *q;
}

extern bool ravi_matrix_estimate_rcond(const ravi_matrix_t* A, double* rcond);

// SVD
// S must be matrix min(m,n) x 1 (vector)
// U must be matrix m x m
// V must be matrix n x n
extern bool ravi_matrix_svd(const ravi_matrix_t* A,
                            ravi_matrix_t* S,
                            ravi_matrix_t* U,
                            ravi_matrix_t* V);

// M must be a square matrix
// V a column vector with rows same as M
extern bool ravi_matrix_solve(ravi_matrix_t* M, ravi_matrix_t* V);

// M must rows > columns
// V a column vector with rows same as M
bool ravi_matrix_lsq_solve(ravi_matrix_t* M,
                           ravi_matrix_t* V,
                           double rcond,
                           bool use_svd);

// A must be square matrix
bool ravi_matrix_inverse(ravi_matrix_t* A);

// transposed must be size nxm where original is sized mxn
// B = transpose(A)
static void ravi_matrix_transpose(ravi_matrix_t* B, const ravi_matrix_t* A) {
#if defined(USE_OPENBLAS)
  double scale = 1.0;
  int lda = A->m > 1 ? A->m : 1;  // max(1, rows);
  int ldb = A->n > 1 ? A->n : 1;  // max(1, cols);
  cblas_domatcopy(CblasColMajor, CblasTrans, A->m, A->n, scale, A->data, lda,
                  B->data, ldb);
#elif defined(USE_MKL)
  double scale = 1.0;
  int lda = A->m > 1 ? A->m : 1;  // max(1, rows);
  int ldb = A->n > 1 ? A->n : 1;  // max(1, cols);
  mkl_domatcopy(CblasColMajor, CblasTrans, A->m, A->n, scale, A->data, lda,
                B->data, ldb);
#else
  /*
  result is column order so we can use
  more optimised iterator
  */
  double* rp = B->data;
  for (int32_t i = 0; i < A->m; i++) {
    for (int32_t j = 0; j < A->n; j++) {
      *rp++ = A->data[j * A->m + i];
    }
  }
#endif
  B->m = A->n;
  B->n = A->m;
}

static inline void ravi_vector_outer_product(int32_t m,
                                             const double* x,
                                             int32_t n,
                                             const double* y,
                                             double* a,
                                             double alpha) {
  int incx = 1, incy = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
  cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, a, m);
#else
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &m);
#endif
}

#ifdef __cplusplus
}
#endif

#endif
