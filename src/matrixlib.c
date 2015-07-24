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

#include <ravimatrix/matrixlib.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER

#include <malloc.h>
#define alloca _alloca

#ifndef __cplusplus
#define inline __inline
#endif

#else

#define DLLEXPORT
#define DLLIMPORT

#include <alloca.h>

#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __APPLE__
extern void dgesv_(int *N, int *nrhs, double *A, int *lda, int *ipiv, double *b,
                   int *ldb, int *info);
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc);
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work,
                    int *lwork, int *info);
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
                   int *lda, double *x, int *incx, double *beta, double *y,
                   int *incy);
extern void dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
                    int *ldb, double *s, double *rcond, int *rank, double *work,
                    int *lwork, int *iwork, int *info);
extern void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
                    int *ldb, int *jpvt, double *rcond, int *rank, double *work,
                    int *lwork, int *info);
extern void dgesvx_(char *fact, char *trans, int *n, int *nrhs, double *a,
                    int *lda, double *af, int *ldaf, int *ipiv, char *equed,
                    double *r__, double *c__, double *b, int *ldb, double *x,
                    int *ldx, double *rcond, double *ferr, double *berr,
                    double *work, int *iwork, int *info);
extern void dgecon_(const char *norm, const int *n, double *a, const int *lda,
                    const double *anorm, double *rcond, double *work,
                    int *iwork, int *info);
extern double dlange_(const char *norm, const int *m, const int *n,
                      const double *a, const int *lda, double *work);
extern void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s,
                    double *u, int *ldu, double *vt, int *ldvt, double *work,
                    int *lwork, int *iwork, int *info);
extern void ddisna_(char *job, int *m, int *n, double *d, double *sep,
                    int *info);
#else /* __APPLE__ */
#include <Accelerate/Accelerate.h>
#endif /* __APPLE__ */

#ifdef __cplusplus
}
#endif

// Internal workhorse for matrix multiplication
bool matrix_multiply(matrix_t *C, matrix_t *A, matrix_t *B, bool transposeA,
                     bool transposeB) {
  int m = A->m;
  int n = B->n;
  int k = A->n;
  if (C->m != m || C->n != n || B->m != k) {
    fprintf(stderr, "Dimensions are unexpected: A.m=%d, A.n=%d, B.m=%d, B.n=%d, C.m=%d, C.n=%d\n", A->m, A->n, B->m, B->n, C->m, C->n);
    assert(false);
    return false;
  }
  double alpha = 1.0;
  double beta = 0.0;

  double *c = &C->data[0];
  double *a = &A->data[0];
  double *b = &B->data[0];

#if __APPLE__
  cblas_dgemm(CblasColMajor, transposeA ? CblasTrans : CblasNoTrans,
              transposeB ? CblasTrans : CblasNoTrans, m, n, k, alpha, a, m, b,
              k, beta, c, m);
#else
  char transa = transposeA ? 'T' : 'N';
  char transb = transposeB ? 'T' : 'N';
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &m, b, &k, &beta, c, &m);
#endif
  return true;
}

double matrix_norm1(matrix_t *A) {
  int m = A->m;
  int n = A->n;
  char norm[1] = { '1' };
  double *a = &A->data[0];
  int lda = m;
  double *work = (double *)alloca(sizeof(double) * m);
  return dlange_(&norm[0], &m, &n, a, &lda, work);
}

int matrix_lufactor(matrix_t *A) {
  int m = A->m;
  int n = A->n;

  int info = 0;
  int size = min(m, n);
  int *ipiv = (int *)alloca(sizeof(int) * size);
  for (int i = 0; i < size; i++) ipiv[i] = 0;

  int lda = max(1, m);
  dgetrf_(&m, &n, &A->data[0], &lda, ipiv, &info);
  if (info != 0)
    fprintf(stderr, "failed LU factorization of input matrix: %d\n", info);
  return info;
}


void matrix_negate(matrix_t *A) {
  const double *e = &A->data[A->m*A->n];
  for (double *p = &A->data[0]; p != e; p++)
    *p = -*p;
}

bool matrix_add(matrix_t *A, matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr, "Matrix addition failed as matrices are not the same size\n");
    assert(false);
    return false;
  }
  const double *e = &A->data[A->m*A->n];
  for (double *p = &A->data[0], *q = &B->data[0]; p != e; p++, q++)
    *p += *q;
  return true;
}

bool matrix_sub(matrix_t *A, matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr, "Matrix subtraction failed as matrices are not the same size\n");
    assert(false);
    return false;
  }
  const double *e = &A->data[A->m*A->n];
  for (double *p = &A->data[0], *q = &B->data[0]; p != e; p++, q++)
    *p -= *q;
  return true;
}

// Copy B into A
void matrix_copy(matrix_t *A, matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr, "Matrix subtraction failed as matrices are not the same size\n");
    assert(false);
  }
  const double *e = &A->data[A->m*A->n];
  for (double *p = &A->data[0], *q = &B->data[0]; p != e; p++, q++)
    *p = *q;
}
