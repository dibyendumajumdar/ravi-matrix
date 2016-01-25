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

#include <ravi_matrixlib.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__APPLE__) && !defined(USE_OPENBLAS)
extern void ddisna_(char *job, int *m, int *n, double *d, double *sep,
                    int *info);
extern void dgesv_(int *N, int *nrhs, double *A, int *lda, int *ipiv, double *b,
                   int *ldb, int *info);
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc);
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
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
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work,
                    int *lwork, int *info);
#elif defined(USE_OPENBLAS)
#include <openblas_config.h>
#include <f77blas.h>
extern void dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
                    int *ldb, double *s, double *rcond, int *rank, double *work,
                    int *lwork, int *iwork, int *info);
extern void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
                    int *ldb, int *jpvt, double *rcond, int *rank, double *work,
                    int *lwork, int *info);
extern void dgecon_(const char *norm, const int *n, double *a, const int *lda,
                    const double *anorm, double *rcond, double *work,
                    int *iwork, int *info);
extern double dlange_(const char *norm, const int *m, const int *n,
                      const double *a, const int *lda, double *work);
extern void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s,
                    double *u, int *ldu, double *vt, int *ldvt, double *work,
                    int *lwork, int *iwork, int *info);
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work,
                    int *lwork, int *info);
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif /* __APPLE__ */

#ifdef __cplusplus
}
#endif

// Internal workhorse for matrix multiplication
// workhorse for matrix multiplication
// C=alpha*A*B + beta*C
// dgemm wrapper
static void matrix_multiply(int32_t rows_A, int32_t cols_A, double *a, int32_t rows_B, int32_t cols_B, double *b,
  int32_t rows_C, int32_t cols_C, double *c, bool transposeA,
  bool transposeB, double alpha, double beta) {
  int m = transposeA ? cols_A : rows_A; /* If transposing then m = columns(A) else rows(A) */
  int n = transposeB ? rows_B : cols_B; /* If transposing then n = rows(B) else columns(B) */ 
  int k = transposeA ? rows_A : cols_A; /* If transposing A then k = rows(A) else columns(A) */
  assert(rows_C == m);
  assert(cols_C == n);
  int lda = rows_A;
  int ldb = rows_B;
  int ldc = rows_C;
  /*
  *  Definition:
  *  ===========
  *
  *       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  *
  *       .. Scalar Arguments ..
  *       DOUBLE PRECISION ALPHA,BETA
  *       INTEGER K,LDA,LDB,LDC,M,N
  *       CHARACTER TRANSA,TRANSB
  *       ..
  *       .. Array Arguments ..
  *       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
  *       ..
  *
  *
  * Purpose:
  * =============
  *
  * DGEMM  performs one of the matrix-matrix operations
  *
  *    C := alpha*op( A )*op( B ) + beta*C,
  *
  * where  op( X ) is one of
  *
  *    op( X ) = X   or   op( X ) = X**T,
  *
  * alpha and beta are scalars, and A, B and C are matrices, with op( A )
  * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  *
  * Arguments:
  * ==========
  *
  * [in] TRANSA
  *          TRANSA is CHARACTER*1
  *           On entry, TRANSA specifies the form of op( A ) to be used in
  *           the matrix multiplication as follows:
  *
  *              TRANSA = 'N' or 'n',  op( A ) = A.
  *
  *              TRANSA = 'T' or 't',  op( A ) = A**T.
  *
  *              TRANSA = 'C' or 'c',  op( A ) = A**T.
  *
  *  [in] TRANSB
  *          TRANSB is CHARACTER*1
  *           On entry, TRANSB specifies the form of op( B ) to be used in
  *           the matrix multiplication as follows:
  *
  *              TRANSB = 'N' or 'n',  op( B ) = B.
  *
  *              TRANSB = 'T' or 't',  op( B ) = B**T.
  *
  *              TRANSB = 'C' or 'c',  op( B ) = B**T.
  *
  *  [in] M
  *          M is INTEGER
  *           On entry,  M  specifies  the number  of rows  of the  matrix
  *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  *
  *  [in] N
  *          N is INTEGER
  *           On entry,  N  specifies the number  of columns of the matrix
  *           op( B ) and the number of columns of the matrix C. N must be
  *           at least zero.
  *
  *  [in] K
  *          K is INTEGER
  *           On entry,  K  specifies  the number of columns of the matrix
  *           op( A ) and the number of rows of the matrix op( B ). K must
  *           be at least  zero.
  *
  *  [in] ALPHA
  *          ALPHA is DOUBLE PRECISION.
  *           On entry, ALPHA specifies the scalar alpha.
  *
  *  [in] A
  *          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
  *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  *           part of the array  A  must contain the matrix  A,  otherwise
  *           the leading  k by m  part of the array  A  must contain  the
  *           matrix A.
  *
  *  [in] LDA
  *          LDA is INTEGER
  *           On entry, LDA specifies the first dimension of A as declared
  *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  *           least  max( 1, k ).
  *
  *  [in] B
  *          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
  *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  *           part of the array  B  must contain the matrix  B,  otherwise
  *           the leading  n by k  part of the array  B  must contain  the
  *           matrix B.
  *
  *  [in] LDB
  *          LDB is INTEGER
  *           On entry, LDB specifies the first dimension of B as declared
  *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  *           least  max( 1, n ).
  *
  *  [in] BETA
  *          BETA is DOUBLE PRECISION.
  *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  *           supplied as zero then C need not be set on input.
  *
  *  [in,out] C
  *          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
  *           Before entry, the leading  m by n  part of the array  C must
  *           contain the matrix  C,  except when  beta  is zero, in which
  *           case C need not be set on entry.
  *           On exit, the array  C  is overwritten by the  m by n  matrix
  *           ( alpha*op( A )*op( B ) + beta*C ).
  *
  *  [in] LDC
  *          LDC is INTEGER
  *           On entry, LDC specifies the first dimension of C as declared
  *           in  the  calling  (sub)  program.   LDC  must  be  at  least
  *           max( 1, m ).
  */
#if __APPLE__
  cblas_dgemm(CblasColMajor, transposeA ? CblasTrans : CblasNoTrans,
              transposeB ? CblasTrans : CblasNoTrans, m, n, k, alpha, a, lda, b,
              ldb, beta, c, ldc);
#else
  char transa = transposeA ? 'T' : 'N';
  char transb = transposeB ? 'T' : 'N';
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}

static double matrix_norm1(ravi_matrix_t *A) {
  int m = A->m;
  int n = A->n;
  char norm[1] = {'1'};
  double *a = &A->data[0];
  int lda = m;
  double *work = (double *)alloca(sizeof(double) * m);
  return dlange_(&norm[0], &m, &n, a, &lda, work);
}

static int matrix_lufactor(ravi_matrix_t *A) {
  int m = A->m;
  int n = A->n;

  int info = 0;
  int size = m < n ? m : n;
  int *ipiv = (int *)alloca(sizeof(int) * size);
  for (int i = 0; i < size; i++)
    ipiv[i] = 0;

  int lda = (1 < m? m: 1);
  dgetrf_(&m, &n, &A->data[0], &lda, ipiv, &info);
  if (info != 0)
    fprintf(stderr, "failed LU factorization of input matrix: %d\n", info);
  return info;
}

static void matrix_negate(ravi_matrix_t *A) {
  const double *e = &A->data[A->m * A->n];
  for (double *p = &A->data[0]; p != e; p++)
    *p = -*p;
}

static bool matrix_add(ravi_matrix_t *A, const ravi_matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr,
            "Matrix addition failed as matrices are not the same size\n");
    assert(false);
    return false;
  }
  const double *e = &A->data[A->m * A->n];
  double *p = &A->data[0];
  const double *q = &B->data[0];
  for (; p != e; p++, q++)
    *p += *q;
  return true;
}

static bool matrix_sub(ravi_matrix_t *A, const ravi_matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr,
            "Matrix subtraction failed as matrices are not the same size\n");
    assert(false);
    return false;
  }
  const double *e = &A->data[A->m * A->n];
  double *p = &A->data[0];
  const double *q = &B->data[0];
  for (; p != e; p++, q++)
    *p -= *q;
  return true;
}

// Copy B into A
static void matrix_copy(ravi_matrix_t *A, const ravi_matrix_t *B) {
  if (A->m != B->m || A->n != B->n) {
    fprintf(stderr, "Matrix copy failed as matrices are not the same size\n");
    assert(false);
  }
  const double *e = &A->data[A->m * A->n];
  double *p = &A->data[0];
  const double *q = &B->data[0];
  for (; p != e; p++, q++)
    *p = *q;
}

static ravi_matrix_t *make_copy(const ravi_matrix_t *src) {
  size_t msize = sizeof(int32_t) * 2 + sizeof(double) * src->m * src->n;
  ravi_matrix_t *A = calloc(1, msize);
  memcpy(A, src, msize);
  return A;
}

// SVD
// S must be matrix min(m,n) x 1 (vector)
// U must be matrix m x m
// V must be matrix n x n
static bool matrix_svd(const ravi_matrix_t *input, ravi_matrix_t *S, ravi_matrix_t *U,
                       ravi_matrix_t *V) {
  double workSize, *work = NULL, *a = NULL;
  int ldu, ldvt, lwork, info = -1, *iwork = NULL;
  char job = 'A';
  ravi_matrix_t *A = NULL;
  int m = input->m, n = input->n, lda = m;
  int ss = m < n ? m : n;
  int s_size = 1 < ss ? ss : 1;
  if (S->m != s_size || S->n != 1) {
    fprintf(stderr, "The vector S must be a column vector of size %d\n",
            s_size);
    goto done;
  }
  if (U->m != m || U->n != m) {
    fprintf(stderr, "The matrix U must be of size %dx%d\n", m, m);
    goto done;
  }
  ldu = U->m;
  if (V->m != n || V->n != n) {
    fprintf(stderr, "The matrix V must be of size %dx%d\n", n, n);
    goto done;
  }
  A = make_copy(input);
  if (!A) {
    fprintf(stderr, "Failed to allocate memory\n");
    goto done;
  }
  a = &A->data[0];
  ldvt = n;
  workSize = 0;
  work = &workSize;
  lwork = -1; // we want to estimate workSize first
  info = 0;
  int min_sz = 8 * (m < n ? m : n); // min(m,n)*8
  int works = 1 < min_sz ? min_sz : 1; // max
  iwork = (int *)calloc(works, sizeof(int));
  if (iwork == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    goto done;
  }
  dgesdd_(&job, &m, &n, a, &lda, &S->data[0], &U->data[0], &ldu, &V->data[0],
          &ldvt, work, &lwork, iwork, &info);
  if (info == 0) {
    lwork = (int)workSize;
    work = (double *)calloc(lwork, sizeof(double));
    if (work) {
      dgesdd_(&job, &m, &n, a, &lda, &S->data[0], &U->data[0], &ldu,
              &V->data[0], &ldvt, work, &lwork, iwork, &info);
    } else {
      fprintf(stderr, "Failed to allocate memory\n");
      info = -1;
    }
  } else {
    fprintf(stderr, "Failed to estimate work size for SVD: info=%d\n", info);
  }
  if (info != 0) {
    fprintf(stderr, "Failed to compute SVD: info=%d\n", info);
  }
done:
  if (A)
    free(A);
  if (iwork)
    free(iwork);
  if (work && work != &workSize)
    free(work);
  return info == 0;
}

static bool matrix_estimate_rcond(const ravi_matrix_t *A, double *rcond) {
  bool ok = false;
  ravi_matrix_t *copy_of_A = make_copy(A);
  double anorm = matrix_norm1(copy_of_A);
  int info = matrix_lufactor(copy_of_A);
  if (info != 0) {
    fprintf(stderr, "failed to estimate rcond (LU factor failed)\n");
    goto done;
  }
  char norm[] = {'1'};
  int n = A->n;
  double *a = &copy_of_A->data[0];
  int lda = A->m;
  if (lda < n) {
    fprintf(stderr, "failed to estimate rcond (LDA < n)\n");
    goto done;
  }
  double *work = (double *)alloca(sizeof(double) * n * 4);
  int *iwork = (int *)alloca(sizeof(int) * n);
  dgecon_(&norm[0], &n, a, &lda, &anorm, rcond, work, iwork, &info);
  if (info != 0) {
    fprintf(stderr, "failed to estimate rcond (DGECON failed)\n");
    goto done;
  }
  ok = true;
done:
  if (copy_of_A)
    free(copy_of_A);
  return ok;
}

static bool matrix_solve(ravi_matrix_t *m, ravi_matrix_t *v) {
  if (m->m == 0 || m->n == 0 || v->m == 0 || v->n != 1) {
    fprintf(stderr, "The matrix A and vector y must have rows > 0\n");
    assert(false);
    return false;
  }
  if (m->m != m->n || m->m != v->m) {
    fprintf(stderr, "The default solver only accepts n x n matrix\n");
    assert(false);
    return false;
  }

  int N = m->n;   // order of the matrix, number of rows of vector v
  int nrhs = 1;   // number of columns of vector v
  int lda = m->m; // leading dimension of the array for A (rows)
  int *ipiv =
      (int *)alloca(N * sizeof(int)); // return value: containing pivot indices
  memset(ipiv, 0, N * sizeof(int));

  int ldb = v->m; // leading dimension of vector v
  int info = 0;   // return value: if 0 successful
  double *A = &m->data[0];
  double *b = &v->data[0];

  dgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);

  return info == 0;
}

static bool matrix_lsq_solve(ravi_matrix_t *A, ravi_matrix_t *y, double rcond, bool svd) {
  if (A->m == 0 || A->n == 0 || y->m == 0 || y->n != 1) {
    fprintf(stderr, "The matrix A and vector y must have rows > 0\n");
    assert(false);
    return false;
  }
  if (A->m < A->n) {
    fprintf(stderr, "The matrix A must have rows >= cols\n");
    assert(false);
    return false;
  }
  if (y->m != A->m) {
    fprintf(stderr, "The vector y must have rows = A.rows\n");
    assert(false);
    return false;
  }
  if (rcond <= 0) {
    if (!matrix_estimate_rcond(A, &rcond))
      return false;
  }
  int m = A->m;
  int n = A->n;
  int nrhs = 1; // number of columns in b
  double *a = &A->data[0];
  int lda = m;
  double *b = &y->data[0];
  int ldb = m;
  int *jpvt = (int *)alloca(sizeof(int) * n);
  memset(jpvt, 0, sizeof(int) * n);
  double *s = (double *)alloca(sizeof(double) * (m > n ? m: n));
  int lwork = -1;
  int liwork = -1;
  int rank = 0;
  double temp = 0;
  int info = 0;
  // First estimate the work space required
  if (svd) {
    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, &temp, &lwork,
            &liwork, &info);
  } else {
    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, &temp, &lwork,
            &info);
  }
  if (info != 0) {
    fprintf(stderr, "failed to estimate work space requirement\n");
    return false;
  }
  // allocate work space
  lwork = (int)temp;
  double *ptr = (double *)calloc(lwork, sizeof(double));
  assert(ptr);
  // solve
  if (svd) {
    double minmn = m < n ? m : n; // min
    minmn = 1 < minmn ? minmn : 1; // max
    double smlsiz_plus_1 = 26.;
    double nlvl = (int)(log(minmn / smlsiz_plus_1) / log(.2)) + 1;
    nlvl = 0.0 < nlvl ? nlvl : 0.0; // max
    liwork = (int)(3.0 * minmn * nlvl + 11.0 * minmn);
    int *iwork_ptr = (int *)calloc(liwork, sizeof(int));
    assert(iwork_ptr);
    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, ptr, &lwork,
            iwork_ptr, &info);
    free(iwork_ptr);
  } else {
    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, ptr, &lwork,
            &info);
  }
  free(ptr);
  if (info != 0) {
    fprintf(stderr, "failed to solve linear least squares problem\n");
    return false;
  }
  return true;
}

static bool matrix_inverse(ravi_matrix_t *a) {
  int m = a->m;
  int n = a->n;
  if (m != n) {
    fprintf(stderr, "matrix is not square");
    return false;
  }
  int info = 0;
  int size = m < n ? m : n; //min
  int *ipiv = (int *)alloca(sizeof(int) * size);
  memset(ipiv, 0, sizeof(int) * size);
  int lda = (1 < m ? m : 1);
  dgetrf_(&m, &n, &a->data[0], &lda, ipiv, &info);
  if (info != 0) {
    fprintf(stderr, "failed LU factorization of input matrix");
    return false;
  }
  n = a->n;
  lda = n;
  int lwork = n * n;
  double *work = (double *)alloca(sizeof(double) * lwork);
  info = 0;
  dgetri_(&n, &a->data[0], &lda, ipiv, work, &lwork, &info);
  if (info != 0) {
    fprintf(stderr, "failed to compute inverse of matrix");
    return false;
  }
  return true;
}

static void matrix_transpose(ravi_matrix_t *result, const ravi_matrix_t *m) {
  assert(result->n == m->m && result->m == m->n);
  if (m->m <= 0 || m->n <= 0)
    return;
#if 0 // USE_OPENBLAS (for some reason link fails)
  char ordering = 'C';
  char transpose = 'T';
  double scale = 1.0;
  int rows = m->m;
  int cols = m->n;
  int lda = max(1, rows);
  int ldb = max(1, cols);
  domatcopy_(&ordering, &transpose, &rows, &cols, &scale, &m->data[0], &lda, &result->data[0], &ldb);
#else
  /*
    result is column order so we can use
    more optimimised iterator
  */
  double *rp = &result->data[0];
  for (int32_t i = 0; i < m->m; i++) {
    for (int32_t j = 0; j < m->n; j++) {
      *rp++ = m->data[j * m->m + i];
    }
  }
#endif
}

/////////////////////////////////////////////////////

static ravi_matrix_ops_t ops = {matrix_multiply, matrix_norm1, matrix_lufactor,
                           matrix_estimate_rcond, matrix_svd, matrix_negate,
                           matrix_add, matrix_sub, matrix_copy, matrix_solve,
                           matrix_lsq_solve, matrix_inverse, matrix_transpose};

const ravi_matrix_ops_t *ravi_matrix_get_implementation() { return &ops; }