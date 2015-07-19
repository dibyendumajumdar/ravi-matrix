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

#include <ravimatrix/matrix.h>

#if __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#ifndef __APPLE__
#ifdef __cplusplus
extern "C" {
#endif
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

#ifdef __cplusplus
}
#endif

#endif

#if LUA_VERSION_NUM < 502

// Forward compatibility with Lua 5.2
// Following is not defined in 5.1

#define lua_absindex(L, i)                                                     \
  ((i) > 0 || (i) <= LUA_REGISTRYINDEX ? (i) : lua_gettop(L) + (i) + 1)

#endif

#if (LUA_VERSION_NUM >= 502)

// backward compatibility with Lua 5.1
#undef lua_equal
#define lua_equal(L, idx1, idx2) lua_compare(L, (idx1), (idx2), LUA_OPEQ)

#undef lua_objlen
#define lua_objlen lua_rawlen

#endif

/* We will cast the Ravi numeric arrays into this format
 * exploiting the fact that there is an extra element in Ravi
 * arrays
 * When not using Ravi this will be our vector or matrix
 * user defined type
 */
typedef struct matrix_t matrix_t;
struct matrix_t {
  int m;
  int n;
  double data[1];
};

// The normal Lua metatable functions in C use string
// keys - these are expensive as the key needs to be
// converted to Lua string, hash code computed etc.
// Following implementations are taken from a post in
// Lua mailing list (http://lua-users.org/lists/lua-l/2010-11/msg00151.html)
// They use lightuserdata instead of strings to speed
// things up
// meta_key is the key assigned to the meta table of the userdata
int l_newmetatable(lua_State *L, const char *meta_key) {
  lua_pushlightuserdata(L, (void *)meta_key);
  lua_rawget(L, LUA_REGISTRYINDEX);
  if (!lua_isnil(L, -1)) // name already in use?
    return 0;            // leave previous value on top, but return 0
  lua_pop(L, 1);         // pop the nil value
  lua_newtable(L);       // create metatable
  lua_pushlightuserdata(L, (void *)meta_key); // meta_key
  lua_pushvalue(L, -2);                       // table
  lua_rawset(L, LUA_REGISTRYINDEX); // assign table to meta_key in the registry
  return 1;
}

// meta_key is the key assigned to the meta table of the userdata
void l_getmetatable(lua_State *L, const char *meta_key) {
  lua_pushlightuserdata(L, (void *)meta_key); // meta_key
  lua_rawget(L, LUA_REGISTRYINDEX); // obtain the value associated with
                                    // meta_key from registry
}

// arg_index is the position of userdata argument on the stack
// meta_key is the key assigned to the meta table of the userdata
void *l_testudata(lua_State *L, int arg_index, const char *meta_key) {
  void *p = lua_touserdata(L, arg_index);
  if (p != NULL) {                                // value is a userdata?
    if (lua_getmetatable(L, arg_index)) {         // does it have a metatable?
      lua_pushlightuserdata(L, (void *)meta_key); // meta_key
      lua_rawget(
          L,
          LUA_REGISTRYINDEX); // get correct metatable associated with meta_key
      if (lua_rawequal(L, -1, -2)) { // compare: does it have the correct mt?
        lua_pop(L, 2);               // remove both metatables
        return p;                    // test ok
      }
    }
  }
  return NULL; /* to avoid warnings */
}

// arg_index is the position of userdata argument on the stack
// meta_key is the key assigned to the meta table of the userdata
void *l_checkudata(lua_State *L, int arg_index, const char *meta_key) {
  void *p = l_testudata(L, arg_index, meta_key);
  if (p == NULL)
    luaL_argerror(L, arg_index, meta_key);
  return p;
}

// We need fixed pointer values for metatable keys
static const char *Ravi_matrix = "Ravi.matrix";
static const char *Ravi_matrixx = "Ravi.matrixx";

static inline matrix_t *alloc_matrix(lua_State *L, int m, int n) {
  matrix_t *vector =
      (matrix_t *)lua_newuserdata(L, sizeof(double) * (m * n + 1));
  luaL_getmetatable(L, Ravi_matrix);
  lua_setmetatable(L, -2);
  return vector;
}

// arg_index is the position of userdata argument on the stack
// meta_key is the key assigned to the meta table of the userdata
matrix_t *test_matrixx(lua_State *L, int arg_index) {
  if (ravi_is_numberarray(L, arg_index)) {    // value is an array?
    if (lua_getmetatable(L, arg_index)) {     // does it have a metatable?
      lua_pushlightuserdata(L, (void *) Ravi_matrixx); // meta_key
      lua_rawget(
          L,
          LUA_REGISTRYINDEX); // get correct metatable associated with meta_key
      if (lua_rawequal(L, -1, -2)) { // compare: does it have the correct mt?
        lua_pop(L, 2);               // remove both metatables
        double *s, *e;
        ravi_get_numberarray_rawdata(L, arg_index, &s, &e);
        return (matrix_t *)s; // test ok
      }
    }
  }
  return NULL; /* to avoid warnings */
}

static int make_vector(lua_State *L) {
  if (lua_isnumber(L, 1)) {
    lua_Integer size = lua_tointeger(L, 1);
    lua_Number initv = 0.0;
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    matrix_t *vector = alloc_matrix(L, 1, (int)size);
    vector->m = 1;
    vector->n = (int)size;
    for (int i = 0; i < vector->n; i++)
      vector->data[i] = initv;
    return 1;
  }
  luaL_argerror(L, 1, "Cannot create to vector from supplied arguments");
  return 0;
}

/**
 * Constructs a matrix
 */
static int make_matrix(lua_State *L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");

    matrix_t *matrix = alloc_matrix(L, rows, cols);
    matrix->m = rows;
    matrix->n = cols;

    if (lua_istable(L, 3)) {
      int size = (int)lua_objlen(L, 3);
      luaL_argcheck(L, size == rows * cols, 3,
                    "data supplied is not of required length");
      for (int x = 0; x < size; x++) {
        lua_rawgeti(L, 3, x + 1);
        double v = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        matrix->data[x] = v;
      }
    } else {
      lua_Number initv = 0.0;
      if (lua_isnumber(L, 3))
        initv = lua_tonumber(L, 3);
      int size = rows * cols;
      for (int x = 0; x < size; x++) {
        matrix->data[x] = initv;
      }
    }
    return 1;
  } else if (lua_istable(L, 1)) {
    matrix_t *m = NULL;
    // we expect an array for each column
    int cols = (int)lua_objlen(L, 1);
    int rows = -1; // don't know how many rows yet
    if (cols == 0)
      rows = 0;
    else {
      // test first col
      lua_rawgeti(L, 1, 1);
      if (lua_istable(L, -1)) {
        // all rows are expected to be same length
        rows = (int)lua_objlen(L, -1);
      } else {
        luaL_argerror(L, 1, "expecting (table) array of arrays");
      }
      lua_pop(L, 1);
    }
    m = alloc_matrix(L, rows, cols);
    for (int j = 1; j <= cols; j++) {
      lua_rawgeti(L, 1, j);
      if (lua_istable(L, -1)) {
        luaL_argcheck(L, rows == lua_objlen(L, -1), 1,
                      "columns are not the same size");
        for (int k = 1; k <= rows; k++) {
          lua_rawgeti(L, -1, k);
          double v = luaL_checknumber(L, -1);
          lua_pop(L, 1);
          int pos = (j - 1) * rows + (k - 1);
          assert(pos >= 0 && pos < (rows * cols));
          m->data[pos] = v;
        }
      } else {
        luaL_argerror(L, 1, "expecting (table) array of arrays");
      }
      lua_pop(L, 1);
    }
    return 1;
  }
  luaL_argcheck(L, false, 1, "Unexpected arguments");
  return 0;
}

static inline matrix_t *set_matrixx_meta(lua_State *L, int rows, int cols) {
  l_getmetatable(L, Ravi_matrixx);
  lua_setmetatable(L, -2);
  double *start, *end;
  ravi_get_numberarray_rawdata(L, -1, &start, &end);
  matrix_t *v = (matrix_t *)start;
  assert(&v->data[0] == start + 1);
  // assert(&v->data[rows*cols] == end);
  v->m = rows;
  v->n = cols;
  return v;
}

static int make_matrixx(lua_State *L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");
    if (lua_gettop(L) == 3 && ravi_is_numberarray(L, 3)) {
      int size = (int)lua_objlen(L, 3);
      luaL_argcheck(L, size == rows * cols, 3,
                    "data supplied is not of required length");
    } else {
      lua_Number initv = 0.0;
      if (lua_isnumber(L, 3))
        initv = lua_tonumber(L, 3);
      int size = rows * cols;
      ravi_createnumberarray(L, size, initv);
    }
    set_matrixx_meta(L, rows, cols);
    lua_pushvalue(L, -1);
    return 1;
  }
  luaL_argcheck(L, false, 1, "Unexpected arguments");
  return 0;
}

/* make a num array vector (matrix with 1 row) */
static int make_vectorx(lua_State *L) {
  if (lua_isnumber(L, 1)) {
    int size = (int)lua_tointeger(L, 1);
    lua_Number initv = 0.0;
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    ravi_createnumberarray(L, size, initv);
    set_matrixx_meta(L, 1, size);
    return 1;
  } else if (lua_gettop(L) == 1 && ravi_is_numberarray(L, 1)) {
    int size = (int)lua_objlen(L, 1);
    set_matrixx_meta(L, 1, size);
    lua_pushvalue(L, -1);
    return 1;
  }
  luaL_argerror(L, 1, "Cannot create to vector from supplied arguments");
  return 0;
}

/* obtain a particular element of a userdata matrix */
static int vector_get(lua_State *L) {
  int nargs = lua_gettop(L);
  if (nargs == 0)
    return 0;
  matrix_t *vector = NULL;
  vector = (matrix_t *)luaL_checkudata(L, 1, Ravi_matrix);
  int pos = (int)luaL_checkinteger(L, 2);
  luaL_argcheck(L, pos >= 1 && pos <= (vector->m * vector->n), 2,
                "access out of bounds");
  lua_pushnumber(L, vector->data[pos - 1]);
  return 1;
}

/* compute array length of a userdata matrix */
static int vector_len(lua_State *L) {
  matrix_t *vector = (matrix_t *)l_testudata(L, 1, Ravi_matrix);
  if (vector) {
    lua_pushinteger(L, vector->m * vector->n);
    return 1;
  }
  luaL_argerror(L, 1, "A vector expected");
  return 0;
}

bool matrix_multiply(matrix_t *C, matrix_t *A, matrix_t *B, bool transposeA,
                     bool transposeB) {
  int m = A->m;
  int n = B->n;
  int k = A->n;
  if (C->m != m || C->n != n || B->m != k) {
    fprintf(stderr, "dimensions are unexpected");
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

/* multipy two userdata matrices */
static int matrix_mult(lua_State *L) {
  matrix_t *A = (matrix_t *)l_testudata(L, 1, Ravi_matrix);
  if (A) {
    matrix_t *B = (matrix_t *)l_testudata(L, 2, Ravi_matrix);
    if (B) {
      matrix_t *matrix = alloc_matrix(L, A->m, B->n);
      if (!matrix_multiply(matrix, A, B, false, false)) {
        luaL_error(L, "matrix multiplication failed");
        return 0;
      }
      return 1;
    }
  }
  luaL_argerror(L, 1, "Bad arguments");
  return 0;
}

/* multiply two number array matrices */
static int matrixx_mult(lua_State *L) {
  matrix_t *A = (matrix_t *)test_matrixx(L, 1);
  if (A) {
    matrix_t *B = (matrix_t *)test_matrixx(L, 2);
    if (B) {
      ravi_createnumberarray(L, A->m * B->n, 0.0);
      matrix_t *matrix = set_matrixx_meta(L, A->m, B->n);
      if (!matrix_multiply(matrix, A, B, false, false)) {
        luaL_error(L, "matrix multiplication failed");
        return 0;
      }
      return 1;
    }
  }
  luaL_argerror(L, 1, "Bad arguments");
  return 0;
}

// adds to an existing table
// table must be on stop of the stack
void add_to_library(lua_State *L, const struct luaL_Reg *regs) {
  for (int i = 0; regs[i].name != NULL; i++) {
    lua_pushcclosure(L, regs[i].func, 0);
    lua_setfield(L, -2, regs[i].name);
  }
}

// creates a table of functions which is a library in Lua
// terms. We use this as a portable way across 5.1 and 5.2 which
// follows the latest conventions - i.e. avoid polluting global
// namespace
void create_library(lua_State *L, const struct luaL_Reg *regs) {
  int count = 0;
  for (int i = 0; regs[i].name != NULL; i++) {
    count++;
  }
  lua_createtable(L, 0, count);
  add_to_library(L, regs);
  // leave table on the stack
}

static const struct luaL_Reg mylib[] = {{"vector", make_vector},
                                        {"vectorx", make_vectorx},
                                        {"matrix", make_matrix},
                                        {"matrixx", make_matrixx},
                                        {NULL, NULL}};

int luaopen_ravimatrix(lua_State *L) {
  fprintf(stderr, "Initializing RaviMatrix\n");
  l_newmetatable(L, Ravi_matrix);
  lua_pushstring(L, "Ravi.matrix");
  lua_setfield(L, -2, "__name");
  lua_pushcfunction(L, vector_get);
  lua_setfield(L, -2, "__index");
  lua_pushcfunction(L, vector_len);
  lua_setfield(L, -2, "__len");
  lua_pushcfunction(L, matrix_mult);
  lua_setfield(L, -2, "__mul");
  lua_pop(L, 1);

  l_newmetatable(L, Ravi_matrixx);
  lua_pushstring(L, "Ravi.matrixx");
  lua_setfield(L, -2, "__name");
  lua_pushcfunction(L, matrixx_mult);
  lua_setfield(L, -2, "__mul");
  lua_pop(L, 1);

  create_library(L, mylib);
  fprintf(stdout, "RaviMatrix initialized successfully\n");
  return 1;
}
