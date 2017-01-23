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

#include <ravi_matrix.h>
#include <ravi_matrixlib.h>

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

// We need fixed pointer values for metatable keys
static const char* Lua_matrix = "LuaMatrix";
static const char* Ravi_matrix = "RaviMatrix";

static bool test_Lua_matrix(lua_State* L, int idx, ravi_matrix_t* m) {
  m = (ravi_matrix_t*)raviL_testudata(L, idx, Lua_matrix);
  return m != NULL;
}
ravi_matrix_t check_Lua_matrix(lua_State* L, int idx) {
  ravi_matrix_t* m = (ravi_matrix_t*)raviL_checkudata(L, idx, Lua_matrix);
  return *m;
}

// Allocate a Lua matrix as a userdata (compatible with Lua 5.3)
static ravi_matrix_t alloc_Lua_matrix(lua_State* L,
                                      int m,
                                      int n,
                                      double initv) {
  (void)initv;
  assert(m >= 0 && n >= 0);
  ravi_matrix_t* matrix =
      (ravi_matrix_t*)lua_newuserdata(L, sizeof(ravi_matrix_t) + m * n * sizeof(double));
  matrix->m = m;
  matrix->n = n;
  matrix->data = (double *)(matrix + 1);
  raviL_getmetatable(L, Lua_matrix);
  lua_setmetatable(L, -2);
  return *matrix;
}

#if RAVI_ENABLED
// Sets the metatable and other values for Ravi matrix
// Stack top must be a Ravi number array
static ravi_matrix_t set_Ravi_matrix_meta(lua_State* L, int m, int n) {
  assert(m >= 0 && n >= 0);
  raviL_getmetatable(L, Ravi_matrix);
  lua_setmetatable(L, -2);
  size_t len = 0;
  double* start = ravi_get_number_array_rawdata(L, -1, &len);
  ravi_matrix_t* dummy = (ravi_matrix_t*)start;
  dummy->m = m;
  dummy->n = n;
  ravi_matrix_t matrix = {.m = m, .n = n, .data = start + 1};
  return matrix;
}

static ravi_matrix_t alloc_Ravi_matrix(lua_State* L,
                                       int m,
                                       int n,
                                       double initv) {
  ravi_create_number_array(L, m * n, initv);
  return set_Ravi_matrix_meta(L, m, n);
}

// Test that the argument is a Ravi_matrix
// arg_index is the position of argument on the stack
static bool test_Ravi_matrix(lua_State* L, int arg_index, ravi_matrix_t* m) {
  if (ravi_is_number_array(L, arg_index)) {  // value is a Ravi array?
    if (lua_getmetatable(L, arg_index)) {    // does it have a metatable?
      raviL_getmetatable(L, Ravi_matrix);    // Get metatable for Ravi matrices
      if (lua_rawequal(L, -1, -2)) {  // compare: does it have the correct mt?
        lua_pop(L, 2);                // remove both metatables
        size_t len = 0;
        double* start = ravi_get_number_array_rawdata(L, arg_index, &len);
        ravi_matrix_t* dummy = (ravi_matrix_t*)start;
        m->m = dummy->m;
        m->n = dummy->n;
        m->data = start + 1;
        return true;
      }
    }
  }
  return false; /* to avoid warnings */
}

static ravi_matrix_t check_Ravi_matrix(lua_State* L, int arg_index) {
  ravi_matrix_t p;
  if (!test_Ravi_matrix(L, arg_index, &p))
    luaL_argerror(L, arg_index, Ravi_matrix);
  return p;
}
#endif

// Create a Lua vector (Lua matrix with one column)
// Interface 1
//   arg1 - size
//   arg2 - initial value (optional)
// Interface 2
//   arg1 - table
static int make_Lua_vector(lua_State* L) {
  if (lua_isnumber(L, 1)) {
    int size = (int)lua_tointeger(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    lua_Number initv = 0.0;  // initial value
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    // allocate matrix of 1 column
    ravi_matrix_t vector = alloc_Lua_matrix(L, (int)size, 1, initv);
    for (int i = 0; i < vector.m; i++)
      vector.data[i] = initv;
    return 1;
  } else if (lua_istable(L, 1)) {
    int size = (int)lua_rawlen(L, 1);
    // allocate matrix of 1 column
    ravi_matrix_t vector = alloc_Lua_matrix(L, (int)size, 1, 0.0);
    for (int x = 0; x < size; x++) {
      lua_rawgeti(L, 1, x + 1);
      double v = luaL_checknumber(L, -1);
      lua_pop(L, 1);
      vector.data[x] = v;
    }
    return 1;
  }
  luaL_argerror(L, 1, "Cannot create to vector from supplied arguments");
  return 0;
}
// Create a Lua matrix
// Interface 1
//   arg1 - rows
//   arg2 - columns
//   arg3 - initial value (optional)
// Interface 2
//   arg1 - rows
//   arg2 - columns
//   arg3 - table of values in column order (optional)
// Interface 3
//   arg1 - table of tables - each inner table represents a column
static int make_Lua_matrix(lua_State* L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");

    ravi_matrix_t matrix = alloc_Lua_matrix(L, rows, cols, 0.0);
    if (lua_istable(L, 3)) {
      int size = (int)lua_rawlen(L, 3);
      luaL_argcheck(L, size == rows * cols, 3,
                    "Table supplied is not of required length");
      for (int x = 0; x < size; x++) {
        lua_rawgeti(L, 3, x + 1);
        double v = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        matrix.data[x] = v;
      }
    } else {
      lua_Number initv = 0.0;
      if (lua_isnumber(L, 3))
        initv = lua_tonumber(L, 3);
      int size = rows * cols;
      for (int x = 0; x < size; x++) {
        matrix.data[x] = initv;
      }
    }
    return 1;
  } else if (lua_istable(L, 1)) {
    ravi_matrix_t m;
    // we expect an array for each column
    int cols = (int)lua_rawlen(L, 1);
    int rows = -1;  // don't know how many rows yet
    luaL_argcheck(L, cols > 0, 1, "Table must have at least one column table");
    // test first col
    lua_rawgeti(L, 1, 1);
    if (lua_istable(L, -1)) {
      // all rows are expected to be same length
      rows = (int)lua_rawlen(L, -1);
    } else {
      luaL_argerror(L, 1, "Expecting (table) array of arrays");
    }
    lua_pop(L, 1);  // pop the test value
    m = alloc_Lua_matrix(L, rows, cols, 0.0);
    for (int j = 1; j <= cols; j++) {
      lua_rawgeti(L, 1, j);
      if (lua_istable(L, -1)) {
        luaL_argcheck(L, rows == (int)lua_rawlen(L, -1), 1,
                      "columns are not the same size");
        for (int k = 1; k <= rows; k++) {
          lua_rawgeti(L, -1, k);
          double v = luaL_checknumber(L, -1);
          lua_pop(L, 1);
          int pos = (j - 1) * rows + (k - 1);
          assert(pos >= 0 && pos < (rows * cols));
          m.data[pos] = v;
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

#if RAVI_ENABLED

// Create a Ravi matrix
// Interface 1
//   arg1 - rows
//   arg2 - cols
//   arg3 - initial value
// Interface 2
//   arg1 - rows
//   arg2 - cols
//   arg3 - Ravi number array that will be converted
// Interface 3
//   arg1 - table of tables - each inner table represents a column
static int make_Ravi_matrix(lua_State* L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");
    if (lua_gettop(L) == 3 && ravi_is_number_array(L, 3)) {
      int size = (int)lua_rawlen(L, 3);
      luaL_argcheck(L, size == rows * cols, 3,
                    "data supplied is not of required length");
    } else {
      lua_Number initv = 0.0;
      if (lua_isnumber(L, 3))
        initv = lua_tonumber(L, 3);
      int size = rows * cols;
      ravi_create_number_array(L, size, initv);
    }
    set_Ravi_matrix_meta(L, rows, cols);
    lua_pushvalue(L, -1);
    return 1;
  } else if (lua_istable(L, 1)) {
    ravi_matrix_t m;
    // we expect an array for each column
    int cols = (int)lua_rawlen(L, 1);
    int rows = -1;  // don't know how many rows yet
    luaL_argcheck(L, cols > 0, 1, "Table must have at least one column table");
    // test first col
    lua_rawgeti(L, 1, 1);
    if (lua_istable(L, -1)) {
      // all rows are expected to be same length
      rows = (int)lua_rawlen(L, -1);
    } else {
      luaL_argerror(L, 1, "Expecting (table) array of arrays");
    }
    lua_pop(L, 1);  // pop the test value
    m = alloc_Ravi_matrix(L, rows, cols, 0.0);
    for (int j = 1; j <= cols; j++) {
      lua_rawgeti(L, 1, j);
      if (lua_istable(L, -1)) {
        luaL_argcheck(L, rows == (int)lua_rawlen(L, -1), 1,
                      "columns are not the same size");
        for (int k = 1; k <= rows; k++) {
          lua_rawgeti(L, -1, k);
          double v = luaL_checknumber(L, -1);
          lua_pop(L, 1);
          int pos = (j - 1) * rows + (k - 1);
          assert(pos >= 0 && pos < (rows * cols));
          m.data[pos] = v;
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

// Create a Ravi vector (ravi matrix with one column)
// Interface 1
//   arg1 - size
//   arg2 - initial value (optional)
// Interface 2
//   arg1 - Ravi number array to be converted
static int make_Ravi_vector(lua_State* L) {
  if (lua_isnumber(L, 1)) {
    int size = (int)lua_tointeger(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    lua_Number initv = 0.0;
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    alloc_Ravi_matrix(L, size, 1, initv);
    return 1;
  } else if (lua_gettop(L) == 1 && ravi_is_number_array(L, 1)) {
    int size = (int)lua_rawlen(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    set_Ravi_matrix_meta(L, size, 1);
    lua_pushvalue(L, -1);
    return 1;
  } else if (lua_istable(L, 1)) {
    int size = (int)lua_rawlen(L, 1);
    // allocate matrix of 1 column
    ravi_matrix_t vector = alloc_Ravi_matrix(L, (int)size, 1, 0.0);
    for (int x = 0; x < size; x++) {
      lua_rawgeti(L, 1, x + 1);
      double v = luaL_checknumber(L, -1);
      lua_pop(L, 1);
      vector.data[x] = v;
    }
    return 1;
  }
  luaL_argerror(L, 1, "Cannot create to vector from supplied arguments");
  return 0;
}
#endif

// obtain a particular element of a Lua matrix
// arg1 - Lua matrix
// arg2 - element position
static int Lua_vector_get(lua_State* L) {
  ravi_matrix_t vector = check_Lua_matrix(L, 1);
  int pos = (int)luaL_checkinteger(L, 2);
  luaL_argcheck(L, pos >= 1 && pos <= (vector.m * vector.n), 2,
                "read access out of bounds");
  lua_pushnumber(L, vector.data[pos - 1]);
  return 1;
}

// Update a particular element of a Lua matrix
// arg1 - Lua matrix
// arg2 - element position
// arg3 - value to set
static int Lua_vector_set(lua_State* L) {
  ravi_matrix_t vector = check_Lua_matrix(L, 1);
  int pos = (int)luaL_checkinteger(L, 2);
  luaL_argcheck(L, pos >= 1 && pos <= (vector.m * vector.n), 2,
                "write access out of bounds");
  lua_Number val = luaL_checknumber(L, 3);
  vector.data[pos - 1] = val;
  return 0;
}

// compute array length of a Lua matrix
static int Lua_vector_len(lua_State* L) {
  ravi_matrix_t vector = check_Lua_matrix(L, 1);
  lua_pushinteger(L, vector.m * vector.n);
  return 1;
}

// multipy two Lua matrices
static int Lua_matrix_mult(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A.n == B.m, 1, "matrices are not multiplicable");
  ravi_matrix_t C = alloc_Lua_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_multiply(&A, &B, &C, false, false, 1.0, 0.0);
  return 1;
}

static int Lua_matrix_add(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A.m == B.m && A.n == B.n, 1,
                "matrices are not the same size");
  ravi_matrix_t matrix = alloc_Lua_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  ravi_matrix_add(&matrix, &B, 1.0);
  return 1;
}

static int Lua_matrix_sub(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A.m == B.m && A.n == B.n, 1,
                "matrices are not the same size");
  ravi_matrix_t matrix = alloc_Lua_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  ravi_matrix_add(&matrix, &B, -1.0);
  return 1;
}

static int Lua_matrix_copy(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t matrix = alloc_Lua_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  return 1;
}

static int Lua_matrix_norm1(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_ONE));
  return 1;
}
static int Lua_matrix_normI(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_INFINITY));
  return 1;
}
static int Lua_matrix_normF(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_FROBENIUS));
  return 1;
}

static int Lua_matrix_lufactor(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  int32_t ipsize = A.m < A.n ? A.m : A.n;
  int* ipiv = (int*)alloca(sizeof(int) * ipsize);
  int info = ravi_matrix_lufactor(&A, ipsize, ipiv);
  if (info < 0)
    luaL_error(L, "Failed to factorize matrix: info %d", info);
  return 0;
}

static int Lua_matrix_solve(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t vector = check_Lua_matrix(L, 2);
  char alg = A.m == A.m ? 'L' : 'Q';
  if (lua_isstring(L, 3)) {
    alg = *lua_tostring(L, 3);
  }
  bool ok = false;
  luaL_argcheck(L, alg == 'L' || alg == 'Q' || alg == 'S', 3, "bad algorithm");
  ravi_matrix_t solution = alloc_Lua_matrix(L, vector.m, vector.n, 0.0);
  ravi_matrix_copy(&solution, &vector);
  ravi_matrix_t M = alloc_Lua_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&M, &A);
  switch (alg) {
    case 'L':
      ok = ravi_matrix_solve(&M, &solution);
      break;
    case 'Q':
      ok = ravi_matrix_lsq_solve(&M, &solution, 0.0, false);
      break;
    default:
      ok = ravi_matrix_lsq_solve(&M, &solution, 0.0, true);
      break;
  }
  lua_pop(L, 1);  // remove the matrix, leaving the vector which is the solution
  if (!ok) {
    luaL_error(L, "failed to solve");
    return 0;
  }
  return 1;
}

static int Lua_matrix_inverse(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t M = alloc_Lua_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&M, &A);
  if (!ravi_matrix_inverse(&M)) {
    luaL_error(L, "matrix is not invertible");
    return 0;
  }
  return 1;
}

static int Lua_matrix_transpose(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  ravi_matrix_t M = alloc_Lua_matrix(L, A.n, A.m, 0.0);
  ravi_matrix_transpose(&M, &A);
  return 1;
}

static int Lua_matrix_dimensions(lua_State* L) {
  ravi_matrix_t A = check_Lua_matrix(L, 1);
  lua_pushinteger(L, A.m);
  lua_pushinteger(L, A.n);
  return 2;
}

#if RAVI_ENABLED
/* multiply two number array matrices */
static int Ravi_matrix_mult(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A.n == B.m, 1, "matrices are not multiplicable");
  ravi_matrix_t C = alloc_Ravi_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_multiply(&A, &B, &C, false, false, 1.0, 0.0);
  return 1;
}

static int Ravi_matrix_add(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A.m == B.m && A.n == B.n, 1,
                "matrices are not the same size");
  ravi_matrix_t matrix = alloc_Ravi_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  ravi_matrix_add(&matrix, &B, 1.0);
  return 1;
}

static int Ravi_matrix_sub(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A.m == B.m && A.n == B.n, 1,
                "matrices are not the same size");
  ravi_matrix_t matrix = alloc_Ravi_matrix(L, A.m, B.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  ravi_matrix_add(&matrix, &B, -1.0);
  return 1;
}

static int Ravi_matrix_copy(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t matrix = alloc_Ravi_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&matrix, &A);
  return 1;
}

static int Ravi_matrix_norm1(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_ONE));
  return 1;
}
static int Ravi_matrix_normI(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_INFINITY));
  return 1;
}
static int Ravi_matrix_normF(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  lua_pushnumber(L, ravi_matrix_norm(&A, RAVI_MATRIX_NORM_FROBENIUS));
  return 1;
}

static int Ravi_matrix_lufactor(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  int32_t ipsize = A.m < A.n ? A.m : A.n;
  int* ipiv = (int*)alloca(sizeof(int) * ipsize);
  int info = ravi_matrix_lufactor(&A, ipsize, ipiv);
  if (info < 0)
    luaL_error(L, "Failed to factorize matrix: info %d", info);
  return 0;
}

static int Ravi_matrix_solve(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t vector = check_Ravi_matrix(L, 2);
  char alg = A.m == A.m ? 'L' : 'Q';
  if (lua_isstring(L, 3)) {
    alg = *lua_tostring(L, 3);
  }
  bool ok = false;
  luaL_argcheck(L, alg == 'L' || alg == 'Q' || alg == 'S', 3, "bad algorithm");
  ravi_matrix_t solution = alloc_Ravi_matrix(L, vector.m, vector.n, 0.0);
  ravi_matrix_copy(&solution, &vector);
  ravi_matrix_t M = alloc_Ravi_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&M, &A);
  switch (alg) {
    case 'L':
      ok = ravi_matrix_solve(&M, &solution);
      break;
    case 'Q':
      ok = ravi_matrix_lsq_solve(&M, &solution, 0.0, false);
      break;
    default:
      ok = ravi_matrix_lsq_solve(&M, &solution, 0.0, true);
      break;
  }
  lua_pop(L, 1);  // remove the matrix, leaving the vector which is the solution
  if (!ok) {
    luaL_error(L, "failed to solve");
    return 0;
  }
  return 1;
}

static int Ravi_matrix_inverse(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t M = alloc_Ravi_matrix(L, A.m, A.n, 0.0);
  ravi_matrix_copy(&M, &A);
  if (!ravi_matrix_inverse(&M)) {
    luaL_error(L, "matrix is not invertible");
    return 0;
  }
  return 1;
}

static int Ravi_matrix_transpose(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  ravi_matrix_t M = alloc_Ravi_matrix(L, A.n, A.m, 0.0);
  ravi_matrix_transpose(&M, &A);
  return 1;
}

static int Ravi_matrix_dimensions(lua_State* L) {
  ravi_matrix_t A = check_Ravi_matrix(L, 1);
  lua_pushinteger(L, A.m);
  lua_pushinteger(L, A.n);
  return 2;
}

#endif

static const struct luaL_Reg mylib[] = {{"vector", make_Lua_vector},
                                        {"matrix", make_Lua_matrix},
                                        {"copy", Lua_matrix_copy},
                                        {"solve", Lua_matrix_solve},
                                        {"inverse", Lua_matrix_inverse},
                                        {"transpose", Lua_matrix_transpose},
                                        {"dim", Lua_matrix_dimensions},
                                        {"norm1", Lua_matrix_norm1},
                                        {"normI", Lua_matrix_normI},
                                        {"normF", Lua_matrix_normF},
                                        {"lufactor", Lua_matrix_lufactor},
#if RAVI_ENABLED
                                        {"vectorR", make_Ravi_vector},
                                        {"matrixR", make_Ravi_matrix},
                                        {"copyR", Ravi_matrix_copy},
                                        {"solveR", Ravi_matrix_solve},
                                        {"inverseR", Ravi_matrix_inverse},
                                        {"transposeR", Ravi_matrix_transpose},
                                        {"dimR", Ravi_matrix_dimensions},
                                        {"norm1R", Ravi_matrix_norm1},
                                        {"normIR", Ravi_matrix_normI},
                                        {"normFR", Ravi_matrix_normF},
                                        {"lufactorR", Ravi_matrix_lufactor},
#endif
                                        {NULL, NULL}};

static const ravi_matrix_lua_api_t userdata_api = {
    test_Lua_matrix, check_Lua_matrix, alloc_Lua_matrix};

#if RAVI_ENABLED

static const ravi_matrix_lua_api_t ravi_array_api = {
    test_Ravi_matrix, check_Ravi_matrix, alloc_Ravi_matrix};

#endif

const ravi_matrix_lua_api_t* ravi_matrix_get_api(bool use_ravi_arrays) {
#if RAVI_ENABLED
  if (use_ravi_arrays)
    return &ravi_array_api;
#endif
  return &userdata_api;
}

// adds to an existing table
// table must be on stop of the stack
static void raviU_add_to_library(lua_State* L, const struct luaL_Reg* regs) {
  int i = 0;
  for (; regs[i].name != NULL; i++) {
    lua_pushcclosure(L, regs[i].func, 0);
    lua_setfield(L, -2, regs[i].name);
  }
}

// creates a table of functions which is a library in Lua
// terms. We use this as a portable way across 5.1 and 5.2 which
// follows the latest conventions - i.e. avoid polluting global
// namespace
void raviU_create_library(lua_State* L, const struct luaL_Reg* regs) {
  int count = 0;
  int i = 0;
  for (; regs[i].name != NULL; i++) {
    count++;
  }
  lua_createtable(L, 0, count);
  raviU_add_to_library(L, regs);
  // leave table on the stack
}

int luaopen_ravimatrix(lua_State* L) {
  fprintf(stderr, "Initializing RaviMatrix\n");
  raviL_newmetatable(L, Lua_matrix, Lua_matrix);
  lua_pushstring(L, "Lua matrix");
  lua_setfield(L, -2, "type");
  lua_pushcfunction(L, Lua_vector_set);
  lua_setfield(L, -2, "__newindex");
  lua_pushcfunction(L, Lua_vector_get);
  lua_setfield(L, -2, "__index");
  lua_pushcfunction(L, Lua_vector_len);
  lua_setfield(L, -2, "__len");
  lua_pushcfunction(L, Lua_matrix_mult);
  lua_setfield(L, -2, "__mul");
  lua_pushcfunction(L, Lua_matrix_add);
  lua_setfield(L, -2, "__add");
  lua_pushcfunction(L, Lua_matrix_sub);
  lua_setfield(L, -2, "__sub");
  lua_pop(L, 1);

#if RAVI_ENABLED
  raviL_newmetatable(L, Ravi_matrix, Ravi_matrix);
  lua_pushstring(L, "Ravi matrix");
  lua_setfield(L, -2, "type");
  lua_pushcfunction(L, Ravi_matrix_mult);
  lua_setfield(L, -2, "__mul");
  lua_pushcfunction(L, Ravi_matrix_add);
  lua_setfield(L, -2, "__add");
  lua_pushcfunction(L, Ravi_matrix_sub);
  lua_setfield(L, -2, "__sub");
  lua_pop(L, 1);
#endif

  raviU_create_library(L, mylib);
  fprintf(stdout, "RaviMatrix initialized successfully\n");
  return 1;
}
