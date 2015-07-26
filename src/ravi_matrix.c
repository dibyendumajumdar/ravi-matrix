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

#include <ravimatrix/ravi_matrix.h>
#include <ravimatrix/matrixlib.h>

#define RAVI_ENABLED 1

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

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

// The normal Lua metatable functions in C use string
// keys - these are expensive as the key needs to be
// converted to Lua string, hash code computed etc.
// Following implementations are taken from a post in
// Lua mailing list (http://lua-users.org/lists/lua-l/2010-11/msg00151.html)
// They use lightuserdata instead of strings to speed
// things up
// meta_key is the key assigned to the meta
// table of the userdata
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
static const char *Lua_matrix = "Lua.matrix";
static const char *Ravi_matrix = "Ravi.matrix";

#define test_Lua_matrix(L, idx) ((matrix_t *)l_testudata(L, idx, Lua_matrix))
#define check_Lua_matrix(L, idx) ((matrix_t *)l_checkudata(L, idx, Lua_matrix))

// Allocate a Lua matrix as a userdata (compatible with Lua 5.3)
static inline matrix_t *alloc_Lua_matrix(lua_State *L, int m, int n) {
  assert(m >= 0 && n >= 0);
  assert(sizeof(matrix_t) == (sizeof(double) * 2));
  matrix_t *matrix =
      (matrix_t *)lua_newuserdata(L, sizeof(double) * (m * n + 1));
  matrix->m = m;
  matrix->n = n;
  l_getmetatable(L, Lua_matrix);
  lua_setmetatable(L, -2);
  return matrix;
}

#if RAVI_ENABLED
// Sets the metatable and other values for Ravi matrix
// Stack top must be a Ravi number array
static inline matrix_t *set_Ravi_matrix_meta(lua_State *L, int m, int n) {
  assert(m >= 0 && n >= 0);
  l_getmetatable(L, Ravi_matrix);
  lua_setmetatable(L, -2);
  double *start = ravi_get_number_array_rawdata(L, -1);
  matrix_t *matrix = (matrix_t *)start;
  assert(&matrix->data[0] == start + 1);
  matrix->m = m;
  matrix->n = n;
  return matrix;
}

static inline matrix_t *alloc_Ravi_matrix(lua_State *L, int m, int n,
                                          double initv) {
  ravi_create_number_array(L, m * n, initv);
  return set_Ravi_matrix_meta(L, m, n);
}

// Test that the argument is a Ravi_matrix
// arg_index is the position of argument on the stack
static matrix_t *test_Ravi_matrix(lua_State *L, int arg_index) {
  if (ravi_is_number_array(L, arg_index)) { // value is a Ravi array?
    if (lua_getmetatable(L, arg_index)) {   // does it have a metatable?
      l_getmetatable(L, Ravi_matrix);       // Get metatable for Ravi matrices
      if (lua_rawequal(L, -1, -2)) { // compare: does it have the correct mt?
        lua_pop(L, 2);               // remove both metatables
        return (matrix_t *)ravi_get_number_array_rawdata(L,
                                                         arg_index); // test ok
      }
    }
  }
  return NULL; /* to avoid warnings */
}

static inline matrix_t *check_Ravi_matrix(lua_State *L, int arg_index) {
  matrix_t *p = test_Ravi_matrix(L, arg_index);
  if (p == NULL)
    luaL_argerror(L, arg_index, Ravi_matrix);
  return p;
}
#endif

// Create a Lua vector (Lua matrix with one column)
// arg1 - size
// arg2 - initial value (optional)
static int make_Lua_vector(lua_State *L) {
  if (lua_isnumber(L, 1)) {
    int size = (int)lua_tointeger(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    lua_Number initv = 0.0; // initial value
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    // allocate matrix of 1 column
    matrix_t *vector = alloc_Lua_matrix(L, (int)size, 1);
    for (int i = 0; i < vector->m; i++)
      vector->data[i] = initv;
    return 1;
  } else if (lua_istable(L, 1)) {
    int size = (int)lua_objlen(L, 1);
    // allocate matrix of 1 column
    matrix_t *vector = alloc_Lua_matrix(L, (int)size, 1);
    for (int x = 0; x < size; x++) {
      lua_rawgeti(L, 1, x + 1);
      double v = luaL_checknumber(L, -1);
      lua_pop(L, 1);
      vector->data[x] = v;
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
//   arg1 - table of tables - outer table is columns
static int make_Lua_matrix(lua_State *L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");

    matrix_t *matrix = alloc_Lua_matrix(L, rows, cols);
    if (lua_istable(L, 3)) {
      int size = (int)lua_objlen(L, 3);
      luaL_argcheck(L, size == rows * cols, 3,
                    "Table supplied is not of required length");
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
    luaL_argcheck(L, cols > 0, 1, "Table must have at least one column table");
    // test first col
    lua_rawgeti(L, 1, 1);
    if (lua_istable(L, -1)) {
      // all rows are expected to be same length
      rows = (int)lua_objlen(L, -1);
    } else {
      luaL_argerror(L, 1, "Expecting (table) array of arrays");
    }
    lua_pop(L, 1); // pop the test value
    m = alloc_Lua_matrix(L, rows, cols);
    for (int j = 1; j <= cols; j++) {
      lua_rawgeti(L, 1, j);
      if (lua_istable(L, -1)) {
        luaL_argcheck(L, rows == (int)lua_objlen(L, -1), 1,
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
static int make_Ravi_matrix(lua_State *L) {
  int top = lua_gettop(L);
  if (top >= 2 && lua_isnumber(L, 1) && lua_isnumber(L, 2)) {
    int rows = (int)lua_tointeger(L, 1);
    int cols = (int)lua_tointeger(L, 2);
    luaL_argcheck(L, rows >= 0, 1, "rows must be >= 0");
    luaL_argcheck(L, cols >= 0, 2, "cols must be >= 0");
    if (lua_gettop(L) == 3 && ravi_is_number_array(L, 3)) {
      int size = (int)lua_objlen(L, 3);
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
    matrix_t *m = NULL;
    // we expect an array for each column
    int cols = (int)lua_objlen(L, 1);
    int rows = -1; // don't know how many rows yet
    luaL_argcheck(L, cols > 0, 1, "Table must have at least one column table");
    // test first col
    lua_rawgeti(L, 1, 1);
    if (lua_istable(L, -1)) {
      // all rows are expected to be same length
      rows = (int)lua_objlen(L, -1);
    } else {
      luaL_argerror(L, 1, "Expecting (table) array of arrays");
    }
    lua_pop(L, 1); // pop the test value
    m = alloc_Ravi_matrix(L, rows, cols, 0.0);
    for (int j = 1; j <= cols; j++) {
      lua_rawgeti(L, 1, j);
      if (lua_istable(L, -1)) {
        luaL_argcheck(L, rows == (int)lua_objlen(L, -1), 1,
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

// Create a Ravi vector (ravi matrix with one column)
// Interface 1
//   arg1 - size
//   arg2 - initial value (optional)
// Interface 2
//   arg1 - Ravi number array to be converted
static int make_Ravi_vector(lua_State *L) {
  if (lua_isnumber(L, 1)) {
    int size = (int)lua_tointeger(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    lua_Number initv = 0.0;
    if (lua_isnumber(L, 2))
      initv = lua_tonumber(L, 2);
    alloc_Ravi_matrix(L, size, 1, initv);
    return 1;
  } else if (lua_gettop(L) == 1 && ravi_is_number_array(L, 1)) {
    int size = (int)lua_objlen(L, 1);
    luaL_argcheck(L, size >= 0, 1, "size must be >= 0");
    set_Ravi_matrix_meta(L, size, 1);
    lua_pushvalue(L, -1);
    return 1;
  } else if (lua_istable(L, 1)) {
    int size = (int)lua_objlen(L, 1);
    // allocate matrix of 1 column
    matrix_t *vector = alloc_Ravi_matrix(L, (int)size, 1, 0.0);
    for (int x = 0; x < size; x++) {
      lua_rawgeti(L, 1, x + 1);
      double v = luaL_checknumber(L, -1);
      lua_pop(L, 1);
      vector->data[x] = v;
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
static int Lua_vector_get(lua_State *L) {
  matrix_t *vector = check_Lua_matrix(L, 1);
  int pos = (int)luaL_checkinteger(L, 2);
  luaL_argcheck(L, pos >= 1 && pos <= (vector->m * vector->n), 2,
                "access out of bounds");
  lua_pushnumber(L, vector->data[pos - 1]);
  return 1;
}

// compute array length of a Lua matrix
static int Lua_vector_len(lua_State *L) {
  matrix_t *vector = check_Lua_matrix(L, 1);
  lua_pushinteger(L, vector->m * vector->n);
  return 1;
}

// multipy two Lua matrices
static int Lua_matrix_mult(lua_State *L) {
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A->n == B->m, 1, "matrices are not multiplicable");
  matrix_t *matrix = alloc_Lua_matrix(L, A->m, B->n);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  if (!ops->multiply(matrix, A, B, false, false)) {
    luaL_error(L, "matrix multiplication failed");
    return 0;
  }
  return 1;
}

static int Lua_matrix_add(lua_State *L) {
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A->m == B->m && A->n == B->n, 1,
                "matrices are not the same size");
  matrix_t *matrix = alloc_Lua_matrix(L, A->m, B->n);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  ops->add(matrix, B);
  return 1;
}

static int Lua_matrix_sub(lua_State *L) {
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *B = check_Lua_matrix(L, 2);
  luaL_argcheck(L, A->m == B->m && A->n == B->n, 1,
                "matrices are not the same size");
  matrix_t *matrix = alloc_Lua_matrix(L, A->m, B->n);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  ops->sub(matrix, B);
  return 1;
}

static int Lua_matrix_copy(lua_State *L) {
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *matrix = alloc_Lua_matrix(L, A->m, A->n);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  return 1;
}

static int Lua_matrix_solve(lua_State *L) {
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *vector = check_Lua_matrix(L, 2);
  char alg = A->m == A->m ? 'L' : 'Q';
  if (lua_isstring(L, 3)) {
    alg = *lua_tostring(L, 3);
  }
  bool ok = false;
  luaL_argcheck(L, alg == 'L' || alg == 'Q' || alg == 'S', 3, "bad algorithm");
  matrix_t *solution = alloc_Lua_matrix(L, vector->m, vector->n);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(solution, vector);
  matrix_t *M = alloc_Lua_matrix(L, A->m, A->n);
  ops->copy(M, A);
  switch (alg) {
  case 'L':
    ok = ops->solve(M, solution);
    break;
  case 'Q':
    ok = ops->lsq_solve(M, solution, 0.0, false);
    break;
  default:
    ok = ops->lsq_solve(M, solution, 0.0, true);
    break;
  }
  lua_pop(L, 1); // remove the matrix, leaving the vector which is the solution
  if (!ok) {
    luaL_error(L, "failed to solve");
    return 0;
  }
  return 1;
}

static int Lua_matrix_inverse(lua_State *L) {
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  matrix_t *A = check_Lua_matrix(L, 1);
  matrix_t *M = alloc_Lua_matrix(L, A->m, A->n);
  ops->copy(M, A);
  if (!ops->inverse(M)) {
    luaL_error(L, "matrix is not invertible");
    return 0;
  }
  return 1;
}

#if RAVI_ENABLED
/* multiply two number array matrices */
static int Ravi_matrix_mult(lua_State *L) {
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A->n == B->m, 1, "matrices are not multiplicable");
  matrix_t *matrix = alloc_Ravi_matrix(L, A->m, B->n, 0.0);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  if (!ops->multiply(matrix, A, B, false, false)) {
    luaL_error(L, "matrix multiplication failed");
    return 0;
  }
  return 1;
}

static int Ravi_matrix_add(lua_State *L) {
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A->m == B->m && A->n == B->n, 1,
                "matrices are not the same size");
  matrix_t *matrix = alloc_Ravi_matrix(L, A->m, B->n, 0.0);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  ops->add(matrix, B);
  return 1;
}

static int Ravi_matrix_sub(lua_State *L) {
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *B = check_Ravi_matrix(L, 2);
  luaL_argcheck(L, A->m == B->m && A->n == B->n, 1,
                "matrices are not the same size");
  matrix_t *matrix = alloc_Ravi_matrix(L, A->m, B->n, 0.0);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  ops->sub(matrix, B);
  return 1;
}

static int Ravi_matrix_copy(lua_State *L) {
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *matrix = alloc_Ravi_matrix(L, A->m, A->n, 0.0);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(matrix, A);
  return 1;
}

static int Ravi_matrix_solve(lua_State *L) {
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *vector = check_Ravi_matrix(L, 2);
  char alg = A->m == A->m ? 'L' : 'Q';
  if (lua_isstring(L, 3)) {
    alg = *lua_tostring(L, 3);
  }
  bool ok = false;
  luaL_argcheck(L, alg == 'L' || alg == 'Q' || alg == 'S', 3, "bad algorithm");
  matrix_t *solution = alloc_Ravi_matrix(L, vector->m, vector->n, 0.0);
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  ops->copy(solution, vector);
  matrix_t *M = alloc_Ravi_matrix(L, A->m, A->n, 0.0);
  ops->copy(M, A);
  switch (alg) {
  case 'L':
    ok = ops->solve(M, solution);
    break;
  case 'Q':
    ok = ops->lsq_solve(M, solution, 0.0, false);
    break;
  default:
    ok = ops->lsq_solve(M, solution, 0.0, true);
    break;
  }
  lua_pop(L, 1); // remove the matrix, leaving the vector which is the solution
  if (!ok) {
    luaL_error(L, "failed to solve");
    return 0;
  }
  return 1;
}

static int Ravi_matrix_inverse(lua_State *L) {
  const matrix_ops_t *ops = ravi_matrix_get_implementation();
  matrix_t *A = check_Ravi_matrix(L, 1);
  matrix_t *M = alloc_Ravi_matrix(L, A->m, A->n, 0.0);
  ops->copy(M, A);
  if (!ops->inverse(M)) {
    luaL_error(L, "matrix is not invertible");
    return 0;
  }
  return 1;
}

#endif

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

static const struct luaL_Reg mylib[] = {{"vector", make_Lua_vector},
                                        {"matrix", make_Lua_matrix},
                                        {"copy", Lua_matrix_copy},
                                        {"solve", Lua_matrix_solve},
                                        {"inv", Lua_matrix_inverse},

#if RAVI_ENABLED
                                        {"vectorx", make_Ravi_vector},
                                        {"matrixx", make_Ravi_matrix},
                                        {"copyx", Ravi_matrix_copy},
                                        {"invx", Ravi_matrix_inverse},
#endif
                                        {NULL, NULL}};

int luaopen_ravimatrix(lua_State *L) {
  fprintf(stderr, "Initializing RaviMatrix\n");
  l_newmetatable(L, Lua_matrix);
  lua_pushstring(L, "Lua matrix");
  lua_setfield(L, -2, "type");
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
  l_newmetatable(L, Ravi_matrix);
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

  create_library(L, mylib);
  fprintf(stdout, "RaviMatrix initialized successfully\n");
  return 1;
}
