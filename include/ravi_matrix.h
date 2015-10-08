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

#ifndef RAVI_MATRIX_H_
#define RAVI_MATRIX_H_

#include <ravi_matrix_conf.h>
#include <ravi_matrixlib.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

typedef struct ravi_matrix_lua_api_t ravi_matrix_lua_api_t;
struct ravi_matrix_lua_api_t {
  ravi_matrix_t *(*test_ismatrix)(lua_State *L, int idx);
  ravi_matrix_t *(*check_ismatrix)(lua_State *L, int idx);
  ravi_matrix_t *(*alloc_matrix)(lua_State *L, int m, int n, double initv);
};

RAVIMATRIX_API const ravi_matrix_lua_api_t *ravi_matrix_get_api(bool use_ravi_array);

RAVIMATRIX_API int luaopen_ravimatrix(lua_State *L);

#ifdef __cplusplus
}
#endif

#endif
