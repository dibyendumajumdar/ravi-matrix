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

// This header file must be C compilable
// No C++ artifacts allowed outside of #ifdef __cplusplus
// Best to avoid

#ifdef _MSC_VER

#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)

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

// When compiling the library DLL this
// must be set, but when call the library
// from another program this must not be set
#ifdef RAVIMATRIX_IMPLEMENTATION
#define API DLLEXPORT
#else
#define API DLLIMPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

DLLEXPORT int luaopen_ravimatrix(lua_State *L);

#ifdef __cplusplus
}
#endif

#endif