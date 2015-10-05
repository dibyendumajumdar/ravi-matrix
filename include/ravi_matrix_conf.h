#ifndef RAVI_MATRIX_CONF_H
#define RAVI_MATRIX_CONF_H

// This header file must be C compilable
// No C++ artifacts allowed outside of #ifdef __cplusplus
// Best to avoid

#ifdef _MSC_VER

#define RAVIMATRIX_DLLEXPORT __declspec(dllexport)
#define RAVIMATRIX_DLLIMPORT __declspec(dllimport)

#include <malloc.h>
#define alloca _alloca

#ifndef __cplusplus
#define inline __inline
#endif

#else

#define RAVIMATRIX_DLLEXPORT extern 
#define RAVIMATRIX_DLLIMPORT extern

#include <alloca.h>

#endif

// When compiling the library DLL this
// must be set, but when call the library
// from another program this must not be set
#ifdef RAVIMATRIX_IMPLEMENTATION
#define RAVIMATRIX_API RAVIMATRIX_DLLEXPORT
#else
#define RAVIMATRIX_API RAVIMATRIX_DLLIMPORT
#endif

#endif
