#ifndef _RAVI_CONF_H
#define _RAVI_CONF_H

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

#ifdef _WIN32
// When compiling the library DLL this
// must be set, but when call the library
// from another program this must not be set
#ifdef RAVIMATRIX_IMPLEMENTATION
#define API DLLEXPORT
#else
#define API DLLIMPORT
#endif

#else

#define API extern

#endif

#endif
