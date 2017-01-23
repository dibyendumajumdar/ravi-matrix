A Matrix Library for Ravi and Lua 5.3
=====================================

This project will create a matrix library for Ravi and lua 5.3


Pre-requisites
--------------
The following projects are pre-requisites:

* Ravi Dist - https://github.com/dibyendumajumdar/ravi-dist
* Ravi - https://github.com/dibyendumajumdar/ravi

Design 
------
The library is a thin wrapper around BLAS and LAPACK routines. The design has following salient features:

* The BLAS and LAPACK functions are wrapped in a C API that uses the following definition of a matrix structure::

    typedef struct ravi_matrix_t ravi_matrix_t;
    struct ravi_matrix_t {
      int32_t m; /* rows */
      int32_t n; /* columns */
      double *data;
    };

* The C API itself is usable in C / C++ programs and has no dependency on Lua or Ravi.

* The C API is exposed as a set of functions . The contents of the API can be seen in the header file `ravi_matrixlib.h <https://github.com/dibyendumajumdar/ravi-matrix/blob/master/include/ravi_matrixlib.h>`_.

* The Lua / Ravi APIs provide two implementations of the matrix object.

* A Lua compatible userdata based implementation is provided. 

* A different implementation is provided for Ravi that makes use of Ravi number arrays. This allows a Ravi number array to be converted to a matrix or vector - by simply a) attaching a metatable, and b) setting the row/column sizes in the slot [0] of the number array. This approach relies upon the fact that Lua programs use array indices >= 1 typically. Of course the slot [0] is accessible from Lua  or Ravi when using direct indexing; so a user of this library needs to be aware of the special use of this slot and not attempt to overwrite the value held there.

* In a Ravi environment both implementations will be available. In a Lua environment only the Lua compatible implementation will be available.

* The idea is that other libraries can build upon the matrix and vector types provided in this library, thereby ensuring that the libraries work for standard Lua as well as providing a faster implementation when used in Ravi.

* The available Lua / Ravi API functions are as follows.

Lua API
-------
The Lua compatible API can be accessed via the 'ravimatrix' module::

  local matrix = require 'ravimatrix'
  
The available functions are:

*matrix.vector(n: integer [,initialvalue: number])*
  Creates a column vector of size 'n', and initializes the vector elements to 'initialvalue' if supplied else to 0.

::

    local v = matrix.vector(2, 2.0)
    assert(#v == 2)
    assert(v[1] == 2.0)
    assert(v[2] == 2.0)

*matrix.vector(t: table)*
  Creates a column vector from supplied table t.
  
::
    
    local bx = matrix.vector { 10, 7, 43 }

*matrix.matrix(m: integer, n: integer [, initialvalue: number])*
  Creates a matrix of 'm' rows, 'n' columns; and the elements are initialized to 'initialvalue' if supplied else to 0.

*matrix.matrix(m: integer, n: integer, t: table)*
  Creates a matrix of 'm' rows, 'n' columns and initializes using the data in table 't'; the table is interpreted in column major order.

*matrix.matrix(t: table)*
  Creates a matrix from supplied table, which must be a table of tables; the each inner table representing a column.

::

    -- 3x3 matrix A
    -- 76 25 11
    -- 27 89 51
    -- 18 60 32

    local A = matrix.matrix { {76,27,18}, {25,89,60}, {11,51,32} }

Ravi API
--------
The Ravi API is not compatible with Lua as it requires the use of number arrays. A number of the Lua API methods have counterparts in Ravi API.

*matrix.vectorR(n: integer [,initialvalue: number])*
  Creates a column vector of size 'n', and initializes the vector elements to 'initialvalue' if supplied else to 0.

::

    local v: number[] = matrix.vectorR (2, 2.0)
    assert(#v == 2)
    assert(v[1] == 2.0)
    assert(v[2] == 2.0)

*matrix.vectorR(t: table)*
  Creates a column vector from supplied table t.
  
::
    
    local bx: number[] = matrix.vectorR { 10, 7, 43 }

*matrix.vectorR(t: number[])*
  Converts the supplied number[] to a vector by attaching a metatable.
  
::
    
    local data: number[] = { 10, 7, 43 }
    local bx: number[] = matrix.vectorR (data) 


*matrix.matrixR(m: integer, n: integer [, initialvalue: number])*
  Creates a matrix of 'm' rows, 'n' columns; and the elements are initialized to 'initialvalue' if supplied else to 0.

*matrix.matrixR(m: integer, n: integer, t: number[])*
  Converts the supplied number[] object to a matrix of 'm' rows, 'n' columns by attaching a metatable.

*matrix.matrixR(t: table)*
  Creates a matrix from supplied table, which must be a table of tables; the each inner table representing a column.

::

    -- 3x3 matrix A
    -- 76 25 11
    -- 27 89 51
    -- 18 60 32

    local A: number[] = matrix.matrixR { {76,27,18}, {25,89,60}, {11,51,32} }

Building on Windows
-------------------

::

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=c:\ravi -G "Visual Studio 14 Win64" -DCMAKE_BUILD_TYPE=Release ..

Then open is Visual Studio 2015 and do the build from there.

Building on UNIX or MAC OSX
---------------------------

::

    mkdir build
    cd build
    cmake  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/ravi ..
    make 
    make install

