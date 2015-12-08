A Matrix Library for Ravi and Lua 5.3
=====================================

This project will create a matrix library for Ravi and lua 5.3


Pre-requisites
--------------
The following projects are pre-requisites:

* Ravi Dist - https://github.com/dibyendumajumdar/ravi-dist
* Ravi - https://github.com/dibyendumajumdar/ravi
* Ravi Lua Utils - https://github.com/dibyendumajumdar/ravi-luautils

Design 
------
The library is a thin wrapper around BLAS and LAPACK routines. The design has following salient features:

* The BLAS and LAPACK functions are wrapped in a C API that uses the following definition of a matrix structure::

    typedef struct ravi_matrix_t ravi_matrix_t;
    struct ravi_matrix_t {
      int32_t m; /* rows */
      int32_t n; /* columns */
      double data[1];
    };

* Above structure is designed to overlay a double array where [0] offset corresponds to the two 4-byte integer fields. This is important as Ravi number arrays are stored in double arrays.

* The C API itself is usable in C / C++ programs and has no dependency on Lua or Ravi.

* The C API is exposed as a set of function pointers to keep the namespace pollution to a minimum. The structure of the API can be seen in the header file `ravi_matrixlib.h <https://github.com/dibyendumajumdar/ravi-matrix/blob/master/include/ravi_matrixlib.h>`_.

* The Lua / Ravi APIs provide two implementations of the matrix object.

* A Lua compatible userdata based implementation is provided. 

* A different implemenatation is provided for Ravi that makes use of Ravi number arrays. This allows a Ravi number array to be converted to a matrix or vector - by simply a) attaching a metatable, and b) setting the row/column sizes in the slot [0] of the number array. This approach relies upon the fact that Lua programs use array indices >= 1 typically. Of course the slot [0] is accessible from Lua  or Ravi when using direct indexing; so a user of this library needs to be aware of the special use of this slot and not attempt to overwrite the value held there.

* In a Ravi environment both implementations will be available. In a Lua environment only the Lua compatible implementation will be available.

* The idea is that other libraries such as ravi-gsl can build upon the matrix and vector types provided in this library, thereby ensuring that the libraries work for standard Lua as well as providing a faster implementation when used in Ravi.

* The available Lua / Ravi API functions are as follows::

  TODO

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

