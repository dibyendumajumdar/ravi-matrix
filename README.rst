A Matrix Library for Ravi and Lua 5.3
=====================================

This project will create a matrix library for Ravi and lua 5.3


Pre-requisites
--------------
The following projects are pre-requisites:

* Ravi Dist - https://github.com/dibyendumajumdar/ravi-dist
* Ravi - https://github.com/dibyendumajumdar/ravi
* Ravi Lua Utils - https://github.com/dibyendumajumdar/ravi-luautils

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

