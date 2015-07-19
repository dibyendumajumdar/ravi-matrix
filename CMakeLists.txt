cmake_minimum_required(VERSION 2.8)
project(RaviMatrix)

enable_language(C)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

#LAPACK
if (NOT MSVC)
  find_package(LAPACK REQUIRED)
endif()

find_package(RAVI REQUIRED)

message( STATUS "LAPACK_LINKER_FLAGS   : " ${LAPACK_LINKER_FLAGS} )
message( STATUS "LAPACK_LIBRARIES      : " ${LAPACK_LIBRARIES} )
message( STATUS "BLAS_LIBRARIES        : " ${BLAS_LIBRARIES} )
message( STATUS "BLAS_LINKER_FLAGS     : " ${BLAS_LINKER_FLAGS} )
message( STATUS "LUA_LIBRARIES         : " ${LUA_LIBRARIES} )
message( STATUS "LUA_INCLUDE_DIRS      : " ${LUA_INCLUDE_DIRS} )

if (MSVC)
  set(CMAKE_C_FLAGS_DEBUG "/Od /D_DEBUG /MDd /Zi /RTC1 /EHsc")
  set(CMAKE_C_FLAGS_RELEASE "/DNDEBUG /O2 /MD /EHsc")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "/DNDEBUG /O2 /MD /Zi /EHsc")
  set(CMAKE_C_FLAGS_MINSIZEREL "/DNDEBUG /O2 /MD /EHsc")    

  set(CMAKE_CXX_FLAGS "/EHsc")
  set(CMAKE_CXX_FLAGS_DEBUG "/Od /D_DEBUG /MDd /Zi /RTC1 /EHsc")
  set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /O2 /MD  /EHsc")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/DNDEBUG /O2 /MD /Zi /EHsc")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "/DNDEBUG /O2 /MD /EHsc")  
  
  SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "/DEBUG /STACK:10000000")
  SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "/RELEASE /INCREMENTAL:NO /STACK:10000000")
  SET(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "/DEBUG /RELEASE /INCREMENTAL:NO /STACK:10000000")
  SET(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "/RELEASE /INCREMENTAL:NO /STACK:10000000")

  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TCMALLOC_LINK_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
  SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
  SET(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")
  SET(CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}")

  # Disable signed/unsigned int conversion warnings.
  add_definitions("/wd4018")
  # Disable warning about using struct/class for the same symobl.
  add_definitions("/wd4099")
  # Disable performance warning about int-to-bool conversion.
  add_definitions("/wd4800")
  # Disable signed/unsigned int conversion warnings.
  add_definitions("/wd4018")
  # Disable warning about the insecurity of using "std::copy".
  add_definitions("/wd4996")
  add_definitions("/wd4291")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS_DEBUG "-Werror -std=c99 -O0 -g3 -fsanitize=address -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
  set(CMAKE_C_FLAGS_RELEASE "-Werror -std=c99 -O3 -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-Werror -std=c99 -O1 -g3 -fsanitize=address -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
  set(CMAKE_CXX_FLAGS_RELEASE "-Werror -fno-rtti -O3 -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
  set(CMAKE_CXX_FLAGS_DEBUG "-Werror -fno-rtti -O0 -fsanitize=address -g3 -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Werror -fno-rtti -O1 -g3 -fsanitize=address -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_LINUX")
endif()

if (APPLE)
  set(CMAKE_C_FLAGS "-std=c99 -O1 -g3 -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_C_FLAGS_DEBUG "-std=c99 -O0 -g3 -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_C_FLAGS_RELEASE "-std=c99 -O3 -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-std=c99 -O1 -g3 -Wall -Wextra -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_CXX_FLAGS "-O1 -g3 -Wall -fno-rtti -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_CXX_FLAGS_RELEASE "-fno-rtti -O3 -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_CXX_FLAGS_DEBUG "-fno-rtti -O0 -g3 -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-fno-rtti -O1 -g3 -Wall -Wno-sign-compare -std=c++11 -fno-exceptions -DLUA_COMPAT_5_2 -DLUA_USE_MACOSX")
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")
add_definitions(-DLUA_COMPAT_MODULE)

if (MSVC)
  # LAPACK defines  
  add_definitions(-DADD_)
  add_definitions(-DHAVE_LAPACK_CONFIG_H)
  add_definitions(-DLAPACK_COMPLEX_STRUCTURE)
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${LUA_INCLUDE_DIRS}")

if (APPLE)
  include_directories("/System/Library/Frameworks")
endif()

if (MSVC)
  include_directories("c:/libs")
  set(LAPACK_LIBRARIES "liblapack.lib;libblas.lib")
  link_directories("c:/libs")
elseif(APPLE)
else()
endif()

set (RAVIMATRIX_SRCS src/matrix.c src/ravi_matrix.c)
    
if (MSVC)
  source_group("Source Files" FILES ${RAVIMATRIX_SRCS})
endif(MSVC)

file(GLOB RAVIMATRIX_HEADERS "${PROJECT_SOURCE_DIR}/include/ravimatrix/*.h")

if (MSVC)
  source_group("Headers" FILES ${RAVIMATRIX_HEADERS})
endif()

#Main library
add_library(ravimatrix SHARED 
  ${RAVIMATRIX_HEADERS} ${LUA_HEADERS} 
  ${RAVIMATRIX_SRCS} )
# RAVIMATRIX_IMPLEMENTATION is set when compiling shared library so that 
# appropriate export declarations can be set 
set_target_properties(ravimatrix 
  PROPERTIES DEFINE_SYMBOL "RAVIMATRIX_IMPLEMENTATION")
target_link_libraries(ravimatrix ${LAPACK_LIBRARIES} ${LUA_LIBRARIES}) 

