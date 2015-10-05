find_path(RAVILUAUTILS_INCLUDE_DIR ravi_lua_utils.h
  PATHS
  c:/ravi/include
  ~/ravi/include
)

find_library(RAVILUAUTILS_LIBRARY
  NAMES raviluautils libraviluautils
  PATHS
  c:/ravi/bin
  ~/ravi/bin
)

set( RAVILUAUTILS_INCLUDE_DIRS "${RAVILUAUTILS_INCLUDE_DIR}" )
set( RAVILUAUTILS_LIBRARIES "${RAVILUAUTILS_LIBRARY}" )
