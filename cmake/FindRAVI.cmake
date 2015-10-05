find_path(LUA_INCLUDE_DIR lua.h
  PATHS
  c:/ravi/include
  ~/ravi/include
)

find_library(LUA_LIBRARY
  NAMES ravi libravi
  PATHS
  c:/ravi/bin
  ~/ravi/bin
)

set( LUA_INCLUDE_DIRS "${LUA_INCLUDE_DIR}" )
set( LUA_LIBRARIES "${LUA_LIBRARY}" )
