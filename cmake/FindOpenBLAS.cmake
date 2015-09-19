find_path(OPENBLAS_INCLUDE_DIR openblas_config.h
  PATHS
  c:/ravi/include/openblas
  ~/ravi/include/openblas
)

find_library(OPENBLAS_LIBRARY
  NAMES openblas libopenblas
  PATHS
  c:/ravi/bin
  ~/ravi/bin
)

set( OPENBLAS_INCLUDE_DIRS "${OPENBLAS_INCLUDE_DIR}" )
set( OPENBLAS_LIBRARIES "${OPENBLAS_LIBRARY}" )
