set(
    daqp_headers
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/api.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/auxiliary.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/bnb.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/codegen.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/constants.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/daqp.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/daqp_prox.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/factorization.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/types.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/daqp/utils.h"
)

set(daqp_headers "${daqp_headers}" PARENT_SCOPE)
