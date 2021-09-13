cmake_minimum_required(VERSION 3.10)
project(DAQP VERSION 0.0.0)
add_subdirectory(src)
add_subdirectory(include)

add_library(daqpstat STATIC ${daqp_src} ${daqp_headers})
add_library(daqp SHARED ${daqp_src} ${daqp_headers})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC
set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES}")      # Include math
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra")
#target_link_libraries(daqpstat blas)
target_include_directories(daqpstat
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/daqp>")

target_include_directories(daqp
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/daqp>")

#add_executable(daqp src/test.c)
target_link_libraries(daqpstat daqp)

