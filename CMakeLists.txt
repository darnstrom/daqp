cmake_minimum_required(VERSION 3.10)
project(daqp VERSION 0.8.4)


include(GNUInstallDirs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/codegen)
add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(codegen)


if(TEST)
  include(CTest)
  enable_testing()
elseif(BENCHMARK)
  include(CTest)
  enable_testing()
endif()

option (PROFILING "Enable profiling" ON)
if(PROFILING)
  add_compile_definitions(PROFILING)
endif()

option (SOFT_WEIGHTS "Use different weights/bounds for soft constraints slacks" OFF)
if(SOFT_WEIGHTS)
  add_compile_definitions(SOFT_WEIGHTS)
endif()

option (FASTER_MATH "Enable some faster math" ON)
if(FASTER_MATH)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fassociative-math -fno-signed-zeros -fno-trapping-math")
endif()

option (BENCHMARK "Enable performance benchmark tests (compares current version with git tag)" OFF)
if(BENCHMARK)
  if(NOT JULIA)
    message(WARNING "BENCHMARK requires Julia interface (-DJULIA=ON), but JULIA is not enabled")
  endif()
endif()


add_library(daqpstat STATIC ${daqp_src} ${daqp_headers} ${daqp_codegen_src} ${daqp_codegen_headers})
target_include_directories(daqpstat PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(daqp SHARED ${daqp_src} ${daqp_headers} ${daqp_codegen_src} ${daqp_codegen_headers})
target_include_directories(daqp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES}") # math for C api
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra")

target_link_libraries(daqpstat daqp)

# Install
install(TARGETS daqpstat EXPORT  ${PROJECT_NAME}Config
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS daqp EXPORT  ${PROJECT_NAME}Config
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES ${daqp_headers} ${daqp_codegen_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/daqp")

install(EXPORT ${PROJECT_NAME}Config DESTINATION cmake)

if(EIGEN)
  add_subdirectory(interfaces/daqp-eigen)
endif()


if(MATLAB)
  find_package(Matlab REQUIRED MAIN_PROGRAM)
  set(MATLAB_INTERFACE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/interfaces/daqp-matlab")
  matlab_add_mex(NAME daqpmex SRC ${MATLAB_INTERFACE_DIR}/daqpmex.c LINK_TO daqpstat)
  set_target_properties(daqpmex PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${MATLAB_INTERFACE_DIR})
  matlab_add_unit_test(NAME MATLAB UNITTEST_FILE ${MATLAB_INTERFACE_DIR}/test/core_test.m ADDITIONAL_PATH ${MATLAB_INTERFACE_DIR} ${MATLAB_INTERFACE_DIR}/utils)
endif()

if(JULIA)
    set(JULIA_INTERFACE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/interfaces/daqp-julia")
    set(JULIA_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/interfaces/daqp-julia")

    add_custom_target(sync_julia_interface ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${JULIA_INTERFACE_DIR}" "${JULIA_BINARY_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:daqp> "${JULIA_BINARY_DIR}"
        DEPENDS daqp
        COMMENT "Syncing Julia files and library..."
    )

    add_dependencies(sync_julia_interface daqp)

    if(TEST)
        add_test(NAME Julia COMMAND sh -c "julia --color=yes --project=${CMAKE_CURRENT_BINARY_DIR}/interfaces/daqp-julia -e \"using Pkg; Pkg.test()\"")
    endif()
    
    if(BENCHMARK)
        # Add benchmark comparison test (git-based)
        add_test(
            NAME JuliaBenchmark
            COMMAND bash "${JULIA_INTERFACE_DIR}/test/benchmark_comparison_git.sh"
                "${JULIA_BINARY_DIR}"
                "v0.8.0"
                "medium"
                "5"
                "${CMAKE_CURRENT_BINARY_DIR}/benchmark_results"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
        set_tests_properties(JuliaBenchmark PROPERTIES
            TIMEOUT 600
            LABELS "benchmark"
        )
    endif()
endif()

if(PYTHON)
    set(PYTHON_INTERFACE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/interfaces/daqp-python")
    set(PYTHON_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/interfaces/daqp-python")
    set(PYTHON_VENV_DIR "${PYTHON_BINARY_DIR}/venv")

    add_custom_target(sync_python_interface ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${PYTHON_INTERFACE_DIR}" "${PYTHON_BINARY_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/src" "${PYTHON_BINARY_DIR}/csrc/src"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/include" "${PYTHON_BINARY_DIR}/csrc/include"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/codegen" "${PYTHON_BINARY_DIR}/csrc/codegen"
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" "${PYTHON_BINARY_DIR}/csrc/"
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" "${PYTHON_BINARY_DIR}/"
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:daqp> "${PYTHON_BINARY_DIR}"
        DEPENDS daqp
        COMMENT "Syncing Python files and library..."
    )

    add_dependencies(sync_python_interface daqp)

    add_test(NAME Python COMMAND sh -c "cd ${PYTHON_BINARY_DIR} && python -m venv venv && ${PYTHON_VENV_DIR}/bin/pip install -e . && ${PYTHON_VENV_DIR}/bin/python -m unittest discover -s test -p '*_test.py' -v")
endif()
