project(daqp_eigen)

find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)
find_package(daqp)

add_executable(00_basic_qp tests/00_basic_qp.cpp daqp.cpp)
add_executable(01_static_func tests/01_static_func.cpp daqp.cpp)
add_executable(02_class tests/02_class.cpp daqp.cpp)
add_executable(03_update tests/03_update.cpp daqp.cpp)

set(TARGETS
    00_basic_qp
    01_static_func
    02_class
    03_update
)

foreach(TARGET ${TARGETS})
    target_link_libraries(${TARGET} PUBLIC Eigen3::Eigen daqp)
    target_include_directories(${TARGET}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_test(NAME test_${TARGET} COMMAND ${TARGET})
endforeach()

install(FILES daqp.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/daqp")
