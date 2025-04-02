if(EIGEN)
    project(daqp_eigen)

    find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)
    find_package(daqp)

    include(CTest)

    add_executable(01_static_func tests/01_static_func.cpp)

    set(TARGETS 01_static_func)

    foreach(TARGET ${TARGETS})
        target_link_libraries(${TARGET} Eigen3::Eigen daqp)
        target_include_directories(${TARGET}
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test(NAME test_${TARGET} COMMAND ${TARGET})
    endforeach()
endif()
