project(daqp_eigen)

find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} SHARED daqp.cpp)
add_library(daqp::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES
  VERSION ${daqp_VERSION}
  SOVERSION ${daqp_VERSION_MAJOR}.${daqp_VERSION_MINOR}
)
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen daqpstat)

if(TEST)
    add_executable(00_basic_qp tests/00_basic_qp.cpp)
    add_executable(01_static_func tests/01_static_func.cpp)
    add_executable(02_class tests/02_class.cpp)
    add_executable(03_update tests/03_update.cpp)
    add_executable(04_slack_sign tests/04_slack_sign.cpp)
    add_executable(05_warmstart tests/05_warmstart.cpp)
    add_executable(06_general_hessian tests/06_general_hessian.cpp)
    add_executable(07_symmetrize_hessian tests/07_symmetrize_hessian.cpp)
    add_executable(08_unconstrained_check tests/08_unconstrained_check.cpp)
    add_executable(09_soft_weights tests/09_soft_weights.cpp)
    add_executable(10_eq_reduction tests/10_eq_reduction.cpp)
    add_executable(11_warm_update_M tests/11_warm_update_M.cpp)
    add_executable(12_update_state tests/12_update_state.cpp)
    add_executable(13_warm_start_drift tests/13_warm_start_drift.cpp)
    add_executable(14_degenerate_lp_warm tests/14_degenerate_lp_warm.cpp)
    add_executable(15_degenerate_lp_skip tests/15_degenerate_lp_skip.cpp)
    add_executable(16_update_M tests/16_update_M.cpp)
    add_executable(17_cold_start_reduction tests/17_cold_start_reduction.cpp)
    add_executable(18_equality_warm_activation tests/18_equality_warm_activation.cpp)

    set(TARGETS
        00_basic_qp
        01_static_func
        02_class
        03_update
        04_slack_sign
        05_warmstart
        06_general_hessian
        07_symmetrize_hessian
        08_unconstrained_check
        09_soft_weights
        10_eq_reduction
        11_warm_update_M
        12_update_state
        13_warm_start_drift
        14_degenerate_lp_warm
        15_degenerate_lp_skip
        16_update_M
        17_cold_start_reduction
        18_equality_warm_activation
    )

    foreach(TARGET ${TARGETS})
        target_link_libraries(${TARGET} PRIVATE ${PROJECT_NAME})
        target_include_directories(${TARGET}
            PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test(NAME test_${TARGET} COMMAND ${TARGET})
    endforeach()
endif()

install(TARGETS ${PROJECT_NAME} EXPORT  ${PROJECT_NAME}Config
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES daqp.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/daqp")
