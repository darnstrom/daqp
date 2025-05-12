project(daqp_eigen)

find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} SHARED daqp.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen daqpstat)

add_executable(00_basic_qp tests/00_basic_qp.cpp)
add_executable(01_static_func tests/01_static_func.cpp)
add_executable(02_class tests/02_class.cpp)
add_executable(03_update tests/03_update.cpp)

set(TARGETS
    00_basic_qp
    01_static_func
    02_class
    03_update
)

foreach(TARGET ${TARGETS})
    target_link_libraries(${TARGET} PRIVATE ${PROJECT_NAME})
    target_include_directories(${TARGET}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_test(NAME test_${TARGET} COMMAND ${TARGET})
endforeach()

install(TARGETS ${PROJECT_NAME} EXPORT  ${PROJECT_NAME}Config
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES daqp.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/daqp")
