project(daqp_eigen)

include(CTest)

file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})

find_package (Eigen3 3.4.0 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

find_package (daqp)
target_link_libraries(${PROJECT_NAME} daqp)

target_include_directories(
    ${PROJECT_NAME}
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
