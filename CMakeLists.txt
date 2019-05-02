cmake_minimum_required(VERSION 2.8)

project(POrbitD)
add_executable(${PROJECT_NAME} "MatlabUtils.c" "MatlabUtilsTest.c" "rpoly.c" "MatlabUtils.h" "rpoly.h")
target_link_libraries(${PROJECT_NAME} m)
