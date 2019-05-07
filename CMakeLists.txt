cmake_minimum_required(VERSION 3.1)

project(POrbitD)

add_compile_options(-Wall -pedantic)
if(coverage)
    add_compile_options(-g -fprofile-arcs -ftest-coverage)
endif()
set(CMAKE_C_OUTPUT_EXTENSION_REPLACE ON)
set(SOURCES
    "MatlabUtils.c"
    "MatlabUtilsTest.c"
    "rpoly.c"
    "MatlabUtils.h"
    "rpoly.h")
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} m gcov)
set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 99)

set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir)
add_custom_target(gcov
    COMMAND gcov -o ${OBJECT_DIR} ${SOURCES})

set_target_properties(gcov PROPERTIES EXCLUDE_FROM_ALL TRUE)

