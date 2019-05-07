cmake_minimum_required(VERSION 3.1)

project(POrbitD)


set(N_LIBS
    m)

set(SOURCES
    "MatlabUtils.c"
    "MatlabUtilsTest.c"
    "rpoly.c"
    "MatlabUtils.h"
    "rpoly.h")

set(CMAKE_C_OUTPUT_EXTENSION_REPLACE ON)

if(coverage)
    add_compile_options(-Wall -Wextra -Werror -pedantic -g -fprofile-arcs -ftest-coverage)
    set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir)
    add_custom_target(coverage
        COMMAND gcov -o ${OBJECT_DIR} ${SOURCES})
    set(N_LIBS
        ${N_LIBS}
        gcov)
elseif(NOT coverage)
    add_compile_options(-Wall -Wextra -Werror -pedantic)
endif()

add_executable(${PROJECT_NAME} ${SOURCES})

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 99)
set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS TRUE)

if(docs)
    add_custom_target(docs
        COMMAND doxygen .)
endif()

target_link_libraries(${PROJECT_NAME} ${N_LIBS})
