cmake_minimum_required(VERSION 3.1)

project(POrbitD)
#set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 99)
set(CMAKE_C_STANDARD 99)
set(p_version 0.2.1)

set(N_LIBS
    m)

set(COMMON_SOURCES
    "MatlabUtils.c"
    "rpoly.c"
    "MatlabUtils.h"
    "rpoly.h"
    "SAT_Const.h"
    "Position.h"
    "Position.c"
    "Mjday.h"
    "Mjday.c"
    "MeanObliquity.c"
    "MeanObliquity.h"
    "NutAngles.h"
    "NutAngles.c"
    "timediff.h"
    "timediff.c"
    "Frac.h"
    "Frac.c"
    "R_x.h"
    "R_x.c"
    "R_y.h"
    "R_y.c"
    "R_z.h"
    "R_z.c"
    "unit.h"
    "unit.c"
    "IERS.h"
    "IERS.c"
    "EqnEquinox.h"
    "EqnEquinox.c"
    "gmst.h"
    "gmst.c"
    "NutMatrix.h"
    "NutMatrix.c"
    "angl.h"
    "angl.c"
    "anglesdr.h"
    "anglesdr.c"
    "doubler.h"
    "doubler.c"
    "gast.h"
    "gast.c"
    "GHAMatrix.h"
    "GHAMatrix.c"
    "PrecMatrix.h"
    "PrecMatrix.c"
    "PoleMatrix.h"
    "PoleMatrix.c"
    "newtonnu.h"
    "newtonnu.c"
)
set(TEST_SOURCES
    ${COMMON_SOURCES}
    "TestUtils.h"
    "TestUtils.c"
    "MatlabUtilsTest.c"
    "MatlabUtilsTest.h"
    "ProjectTest.c")

set(CMAKE_C_OUTPUT_EXTENSION_REPLACE ON)

if(coverage)
    add_compile_options(-Wall -Wextra -pedantic -g -fprofile-arcs -ftest-coverage)
    set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}_test.dir)
    add_custom_target(coverage
        COMMAND gcov -o ${OBJECT_DIR} ${TEST_SOURCES})
    set(N_LIBS
        ${N_LIBS}
        gcov)
elseif(NOT coverage)
    add_compile_options(-Wall -Wextra -pedantic)
endif()

#add_executable(${PROJECT_NAME} ${COMMON_SOURCES})
add_executable(${PROJECT_NAME}_test ${TEST_SOURCES})
#target_link_libraries(${PROJECT_NAME} ${N_LIBS})
target_link_libraries(${PROJECT_NAME}_test ${N_LIBS})
#set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS TRUE)

if(docs)
    if(NOT p_version)
        execute_process(
            COMMAND git rev-parse HEAD
            OUTPUT_VARIABLE rev)
        string(STRIP ${rev} rev)
        set(p_version ${rev})
    endif()
    add_custom_target(docs
        COMMAND PROJECT_NUMBER=${p_version} doxygen Doxyfile)
endif()


