cmake_minimum_required(VERSION 3.1)

project(POrbitD)
set(CMAKE_C_STANDARD 99)
set(p_version 1.0.2)

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
    "rv2coe.h"
    "rv2coe.c"
    "gibbs.h"
    "gibbs.c"
    "hgibbs.h"
    "hgibbs.c"
    "lambert_gooding.h"
    "lambert_gooding.c"
)

set(DOUBLER_SOURCES
    ${COMMON_SOURCES}
    "anglesdr.h"
    "anglesdr.c"
    "doubler.h"
    "doubler.c")

set(GAUSS_SOURCES
    ${COMMON_SOURCES}
    "anglesg.h"
    "anglesg.c")
set(TEST_SOURCES
    ${COMMON_SOURCES}
    "anglesdr.h"
    "anglesdr.c"
    "doubler.h"
    "doubler.c"
    "anglesg.h"
    "anglesg.c"
    "TestUtils.h"
    "TestUtils.c"
    "MatlabUtilsTest.c"
    "MatlabUtilsTest.h"
    "ProjectTest.c"
    )

set(CMAKE_C_OUTPUT_EXTENSION_REPLACE ON)

if(coverage)
    add_compile_options(-g -fprofile-arcs -ftest-coverage)
    set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/)
    add_custom_target(coverage
        COMMAND gcov -o ${OBJECT_DIR}/test-orbit.dir ${TEST_SOURCES}
        COMMAND gcov -o ${OBJECT_DIR}/example1.dir ${DOUBLER_SOURCES} "example1.c"
        COMMAND gcov -o ${OBJECT_DIR}/example2.dir ${DOUBLER_SOURCES} "example2.c"
        COMMAND gcov -o ${OBJECT_DIR}/example3.dir ${DOUBLER_SOURCES} "example3.c"
        COMMAND gcov -o ${OBJECT_DIR}/example5.dir ${GAUSS_SOURCES} "example5.c"
        COMMAND gcov -o ${OBJECT_DIR}/example6.dir ${DOUBLER_SOURCES} "example6.c"
        COMMAND gcov -o ${OBJECT_DIR}/example7.dir ${DOUBLER_SOURCES} "example7.c")
    set(N_LIBS
        ${N_LIBS}
        gcov)
endif()

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

configure_file(eop19620101.txt eop19620101.txt COPYONLY)
configure_file(sat1.txt sat1.txt COPYONLY)
configure_file(sat2.txt sat2.txt COPYONLY)
configure_file(sat3.txt sat3.txt COPYONLY)
configure_file(sat5.txt sat5.txt COPYONLY)
configure_file(sat6.txt sat6.txt COPYONLY)
configure_file(sat7.txt sat7.txt COPYONLY)

add_executable(test-orbit ${TEST_SOURCES})
target_link_libraries(test-orbit ${N_LIBS})

add_executable(example1 ${DOUBLER_SOURCES} "example1.c")
target_link_libraries(example1 ${N_LIBS})

add_executable(example2 ${DOUBLER_SOURCES} "example2.c")
target_link_libraries(example2 ${N_LIBS})

add_executable(example3 ${DOUBLER_SOURCES} "example3.c")
target_link_libraries(example3 ${N_LIBS})

add_executable(example5 ${GAUSS_SOURCES} "example5.c")
target_link_libraries(example5 ${N_LIBS})

add_executable(example6 ${DOUBLER_SOURCES} "example6.c")
target_link_libraries(example6 ${N_LIBS})

add_executable(example7 ${DOUBLER_SOURCES} "example7.c")
target_link_libraries(example7 ${N_LIBS})

