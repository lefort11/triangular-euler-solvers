cmake_minimum_required(VERSION 3.7)
project(Triangular_Solvers)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE FILEPATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE FILEPATH "C++ compiler" FORCE)


FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib)
IF(FADE2D_LIBRARY)
    message("kek")
ENDIF(FADE2D_LIBRARY)

set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)

include_directories(${FADE2D_INCLUDE_DIRS})

set(CMAKE_BUILD_TYPE Debug)


set(SOURCE_FILES main.cpp SpaceMesh/Area.cpp SpaceMesh/Area.h SpaceMesh/TriangularMesh.cpp SpaceMesh/TriangularMesh.h FirstOrderSolver/Solver.h FirstOrderSolver/Solver.cpp Maths/Algebra.cpp Maths/Algebra.h)
add_executable(Triangular_Solvers ${SOURCE_FILES})
target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY})