cmake_minimum_required(VERSION 3.7)
project(Triangular_Solvers)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE FILEPATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE FILEPATH "C++ compiler" FORCE)

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    include_directories(${OpenMP_INCLUDE_DIRS})
endif()

FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib)
IF(FADE2D_LIBRARY)
    message("kek")
ENDIF(FADE2D_LIBRARY)

set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)

include_directories(${FADE2D_INCLUDE_DIRS})

#FIND_LIBRARY(ARMA_LIBRARY armadillo.8.10.0 ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib)
#IF(ARMA_LIBRARY)
#    message("kekasik")
#ENDIF(ARMA_LIBRARY)

#set(ARMA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/armadillo-8.100.0/include)

#include_directories(${ARMA_INCLUDE_DIRS})

FIND_PACKAGE(Armadillo REQUIRED)
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_BUILD_TYPE Debug)


set(SOURCE_FILES main.cpp SpaceMesh/Area.cpp SpaceMesh/Area.h SpaceMesh/TriangularMesh.cpp SpaceMesh/TriangularMesh.h FirstOrderSolver/Solver.h FirstOrderSolver/Solver.cpp Maths/Algebra.cpp Maths/Algebra.h FirstOrderSolver/LaxFriedrichSolver.cpp FirstOrderSolver/LaxFriedrichSolver.h WENO/WENOLF.cpp WENO/WENOLF.h)
add_executable(Triangular_Solvers ${SOURCE_FILES})
target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY} ${ARMADILLO_LIBRARIES})