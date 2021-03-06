cmake_minimum_required(VERSION 3.9)
project(Triangular_Solvers)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMAKE_MODULES")

IF(APPLE)

    set(CMAKE_C_COMPILER /usr/local/bin/clang-omp)
    set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++)

    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    SET(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
    set(OpenMP_EXE_LINKER_LFLAGS "-lomp")

#    FIND_PACKAGE(OpenMP REQUIRED)
#    if(OPENMP_FOUND)
#        message("OPENMP FOUND")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#        include_directories(${OpenMP_INCLUDE_DIRS})
 #   endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF(FADE2D_LIBRARY)

    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)

    include_directories(${FADE2D_INCLUDE_DIRS})

    FIND_PACKAGE(Armadillo REQUIRED)
    INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

ELSEIF(UNIX)



    set(CMAKE_CXX_COMPILER "clang")
    set(CMAKE_CXX_COMPILER "clang++")

    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
    endif()

    FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_ubuntu)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF(FADE2D_LIBRARY)

    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)
    include_directories(${FADE2D_INCLUDE_DIRS})

    FIND_PACKAGE(Armadillo REQUIRED)
    INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

    FIND_PACKAGE(MKL REQUIRED)


ELSE()


    FIND_LIBRARY(FADE2D_LIBRARY fade2D_x64_v141_Release ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_win64)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF()


    #FIND_LIBRARY(GMP_LIBRARY libgmp-10 ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_win64)


    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)
    include_directories(${FADE2D_INCLUDE_DIRS})
	set(ARMA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/armadillo-8.200.0/include)
	include_directories(${ARMA_INCLUDE_DIRS})


#    FIND_LIBRARY(LAPACK_LIBRARY liblapack ${CMAKE_SOURCE_DIR}/lib_source/lib/lapack_win64)
#    IF(LAPACK_LIBRARY)
#        message("lapack")
#    ENDIF()

#    FIND_LIBRARY(BLAS_LIBRARY libblas ${CMAKE_SOURCE_DIR}/lib_source/lib/lapack_win64)
#    IF(LAPACK_LIBRARY)
#        message("lapack")
#    ENDIF()

    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
    endif()



ENDIF()




set(SOURCE_FILES main.cpp SpaceMesh/Area.cpp SpaceMesh/Area.h
                 SpaceMesh/TriangularMesh.cpp SpaceMesh/TriangularMesh.h
                 FirstOrderSolver/Solver.h FirstOrderSolver/Solver.cpp Maths/Algebra.cpp
                 Maths/Algebra.h FirstOrderSolver/LaxFriedrichSolver.cpp FirstOrderSolver/LaxFriedrichSolver.h
                 WENO/WENOSolver.cpp WENO/WENOSolver.h FirstOrderSolver/RoeSolver.h FirstOrderSolver/RoeSolver.cpp HYBRID/HYBRIDSolver.cpp HYBRID/HYBRIDSolver.h)
add_executable(Triangular_Solvers ${SOURCE_FILES})
IF(APPLE)
    target_link_libraries(Triangular_Solvers ${OpenMP_omp_LIBRARY} ${FADE2D_LIBRARY} ${ARMADILLO_LIBRARIES})
ELSEIF(UNIX)
    target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY} ${ARMADILLO_LIBRARIES} ${MKL_LIBRARIES})
ELSE()
    #target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY} ${LAPACK_LIBRARY} ${GMP_LIBRARY} ${BLASS_LIBRARY})
	target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY})
ENDIF()
