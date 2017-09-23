# - Config file for the Armadillo package
# It defines the following variables
#  ARMADILLO_INCLUDE_DIRS - include directories for Armadillo
#  ARMADILLO_LIBRARY_DIRS - library directories for Armadillo (normally not used!)
#  ARMADILLO_LIBRARIES    - libraries to link against

# Tell the user project where to find our headers and libraries
set(ARMADILLO_INCLUDE_DIRS "/Users/Gleb/Desktop/Triangular Solvers/lib_source/armadillo-8.100.0/tmp/include")
set(ARMADILLO_LIBRARY_DIRS "/Users/Gleb/Desktop/Triangular Solvers/lib_source/armadillo-8.100.0")

# Our library dependencies (contains definitions for IMPORTED targets)
include("/Users/Gleb/Desktop/Triangular Solvers/lib_source/armadillo-8.100.0/ArmadilloLibraryDepends.cmake")

# These are IMPORTED targets created by ArmadilloLibraryDepends.cmake
set(ARMADILLO_LIBRARIES armadillo)
