cmake_minimum_required(VERSION 2.8.11)
project( MEPP2 )

### Building configuration OPTIONs
option( BUILD_DOCUMENTATION  "Build the documentation." OFF )
option( BUILD_USE_AIF        "Adjacency and Incidence Framework toggle." OFF)
option( BUILD_USE_CGAL       "CGAL usage toggle." ON )
option( BUILD_USE_GUI        "Using User Interface." OFF)
option( BUILD_USE_OPENMESH   "OpenMesh usage toggle." ON )
option( BUILD_USE_PCL        "PCL usage toggle." OFF )
option( BUILD_USE_VTK        "VTK usage toggle." OFF )
option( BUILD_USE_FBX        "FBX usage toggle." OFF )
option( BUILD_EXAMPLES       "Build the examples." OFF )
option( BUILD_TESTS          "Build the tests." ON )
option( BUILD_ADVANCED_TESTS_COMPR_VALENCE "Build advanced compression-decompression valence tests." OFF )

if( BUILD_USE_AIF )
	option( BUILD_USE_GENERATE_TRIANGLE		"2D (and soon 3D) triangle meshes generator." ON )
endif( BUILD_USE_AIF )

#####
# Set CMAKE module path unless provided by user
if( NOT CMAKE_MODULE_PATH )
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Cmake" )
  message( STATUS "Cmake module path set to " ${PROJECT_SOURCE_DIR}/Cmake )
endif()

# Set CMake Policies
include( "${PROJECT_SOURCE_DIR}/Cmake/CMakePolicies.cmake" )

# Set default build type
if( NOT CMAKE_BUILD_TYPE )
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Set some compiler options
include( "${PROJECT_SOURCE_DIR}/Cmake/CompilerOptions.cmake" )

# Set OS options
include( "${PROJECT_SOURCE_DIR}/Cmake/OSOptions.cmake" )

# Find Dependencies
include( "${PROJECT_SOURCE_DIR}/Cmake/FindDependencies.cmake" )

#####
# settings and libs common to Testing and Example

# Make sure the compiler can find include files.
include_directories (${PROJECT_SOURCE_DIR}/)
include_directories (${PROJECT_BINARY_DIR}/)

# Make sure the linker can find the library once it is built.
link_directories (${PROJECT_BINARY_DIR})

# Set location of compiled libraries
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/Libs)

### Tests
enable_testing()

### Recurse
set( TESTING_DATA_DIR "${PROJECT_SOURCE_DIR}/Testing/Data" )

add_subdirectory( Testing )
add_subdirectory( Examples )
add_subdirectory( Visualization )
add_subdirectory( Doc )
