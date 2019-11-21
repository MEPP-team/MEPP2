cmake_minimum_required(VERSION 2.8.11)
project( MEPP2 )

### Building configuration OPTIONs
option( BUILD_DOCUMENTATION                  "Build the documentation (MEPP2)." OFF )
option( BUILD_USE_CGAL                       "CGAL usage toggle." ON )
option( BUILD_USE_OPENMESH                   "OpenMesh usage toggle." ON )
option( BUILD_USE_AIF                        "Adjacency and Incidence Framework toggle." OFF )
option( BUILD_USE_IMG-3RDPARTY               "Img-3rdparty (jpeg, png, tiff) usage toggle." ON )
option( BUILD_USE_GUI                        "Using User Interface." ON )
option( BUILD_USE_PCL                        "PCL usage toggle." OFF )
option( BUILD_USE_VTK                        "VTK usage toggle." OFF )
option( BUILD_USE_FBX                        "FBX usage toggle." OFF )
option( BUILD_USE_DRACO                      "Draco usage toggle." OFF )
option( BUILD_EXAMPLES                       "Build the examples (MEPP2)." ON )
option( BUILD_TESTS                          "Build the tests (MEPP2)." ON )
option( BUILD_ADVANCED_TESTS_COMPR_VALENCE   "Build advanced compression-decompression valence tests (MEPP2)." OFF )

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

# Make sure the linker can find the library once it is built.
link_directories (${PROJECT_BINARY_DIR})

# Set location of compiled libraries
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/Libs)

### Tests
enable_testing()

### Recurse
set( TESTING_DATA_DIR "${PROJECT_SOURCE_DIR}/Testing/Data" )
set( TESTING_REFOUTPUT_DIR "${PROJECT_SOURCE_DIR}/Testing/RefOutput" )

add_subdirectory( Testing )
add_subdirectory( Examples )
add_subdirectory( Visualization )
add_subdirectory( Doc )


# create UseMEPP2.cmake file

function(init_use_file  use_file_name)
  file( WRITE "${use_file_name}"
        "#\n"
        "# This file can be included in a project to set the needed compiler and linker\n"
        "# settings to use MEPP2 in a program.\n"
        "#\n"
        )
  set( USE_MEPP2_FILE "${use_file_name}" PARENT_SCOPE )
endFunction()

function(append_var_to_use_file  variable_name  variable_value)
  file( APPEND "${USE_MEPP2_FILE}"
        "\nset( ${variable_name} \"${variable_value}\" )\n"
        )
endFunction()

init_use_file( "${PROJECT_BINARY_DIR}/UseMEPP2.cmake" )

append_var_to_use_file( "MEPP2_DIR" "${PROJECT_SOURCE_DIR}" )
append_var_to_use_file( "MEPP2_BUILD_DIR" "${PROJECT_BINARY_DIR}" )
append_var_to_use_file( "CMAKE_BUILD_TYPE" "${CMAKE_BUILD_TYPE}" )

# to know which DS is present and
# to know if plugins can be built (need GUI)
append_var_to_use_file( "BUILD_USE_CGAL" "${BUILD_USE_CGAL}" )
append_var_to_use_file( "BUILD_USE_OPENMESH" "${BUILD_USE_OPENMESH}" )
append_var_to_use_file( "BUILD_USE_AIF" "${BUILD_USE_AIF}" )
append_var_to_use_file( "BUILD_USE_GUI" "${BUILD_USE_GUI}" )
append_var_to_use_file( "BUILD_USE_QT5" "${BUILD_USE_QT5}" )

append_var_to_use_file( "CMAKE_CXX_FLAGS" "${CMAKE_CXX_FLAGS}" )
append_var_to_use_file( "CMAKE_CXX_FLAGS_DEBUG" "${CMAKE_CXX_FLAGS_DEBUG}" )
append_var_to_use_file( "CMAKE_CXX_FLAGS_MINSIZEREL"
                        "${CMAKE_CXX_FLAGS_MINSIZEREL}" )
append_var_to_use_file( "CMAKE_CXX_FLAGS_RELEASE"
                        "${CMAKE_CXX_FLAGS_RELEASE}" )
append_var_to_use_file( "CMAKE_CXX_FLAGS_RELWITHDEBINFO"
                        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
append_var_to_use_file( "CMAKE_EXE_LINKER_FLAGS"
                        "${CMAKE_EXE_LINKER_FLAGS}" )
append_var_to_use_file( "CMAKE_MODULE_LINKER_FLAGS"
                        "${CMAKE_MODULE_LINKER_FLAGS}" )
append_var_to_use_file( "CMAKE_SHARED_LINKER_FLAGS"
                        "${CMAKE_SHARED_LINKER_FLAGS}" )
append_var_to_use_file( "CMAKE_STATIC_LINKER_FLAGS"
                        "${CMAKE_STATIC_LINKER_FLAGS}" )
append_var_to_use_file( "CMAKE_EXE_LINKER_FLAGS_RELEASE"
                        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}" )
append_var_to_use_file( "CMAKE_MODULE_LINKER_FLAGS_RELEASE"
                        "${CMAKE_MODULE_LINKER_FLAGS_RELEASE}" )
append_var_to_use_file( "CMAKE_SHARED_LINKER_FLAGS_RELEASE"
                        "${CMAKE_SHARED_LINKER_FLAGS_RELEASE}" )
append_var_to_use_file( "CMAKE_STATIC_LINKER_FLAGS_RELEASE"
                        "${CMAKE_STATIC_LINKER_FLAGS_RELEASE}" )
#!!! see Visu/App/CMakeLists.txt for other CMAKE_.*_FLAGS changed for GUI apps

append_var_to_use_file( "MEPP2_CGAL_INCLUDE_DIRECTORIES"
                        "${MEPP2_CGAL_INCLUDE_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_CGAL_LIBRARIES"
                        "${MEPP2_CGAL_LIBRARIES}" )
append_var_to_use_file( "MEPP2_CGAL_LINK_LIBRARIES"
                        "${MEPP2_CGAL_LINK_LIBRARIES}" )
append_var_to_use_file( "MEPP2_CGAL_LINK_DIRECTORIES"
                        "${MEPP2_CGAL_LINK_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_CGAL_COMPILE_DEFINITIONS"
                        "${MEPP2_CGAL_COMPILE_DEFINITIONS}" )

append_var_to_use_file( "MEPP2_AIF_INCLUDE_DIRECTORIES"
                        "${MEPP2_AIF_INCLUDE_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_AIF_LIBRARIES"
                        "${MEPP2_AIF_LIBRARIES}" )
append_var_to_use_file( "MEPP2_AIF_LINK_LIBRARIES"
                        "${MEPP2_AIF_LINK_LIBRARIES}" )
append_var_to_use_file( "MEPP2_AIF_LINK_DIRECTORIES"
                        "${MEPP2_AIF_LINK_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_AIF_COMPILE_DEFINITIONS"
                        "${MEPP2_AIF_COMPILE_DEFINITIONS}" )

append_var_to_use_file( "MEPP2_OPENMESH_INCLUDE_DIRECTORIES"
                        "${MEPP2_OPENMESH_INCLUDE_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_OPENMESH_LIBRARIES"
                        "${MEPP2_OPENMESH_LIBRARIES}" )
append_var_to_use_file( "MEPP2_OPENMESH_LINK_LIBRARIES"
                        "${MEPP2_OPENMESH_LINK_LIBRARIES}" )
append_var_to_use_file( "MEPP2_OPENMESH_LINK_DIRECTORIES"
                        "${MEPP2_OPENMESH_LINK_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_OPENMESH_COMPILE_DEFINITIONS"
                        "${MEPP2_OPENMESH_COMPILE_DEFINITIONS}" )

append_var_to_use_file( "MEPP2_GUI_INCLUDE_DIRECTORIES"
                        "${MEPP2_GUI_INCLUDE_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_GUI_LIBRARIES"
                        "${MEPP2_GUI_LIBRARIES}" )
append_var_to_use_file( "MEPP2_GUI_LINK_LIBRARIES"
                        "${MEPP2_GUI_LINK_LIBRARIES}" )
append_var_to_use_file( "MEPP2_GUI_LINK_DIRECTORIES"
                        "${MEPP2_GUI_LINK_DIRECTORIES}" )
append_var_to_use_file( "MEPP2_GUI_COMPILE_DEFINITIONS"
                        "${MEPP2_GUI_COMPILE_DEFINITIONS}" )

append_var_to_use_file( "MEPP2_osgQt_SRC"
                        "${MEPP2_osgQt_SRC}" )

message("++ ${USE_MEPP2_FILE} file created.")
