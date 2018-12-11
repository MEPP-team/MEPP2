if( NOT BUILD_TESTS )
  return()
endif()

# Disabling the tests whose names are found in CTestCustom.cmake
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake
               ${CMAKE_BINARY_DIR})

set( TESTING_DIR ${PROJECT_SOURCE_DIR}/Testing )

#------ TESTS subdirectories ------
add_subdirectory( AIF )
add_subdirectory( CGAL )
add_subdirectory( CImg )
add_subdirectory( Draco )
add_subdirectory( OpenMesh )
add_subdirectory( PCL )
add_subdirectory( Utils )
