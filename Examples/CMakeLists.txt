if( NOT BUILD_EXAMPLES )
  return()
endif()

#------ Examples subdirectories ------
add_subdirectory( OpenMesh )
add_subdirectory( CGAL )
add_subdirectory( AIF )
add_subdirectory( PCL )
