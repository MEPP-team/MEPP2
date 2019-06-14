if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_CGAL )
  return()
endif()

include(${CGAL_USE_FILE})

#####################################################################

add_subdirectory( LCC )
add_subdirectory( Point_set )
add_subdirectory( Polyhedron )
add_subdirectory( Surface_mesh )
