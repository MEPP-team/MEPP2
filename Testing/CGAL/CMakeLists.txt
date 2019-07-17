if( NOT BUILD_USE_CGAL )
  return()
endif()

include(${CGAL_USE_FILE})

#####################################################################

add_subdirectory( LCC )
add_subdirectory( Point_set )
add_subdirectory( Polyhedron )
add_subdirectory( Surface_mesh )

#####################################################################

## variables to be exported to external projects

# include dirs
get_property( INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY INCLUDE_DIRECTORIES
              )
set( MEPP2_CGAL_INCLUDE_DIRECTORIES
     "${INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# libraries to link with
set( MEPP2_CGAL_LIBRARIES
     ${CGAL_LIBRARY}
     ${CGAL_3RD_PARTY_LIBRARIES}
     ${Boost_LIBRARIES}
     ${FBX_LIBRARY}
     ${VTK_LIBRARIES}
     CACHE INTERNAL "no-desc" FORCE
     )

# link libraries (to link all target with, silently!)
get_property( LINK_LIBRARIES_TMP
              DIRECTORY
              PROPERTY LINK_LIBRARIES
              )
set( MEPP2_CGAL_LINK_LIBRARIES
     "${LINK_LIBRARIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# link dirs
get_property( LINK_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY LINK_DIRECTORIES
              )
set( MEPP2_CGAL_LINK_DIRECTORIES
     "${LINK_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# compilation flags and define
get_property( COMPILE_DEFINITIONS_TMP
              DIRECTORY
              PROPERTY COMPILE_DEFINITIONS
              )
set( MEPP2_CGAL_COMPILE_DEFINITIONS
     ${COMPILE_DEFINITIONS_TMP}
     CACHE INTERNAL "no-desc" FORCE
     )