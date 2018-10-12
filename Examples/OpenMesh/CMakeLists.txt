if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_OPENMESH )
  return()
endif()

include_directories( ${OPENMESH_INCLUDE_DIR} )

if( MSVC )
  # AUTOLINK is used for CGAL on Windows. Disable it for OpenMesh in order
  # to prevent false link error.
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
endif()

# Refer to DevelopersDesignCMakeFiles.dox entry 001
if( BUILD_USE_CGAL )
  include(${CGAL_USE_FILE})
else()
  include_directories( BEFORE "${PROJECT_SOURCE_DIR}/External/CGAL-4.11" )
endif()

# Refer to Doc/Devel/CMakeFiles.md entry 002
add_definitions( -DCGAL_NDEBUG )

#####################################################################

#####
add_executable( example_accessing_geometry
                accessing_geometry.cpp )
target_link_libraries( example_accessing_geometry
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )

add_subdirectory( HelloWorld )
