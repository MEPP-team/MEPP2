if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_AIF )
  return()
endif()

if( WIN32 )
  # AUTOLINK is used for CGAL on Windows. Disable it for AIF tests in
  # order to prevent false link error.
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
endif()

# Refer to DevelopersDesignCMakeFiles.dox entry 001
if( BUILD_USE_CGAL )
  include(${CGAL_USE_FILE})
else()
  include_directories( BEFORE "${PROJECT_SOURCE_DIR}/External/CGAL-4.14" )
endif()

# Refer to Doc/Devel/CMakeFiles.md entry 002
add_definitions( -DCGAL_NDEBUG )

#####################################################################

######
add_executable( helloworld_filter_aif
                helloworld_filter_aif.cpp )
target_link_libraries( helloworld_filter_aif
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
