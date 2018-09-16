if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_AIF )
  return()
endif()

if( WIN32 )
  # AUTOLINK is used for CGAL on Windows. Disable it for AIF tests in
  # order to prevent false link error.
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
endif()

# Refer to Doc/Devel/CMakeFiles.md entry 001
include_directories( BEFORE "${PROJECT_SOURCE_DIR}/External/CGAL-4.11" )

# Refer to Doc/Devel/CMakeFiles.md entry 002
add_definitions( -DCGAL_NDEBUG )
#####################################################################

add_subdirectory( HelloWorld )
