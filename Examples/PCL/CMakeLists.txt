if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_PCL )
  return()
endif()

#####################################################################

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

link_directories(${PCL_LIBRARY_DIRS})
include_directories( ${PCL_INCLUDE_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )

#####################################################################

#####
add_executable( compute_normals_pcl
                compute_normals_pcl.cpp )
target_link_libraries( compute_normals_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES}
                       ${PCL_KDTREE} ${PCL_SEARCH_LIBRARIES} ${PCL_FEATURES_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${FBX_LIBRARY} 
                       ${FLANN_LIBRARIES}
                       )

#####
add_executable( shot_color_estimation_pcl
                shot_color_estimation_pcl.cpp )
target_link_libraries( shot_color_estimation_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES}
                       ${PCL_KDTREE} ${PCL_SEARCH_LIBRARIES} ${PCL_FEATURES_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${FBX_LIBRARY}
                       ${FLANN_LIBRARIES}
                       )

#####
add_executable( example_prop_maps_pcl
                example_prop_maps_pcl.cpp )
target_link_libraries( example_prop_maps_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES}
                       ${PCL_KDTREE} ${PCL_SEARCH_LIBRARIES} ${PCL_FEATURES_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${FBX_LIBRARY} )