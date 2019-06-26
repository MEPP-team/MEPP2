if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_PCL )
  return()
endif()

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
                       ${Boost_LIBRARIES} )
