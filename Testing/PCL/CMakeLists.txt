if( NOT BUILD_USE_PCL )
  return()
endif()

#link_directories(${PCL_LIBRARY_DIRS})
include_directories( ${PCL_INCLUDE_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )

if( MSVC )
  # AUTOLINK is used for CGAL on Windows. Disable it for OpenMesh-only tests
  # in order to prevent false link error.
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

#####
add_executable(pcd_write_test pcd_write.cpp)
target_link_libraries(pcd_write_test ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${Boost_LIBRARIES})
add_test(pcd_write_test pcd_write_test)

#####
add_executable( test_complying_concepts_pcl
                test_complying_concepts_pcl.cpp )
target_link_libraries( test_complying_concepts_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${Boost_LIBRARIES}
                       )
add_test( Test_complying_concepts_pcl
          test_complying_concepts_pcl )

#####
add_executable( test_property_maps_pcl
                test_property_maps_pcl.cpp )
target_link_libraries( test_property_maps_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${Boost_LIBRARIES}
                       )
add_test( Test_property_maps_pcl
          test_property_maps_pcl )

#####
add_executable( test_reader_writer_pcl
                test_reader_writer_pcl.cpp )
target_link_libraries( test_reader_writer_pcl
                       ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${Boost_LIBRARIES}
                       ${VTK_LIBRARIES} ${FBX_LIBRARY}
                       )
add_test( Test_reader_writer_pcl
          test_reader_writer_pcl
          ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
          test_reader_writer_pcl.ply
          ${TESTING_REFOUTPUT_DIR}/PCL/test_reader_writer_pcl.ref.ply
          )


#------

add_subdirectory( flann )
