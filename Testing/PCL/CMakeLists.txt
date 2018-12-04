if( NOT BUILD_USE_PCL )
  return()
endif()

#link_directories(${PCL_LIBRARY_DIRS})

include_directories( ${PCL_INCLUDE_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )
add_executable(pcd_write_test pcd_write.cpp)
target_link_libraries(pcd_write_test ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${Boost_LIBRARIES})

add_test(pcd_write_test pcd_write_test)

#------

add_subdirectory( flann )
