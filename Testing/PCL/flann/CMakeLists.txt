if( NOT BUILD_USE_PCL OR NOT FLANN_FOUND )
  return()
endif()

include_directories( ${FLANN_INCLUDE_DIRS} )
add_definitions( ${FLANN_DEFINITIONS} )
add_executable(flann_example flann_example.cpp)
target_link_libraries(flann_example ${FLANN_LIBRARIES})

add_test(flann_example flann_example ${TESTING_DATA_DIR}/flann/dataset.dat ${TESTING_DATA_DIR}/flann/testset.dat)
