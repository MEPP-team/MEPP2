if( NOT BUILD_USE_IMG-3RDPARTY )
  return()
endif()

include_directories( "${PROJECT_SOURCE_DIR}/External/CImg/CImg-2.4.2" )

add_executable(cimg_tests cimg_tests.cpp)
target_link_libraries(cimg_tests ${IMG-3RDPARTY_LIBRARIES})

add_test(cimg_tests cimg_tests ${TESTING_DATA_DIR}/cimg/cnrs.jpg ${TESTING_DATA_DIR}/cimg/liris.png ${TESTING_DATA_DIR}/cimg/mepp.tif)
