if( NOT BUILD_USE_IMG-3RDPARTY )
  return()
endif()

include_directories( ${PROJECT_SOURCE_DIR}/External/CImg-2.4.1 )

if ( JPEG_FOUND )
include_directories( ${JPEG_INCLUDE_DIR} )
endif ()

if ( PNG_FOUND )
include_directories( ${PNG_PNG_INCLUDE_DIR} )
endif ()

if ( TIFF_FOUND )
include_directories( ${TIFF_INCLUDE_DIR} )
endif ()

add_executable(cimg_tests cimg_tests.cpp)

if ( JPEG_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${JPEG_LIBRARIES} )
endif ()

if ( PNG_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${PNG_LIBRARIES} )
endif ()

if ( TIFF_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${TIFF_LIBRARIES} )
endif ()

target_link_libraries(cimg_tests ${IMG-3RDPARTY_LIBRARIES})

add_test(cimg_tests cimg_tests ${TESTING_DATA_DIR}/cimg/cnrs.jpg ${TESTING_DATA_DIR}/cimg/liris.png ${TESTING_DATA_DIR}/cimg/mepp.tif)
