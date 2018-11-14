if( NOT BUILD_USE_IMG-3RDPARTY )
  return()
endif()

# ---

include_directories( ${PROJECT_SOURCE_DIR}/External/CImg/CImg-2.4.2_pre111418 )

if ( DEFINED ENV{TRAVIS} OR DEFINED ENV{APPVEYOR} )
add_definitions( -Dcimg_display=0 ) # Flags used to disable display capablities of CImg (for all CIs)
endif ()

# ---

if ( JPEG_FOUND )
include_directories( ${JPEG_INCLUDE_DIR} )
#message( STATUS "-----> JPEG_INCLUDE_DIR: ${JPEG_INCLUDE_DIR}" )
endif ()

if ( PNG_FOUND )
include_directories( ${PNG_PNG_INCLUDE_DIR} )
#include_directories( ${PNG_INCLUDE_DIR} )
#message( STATUS "-----> PNG_PNG_INCLUDE_DIR: ${PNG_PNG_INCLUDE_DIR}" )
#message( STATUS "-------> PNG_INCLUDE_DIR: ${PNG_INCLUDE_DIR}" )
endif ()

if ( TIFF_FOUND )
include_directories( ${TIFF_INCLUDE_DIR} )
#message( STATUS "-----> TIFF_INCLUDE_DIR: ${TIFF_INCLUDE_DIR}" )
endif ()

# ---

if ( JPEG_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${JPEG_LIBRARIES} )
#message( STATUS "-----> JPEG_LIBRARIES: ${JPEG_LIBRARIES}" )
endif ()

if ( PNG_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${PNG_LIBRARIES} )
#message( STATUS "-----> PNG_LIBRARIES: ${PNG_LIBRARIES}" )
endif ()

if ( TIFF_FOUND )
set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${TIFF_LIBRARIES} )
#message( STATUS "-----> TIFF_LIBRARIES: ${TIFF_LIBRARIES}" )
endif ()

# ---

add_executable(cimg_tests cimg_tests.cpp)
target_link_libraries(cimg_tests ${IMG-3RDPARTY_LIBRARIES})
add_test(cimg_tests cimg_tests ${TESTING_DATA_DIR}/cimg/cnrs.jpg ${TESTING_DATA_DIR}/cimg/liris.png ${TESTING_DATA_DIR}/cimg/mepp.tif)
